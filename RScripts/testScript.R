# rooftop_solar_tile.R
# Purpose: hourly rooftop-focused solar potential for a 10km DSM tile (EPSG:27700)


library(terra)
library(sf)
library(rgrass)

library(reticulate)
library(data.table)

# Download GRASS from https://grass.osgeo.org/download/windows/

# -------------------------
# User inputs / file paths
# -------------------------
dsm_tif      <- "sampleData/DSM/SE23.tiff"     # 2 m DSM, EPSG:27700
buildings_shp <- "sampleData/Building_heights/SE23.Rds"  # building polygons (optional but recommended)
out_dir      <- "outputs_tile"
dir.create(out_dir, showWarnings = FALSE)

# Weather data (you must provide one of these)
era5_dir     <- "sampleData/ERA5/byGrid/"   
era5 = readRDS(file.path(era5_dir,"SE23.Rds"))

# Simulation period (typical year)
start_time <- as.POSIXct("2020-01-01 00:00", tz="UTC")
end_time   <- as.POSIXct("2020-12-31 23:00", tz="UTC")
times      <- seq(start_time, end_time, by = "1 hour")   # hourly (8760)

# Default tilt for non-building pixels (degrees)
default_tilt_deg <- 30

# PV system defaults (simple)
pv_module_tilt_source <- "per_pixel"  # "per_pixel" or "fixed"
pv_module_azimuth_source <- "per_pixel"
system_losses_fraction <- 0.14  # 14% total system losses (soiling, wiring, inverter, etc.)

# -------------------------
# 1. Load DSM and buildings
# -------------------------
dsm <- rast(dsm_tif)
#crs(dsm) <- "EPSG:27700"

buildings <- NULL
if (file.exists(buildings_shp)) {
  buildings <- readRDS(buildings_shp)
}

# -------------------------
# 2. Derive slope (tilt) and aspect (azimuth) rasters
# -------------------------
# Use terra::terrain to compute slope/aspect in degrees
slope_r <- terrain(dsm, v = "slope", unit = "degrees", neighbors = 8)
aspect_r <- terrain(dsm, v = "aspect", unit = "degrees", neighbors = 8)

# For flat/no-data areas, set defaults later
# Save intermediate rasters
writeRaster(slope_r, filename = file.path(out_dir, "slope_deg.tif"), overwrite = TRUE)
writeRaster(aspect_r, filename = file.path(out_dir, "aspect_deg.tif"), overwrite = TRUE)

# -------------------------
# 3. Create roof mask and per-pixel tilt/azimuth
# -------------------------
# Rasterize building polygons to DSM resolution to get roof mask
if (!is.null(buildings)) {
  b_r <- rasterize(vect(buildings), dsm, field = 1, background = NA)
  writeRaster(b_r, filename = file.path(out_dir, "roof_mask.tif"), overwrite = TRUE)
} else {
  b_r <- rast(dsm) * NA
}

# Per-pixel tilt: use slope where building mask == 1, else default tilt
tilt_r <- slope_r
tilt_r[is.na(tilt_r)] <- default_tilt_deg
if (!is.null(buildings)) {
  # For non-building pixels, set default tilt
  nonroof_idx <- is.na(b_r)
  tilt_r[nonroof_idx] <- default_tilt_deg
}

# Per-pixel azimuth: use aspect for roofs; for non-roofs set to south (180)
azimuth_r <- aspect_r
azimuth_r[is.na(azimuth_r)] <- 180
if (!is.null(buildings)) {
  azimuth_r[is.na(b_r)] <- 180
}

writeRaster(tilt_r, filename = file.path(out_dir, "tilt_deg.tif"), overwrite = TRUE)
writeRaster(azimuth_r, filename = file.path(out_dir, "azimuth_deg.tif"), overwrite = TRUE)

# -------------------------
# 4. Initialize GRASS and import DSM
# -------------------------
# Adjust gisBase to your GRASS installation path
gisBase <- "C:/Program Files/GRASS GIS 8.4"   # <-- change to your GRASS path


loc = initGRASS(gisBase = gisBase, home = tempdir(), 
                gisDbase = file.path(tempdir(),"grassdb"),
                mapset = "PERMANENT",
                override = TRUE)

# Create location with EPSG:27700
execGRASS("g.proj", flags = "c", epsg = 27700)

# Export DSM to a temporary GeoTIFF and import into GRASS
tmp_dsm <- file.path(tempdir(), "dsm_for_grass.tif")
writeRaster(dsm, tmp_dsm, overwrite = TRUE)
execGRASS("r.in.gdal", flags = c("o","overwrite"), input = tmp_dsm, output = "dsm")

# Import slope/aspect if desired
tmp_slope <- file.path(tempdir(), "slope_for_grass.tif")
writeRaster(slope_r, tmp_slope, overwrite = TRUE)
execGRASS("r.in.gdal", flags = c("o","overwrite"), input = tmp_slope, output = "slope")

tmp_aspect <- file.path(tempdir(), "aspect_for_grass.tif")
writeRaster(aspect_r, tmp_aspect, overwrite = TRUE)
execGRASS("r.in.gdal", flags = c("o","overwrite"), input = tmp_aspect, output = "aspect")


# -------------------------
# 7. Run r.sun hourly to compute beam/diffuse/global on horizontal surface
# -------------------------
# r.sun expects day-of-year and time (decimal hours). We'll loop hourly.
# Note: r.sun can compute beam/diffuse/reflected using horizon raster.
# For speed, run only daylight hours per day; here we run full 0-23 for simplicity.

# Return to 2m resolution
execGRASS("g.region", raster = "dsm", flags = "p") # back to fine region

# Create a vector of day-of-year and hour
doy <- as.integer(format(times, "%j"))
hour_decimal <- as.numeric(format(times, "%H")) + as.numeric(format(times, "%M"))/60

# We'll store hourly global horizontal irradiance (GHI) raster outputs in a temp mapset
# For large runs, parallelize by day or tile and write to disk incrementally.
for (i in seq_along(times)) {
  t <- times[i]
  day_i <- doy[i]
  hour_i <- hour_decimal[i]
  out_prefix <- paste0("sun_", format(t, "%Y%m%d%H"))
  # r.sun parameters: day, time, horizon, beam_rad, diff_rad, glob_rad
  execGRASS("r.sun",
            flags = c("overwrite"),
            parameters = list(elevation = "dsm",
                              aspect = "aspect",
                              slope = "slope",
                              #horizon_basename = "horizon",
                              #horizon_step = 1,
                              day = day_i,
                              time = hour_i,
                              nprocs = 1, # To enable multicore processing
                              beam_rad = paste0(out_prefix, "_beam"),
                              diff_rad = paste0(out_prefix, "_diff"),
                              glob_rad = paste0(out_prefix, "_glob")
                             ))
  
  # Export glob to GeoTIFF for later transposition
  out_tif <- file.path(out_dir, paste0("ghi_", format(t, "%Y%m%d%H"), ".tif"))
  execGRASS("r.out.gdal", flags = c("overwrite"), input = paste0(out_prefix, "_glob"),
            output = out_tif, format = "GTiff", createopt = "COMPRESS=LZW",
            type = "Float32",  nodata = -9999)
  # Optional: remove GRASS maps to save memory (execGRASS("g.remove", ...))
}

# Process ERA5 Data
era5 = era5[,c("timestamp","SSRC","SSRD","2T","SP","TCC","10V","10U")]
names(era5) = c("timestamp",
                "surface_net_solar_radiation_clear_sky",
                "surface_solar_radiation_downwards",
                "2m_temperature",
                "surface_pressure",
                "total_cloud_cover",
                "10m_u_component_of_wind",
                "10m_v_component_of_wind"
                )
# undefined [-]                                     var22
# Surface net solar radiation clear sky [W/m^2]     SSRC
# undefined [-]                                     var129 
# Surface solar radiation downwards [W*s/m^2]       SSRD 
# undefined [-]                                     var21
# 2 metre temperature [C]                           2T
# Surface pressure [Pa]                             SP 
# Total cloud cover (0 - 1) [-]                     TCC 
# 10 metre v wind component [m/s]                   10V 
# 10 metre u wind component [m/s]                   10U  


