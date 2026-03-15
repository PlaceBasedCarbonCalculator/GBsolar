# rooftop_solar_tile.R
# Purpose: hourly rooftop-focused solar potential for a 10km DSM tile (EPSG:27700)
# Requires: terra, sf, rgrass7, RSAGA, reticulate, data.table
# System: GRASS GIS (r.sun, r.horizon), SAGA (optional SVF), Python pvlib

library(terra); library(sf); library(rgrass); library(RSAGA)
library(reticulate); library(data.table)

# -------------------------
# User inputs / file paths
# -------------------------
dsm_tif      <- "sampleData/DSM/SE23.tiff"     # 2 m DSM, EPSG:27700
buildings_shp <- "sampleData/Building_heights/SE23.Rds"  # building polygons (optional but recommended)
out_dir      <- "outputs_tile"
dir.create(out_dir, showWarnings = FALSE)

# Weather data (you must provide one of these)
# Option A: ERA5 hourly gridded fields (downscaled/bias-corrected) covering tile
era5_dir     <- "era5_hourly_tile/"   # folder with hourly netCDF or CSV per variable
# Option B: MIDAS station CSV for bias-correction/validation
midas_csv    <- "midas_station_hourly.csv"

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
gisBase <- "C:/OSGeo4W/lib/grass78"   # <-- change to your GRASS path
initGRASS(gisBase = gisBase, home = tempdir(), gisDbase = file.path(tempdir(),"grassdb"),
          override = TRUE)
# Create location with EPSG:27700
execGRASS("g.proj", flags = "c", proj4 = paste0("+init=epsg:27700"))

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
# 5. Compute horizon / sky-view factor (expensive) and cache
# -------------------------
# r.horizon computes horizon angle for each azimuth step; r.sun can use horizon
# Example: compute horizon with step=1 degree (adjust for speed/accuracy)
execGRASS("r.horizon", flags = c("overwrite"), input = "dsm", output = "horizon", step = 1.0)

# Optionally compute Sky View Factor (SVF) using SAGA via RSAGA (faster alternative)
saga_env <- rsaga.env()
# Convert DSM to SAGA grid first (RSAGA expects file paths)
saga_dem <- file.path(out_dir, "dsm.sgrd")
rsaga.import.gdal(in.file = tmp_dsm, out.grid = "dsm_saga", env = saga_env)
# Compute SVF (SAGA module name may vary by version)
# rsaga.svf(in.dem = "dsm_saga", out.svf = "svf.sgrd", env = saga_env)  # uncomment if available

# -------------------------
# 6. Prepare weather inputs (hourly) and bias correction
# -------------------------
# You must provide hourly GHI/DNI/DHI or at least GHI + clearness model inputs.
# Here we show a placeholder to load ERA5 hourly GHI (surface_solar_radiation_downwards)
# and temperature; if you only have GHI, we will estimate DNI/DHI using a clearness model.

# Placeholder: load ERA5 CSV with columns: time, ghi, dni, dhi, temp_air, wind_speed
# If you only have ghi, compute dni/dhi using a decomposition model (not implemented here).
weather_csv <- file.path(era5_dir, "era5_tile_hourly.csv")
if (!file.exists(weather_csv)) {
  stop("Provide hourly weather CSV (era5_tile_hourly.csv) with columns: time, ghi, dni, dhi, temp_air")
}
weather <- fread(weather_csv)
weather[, time := as.POSIXct(time, tz = "UTC")]
setkey(weather, time)

# Optional: bias-correct ERA5 GHI using MIDAS station observations (if available)
if (file.exists(midas_csv)) {
  midas <- fread(midas_csv)
  midas[, time := as.POSIXct(time, tz = "UTC")]
  # Simple bias correction: compute monthly ratio of station GHI / ERA5 GHI at nearest grid cell
  # (Detailed spatial bias correction is recommended for national runs)
  # ... (user can implement more advanced correction here)
}

# -------------------------
# 7. Run r.sun hourly to compute beam/diffuse/global on horizontal surface
# -------------------------
# r.sun expects day-of-year and time (decimal hours). We'll loop hourly.
# Note: r.sun can compute beam/diffuse/reflected using horizon raster.
# For speed, run only daylight hours per day; here we run full 0-23 for simplicity.

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
            parameters = list(elevin = "dsm",
                              aspin = "aspect",
                              slopein = "slope",
                              horizon_basename = "horizon",
                              day = day_i,
                              time = hour_i,
                              beam_rad = paste0(out_prefix, "_beam"),
                              diff_rad = paste0(out_prefix, "_diff"),
                              glob_rad = paste0(out_prefix, "_glob")))
  # Export glob to GeoTIFF for later transposition
  out_tif <- file.path(out_dir, paste0("ghi_", format(t, "%Y%m%d%H"), ".tif"))
  execGRASS("r.out.gdal", flags = c("overwrite"), input = paste0(out_prefix, "_glob"),
            output = out_tif, format = "GTiff", createopt = "COMPRESS=LZW")
  # Optional: remove GRASS maps to save memory (execGRASS("g.remove", ...))
}

# -------------------------
# 8. Read hourly GHI rasters and prepare inputs for pvlib transposition
# -------------------------
# For each hour, we will:
#  - read GHI raster (from r.sun output)
#  - compute solar position (pvlib) for tile centroid
#  - compute POA using pvlib.irradiance.get_total_irradiance (Perez)
#  - apply system losses and accumulate energy

# Initialize Python pvlib via reticulate
use_python("/usr/bin/python3", required = TRUE)  # adjust path as needed
pvlib <- import("pvlib")
np <- import("numpy")

# Get tile centroid lat/lon (convert from EPSG:27700)
tile_centroid <- as.vector(ext(dsm)[1:2])  # not exact; better compute centroid
centroid_xy <- crds(centroid(vect(dsm)))
# Convert centroid to lat/lon
centroid_sf <- st_as_sf(as.data.frame(t(centroid_xy)), coords = c("V1","V2"), crs = 27700)
centroid_ll <- st_transform(centroid_sf, 4326)
lat <- st_coordinates(centroid_ll)[1,2]
lon <- st_coordinates(centroid_ll)[1,1]

# Prepare output rasters for accumulated annual POA energy (kWh/m2)
accum_poa <- rast(dsm)
values(accum_poa) <- 0

# Loop through times (for large runs, process in chunks)
for (i in seq_along(times)) {
  t <- times[i]
  ghi_tif <- file.path(out_dir, paste0("ghi_", format(t, "%Y%m%d%H"), ".tif"))
  if (!file.exists(ghi_tif)) next
  ghi_r <- rast(ghi_tif)
  ghi_vals <- values(ghi_r)
  # If r.sun produced beam/diff rasters, you can read them similarly and use them
  # For pvlib transposition we need: dni, dhi, ghi. If only ghi available, decompose (not shown).
  # Here we assume weather CSV provided dni/dhi for the tile centroid
  wrow <- weather[time == t]
  if (nrow(wrow) == 0) next
  ghi <- wrow$ghi
  dni <- wrow$dni
  dhi <- wrow$dhi
  temp_air <- wrow$temp_air
  
  # Compute solar position for timestamp at tile centroid
  solarpos <- pvlib$solarposition$get_solarposition(time = as.POSIXct(t, tz="UTC"),
                                                    latitude = lat, longitude = lon)
  solar_zenith <- as.numeric(solarpos$zenith)
  solar_azimuth <- as.numeric(solarpos$azimuth)
  
  # Prepare per-pixel tilt/azimuth arrays (flatten)
  tilt_vals <- values(tilt_r)
  az_vals   <- values(azimuth_r)
  
  # Use pvlib.irradiance.get_total_irradiance (Perez model)
  # pvlib expects scalar or arrays; we'll loop per-pixel in R is slow.
  # Efficient approach: aggregate by unique tilt/azimuth combinations or vectorize in Python.
  # Here we show a simple vectorized call in Python via reticulate:
  get_total_irradiance <- pvlib$irradiance$get_total_irradiance
  # Convert R vectors to numpy arrays
  np_tilt <- np$array(tilt_vals)
  np_az   <- np$array(az_vals)
  # Call get_total_irradiance with arrays (assumes scalar dni/dhi/ghi and scalar solarpos)
  poa_dict <- get_total_irradiance(surface_tilt = np_tilt,
                                   surface_azimuth = np_az,
                                   dni = as.numeric(dni),
                                   ghi = as.numeric(ghi),
                                   dhi = as.numeric(dhi),
                                   solar_zenith = as.numeric(solar_zenith),
                                   solar_azimuth = as.numeric(solar_azimuth),
                                   model = "perez")
  poa_global <- as.numeric(poa_dict["poa_global"])
  # Apply system losses
  poa_after_losses <- poa_global * (1 - system_losses_fraction)
  # Convert W/m2 to Wh for the hour (W/m2 * 1 hour = Wh/m2)
  poa_wh <- poa_after_losses  # since units are W/m2 and timestep is 1 hour
  # Accumulate into raster (kWh/m2)
  current_vals <- values(accum_poa)
  current_vals[!is.na(current_vals)] <- current_vals[!is.na(current_vals)] + (poa_wh[!is.na(current_vals)] / 1000)
  values(accum_poa) <- current_vals
}

# -------------------------
# 9. Save annual results and summary
# -------------------------
writeRaster(accum_poa, filename = file.path(out_dir, "annual_poa_kwh_m2.tif"), overwrite = TRUE)

# Compute summary stats for roofs vs non-roofs
if (!is.null(buildings)) {
  roof_vals <- mask(accum_poa, b_r)
  nonroof_vals <- mask(accum_poa, is.na(b_r))
  cat("Mean annual POA (kWh/m2) - roofs:", global(roof_vals, "mean", na.rm=TRUE)[1], "\n")
  cat("Mean annual POA (kWh/m2) - non-roofs:", global(nonroof_vals, "mean", na.rm=TRUE)[1], "\n")
}

# -------------------------
# End of script
# -------------------------
