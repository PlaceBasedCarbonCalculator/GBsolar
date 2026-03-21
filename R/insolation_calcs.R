# Function to produce annual isolation for a single tile.

insolation = function(grid = "SE23",
                      day_of_month = 15,
                      year = 2020,
                      dsm_dir = "", 
                      ear5_dir = "sampleData/ERA5/byGrid/" ,
                      out_dir = "C:/rastTemp/solar",
                      gisBase = "C:/Program Files/GRASS GIS 8.4"
                      ){
  
  # Check Paths
  stopifnot(dir.exists(out_dir))
  stopifnot(file.exists(file.path(era5_dir,paste0(grid,".Rds"))))
  stopifnot(file.exists(file.path(dsm_dir,paste0(grid,".tiff"))))
  
  day_of_month <- stringr::str_pad(day_of_month,2,pad = "0")
  
  # Load Weather Data
  era5 = readRDS(file.path(era5_dir,paste0(grid,".Rds")))
  
  # Process ERA5 Data
  era5 = era5[,c("timestamp","SSRD","TCC")]
 
  era5$era5_ghi_wm2 <- era5$SSRD / 3600
  era5$era5_tcc <- era5$TCC  # 0..1
  
  # Load DSM
  dsm <- terra::rast(file.path(dsm_dir,paste0(grid,".tiff")))
  
  # Process Times to run simulation
  # Adjust time for daylight hours
  # Simulation period (typical year)
  extent <- terra::ext(dsm)
  centre <- sf::st_point(c((extent[1] + extent[2]) / 2, (extent[3] + extent[4]) / 2))
  centre <- sf::st_as_sf(sf::st_sfc(list(centre), crs = 27700))
  centre <- sf::st_transform(centre, 4326)
  
  selected_days <- unlist(lapply(1:12, function(month) {
    month_days <- lubridate::ymd(paste0(year,"-",month,"-",day_of_month))
  }))
  selected_days <- as.Date(selected_days)
  selected_days <- as.POSIXct(selected_days)
  
  dawn <- crepuscule(centre, selected_days, solarDep = 6, direction = "dawn", POSIXct.out = TRUE)
  dusk <- crepuscule(centre, selected_days, solarDep = 6, direction = "dusk", POSIXct.out = TRUE)
  dawn$time <- lubridate::ceiling_date(dawn$time, "hour")
  dusk$time <- lubridate::floor_date(dusk$time, "hour")
  
  times <- list()
  for(i in 1:nrow(dawn)){
    times[[i]] <- seq(dawn$time[i], dusk$time[i], by = "1 hour")
  }
  times <- unlist(times)
  times <- as.POSIXct(times, tz = "UTC")
  
  # -------------------------
  # 2. Derive slope (tilt) and aspect (azimuth) rasters
  # -------------------------
  # Use terra::terrain to compute slope/aspect in degrees
  slope_r <- terra::terrain(dsm, v = "slope", unit = "degrees", neighbors = 8)
  aspect_r <- terra::terrain(dsm, v = "aspect", unit = "degrees", neighbors = 8)

  # -------------------------
  # 4. Initialize GRASS and import DSM
  # -------------------------
  loc = rgrass::initGRASS(gisBase = gisBase, 
                  home = tempdir(), 
                  gisDbase = file.path(tempdir(),"grassdb"),
                  mapset = "PERMANENT",
                  override = TRUE)
  
  # Create location with EPSG:27700
  rgrass::execGRASS("g.proj", flags = "c", epsg = 27700)
  
  # Export DSM to a temporary GeoTIFF and import into GRASS
  tmp_dsm <- file.path(tempdir(), "dsm_for_grass.tif")
  terra::writeRaster(dsm, tmp_dsm, overwrite = TRUE)
  rgrass::execGRASS("r.in.gdal", flags = c("o","overwrite"), input = tmp_dsm, output = "dsm")
  
  # Import slope/aspect if desired
  tmp_slope <- file.path(tempdir(), "slope_for_grass.tif")
  terra::writeRaster(slope_r, tmp_slope, overwrite = TRUE)
  rgrass::execGRASS("r.in.gdal", flags = c("o","overwrite"), input = tmp_slope, output = "slope")
  
  tmp_aspect <- file.path(tempdir(), "aspect_for_grass.tif")
  terra::writeRaster(aspect_r, tmp_aspect, overwrite = TRUE)
  rgrass::execGRASS("r.in.gdal", flags = c("o","overwrite"), input = tmp_aspect, output = "aspect")
  
  # -------------------------
  # 7. Run r.sun hourly to compute beam/diffuse/global on horizontal surface
  # -------------------------
  # r.sun expects day-of-year and time (decimal hours). We'll loop hourly.
  # Note: r.sun can compute beam/diffuse/reflected using horizon raster.
  # For speed, run only daylight hours per day; here we run full 0-23 for simplicity.
  
  # Return to 2m resolution
  rgrass::execGRASS("g.region", raster = "dsm", flags = "p") # back to fine region
  
  # Create a vector of day-of-year and hour
  doy <- as.integer(format(times, "%j"))
  hour_decimal <- as.numeric(format(times, "%H")) + as.numeric(format(times, "%M"))/60
  
  # We'll store hourly global horizontal irradiance (GHI) raster outputs in a temp mapset
  # For large runs, parallelize by day or tile and write to disk incrementally.
  for (i in seq_along(times)) {
    t <- times[i]
    day_i <- doy[i]
    hour_i <- hour_decimal[i]
    #message(Sys.time()," ",t)
    out_prefix <- paste0("sun_", format(t, "%Y%m%d%H"))
    # r.sun parameters: day, time, horizon, beam_rad, diff_rad, glob_rad
    rgrass::execGRASS("r.sun",
              flags = c("overwrite"),
              parameters = list(elevation = "dsm",
                                aspect = "aspect",
                                slope = "slope",
                                #horizon_basename = "horizon",
                                #horizon_step = 1,
                                day = day_i,
                                time = hour_i,
                                nprocs = 30, # Use multiple cores
                                beam_rad = paste0(out_prefix, "_beam"),
                                diff_rad = paste0(out_prefix, "_diff"),
                                glob_rad = paste0(out_prefix, "_glob")
              ))
    
    # Export glob to GeoTIFF for later transposition
    out_tif <- file.path(out_dir, paste0(grid,"_ghi_", format(t, "%Y%m%d%H"), ".tif"))
    rgrass::execGRASS("r.out.gdal", flags = c("overwrite"), input = paste0(out_prefix, "_glob"),
              output = out_tif, format = "GTiff", createopt = "COMPRESS=LZW",
              type = "Float32",  nodata = -9999)
    
    
    era5_row <- era5[era5$timestamp == t, ]
    if (nrow(era5_row) == 1) {
      era5_ghi <- era5_row$era5_ghi_wm2
      # optional: include cloud effect multiplier
      cloud_factor <- pmax(0.1, 1 - era5_row$TCC) # fallback
      ghi_r <- terra::rast(out_tif)
      mean_ghi <- terra::global(ghi_r, fun='mean', na.rm=TRUE)[1,1]
      if (!is.na(mean_ghi) && mean_ghi > 0) {
        # bias-corrected raster
        scale <- era5_ghi / mean_ghi
        ghi_corr <- ghi_r * scale * cloud_factor
        out_tif_corr <- file.path(out_dir, paste0(grid,"_ghi_", format(t, "%Y%m%d%H"), "_era5_adj.tif"))
        terra::writeRaster(ghi_corr, filename=out_tif_corr, overwrite=TRUE, datatype="FLT4S", 
                    gdal="COMPRESS=LZW", NAflag=-9999)
      }
    } else {
      stop("Muliple or missing rows in ERA5")
    }
  }
  
  
  # Clean Up
  unlink(tmp_dsm, tmp_slope, tmp_aspect)
  unlink(file.path(tempdir(), "grassdb"), recursive = TRUE)

  return(invisible(NULL))
  
}