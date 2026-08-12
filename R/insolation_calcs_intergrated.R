insolation_annual_strategy <- function(
    grid = "SE23",
    day_of_month = 15,
    year = 2020,
    dsm_dir = "F:/DTM_DSM/GB_10k/DSM",
    era5_dir = "sampleData/ERA5/byGrid/",
    out_dir = "F:/DTM_DSM/GB_10k/solarAnnual",
    gisBase = "C:/Program Files/GRASS GIS 8.4",
    nprocs = 35
) {
  
  # -------------------------
  # Checks
  # -------------------------
  stopifnot(dir.exists(out_dir))
  stopifnot(file.exists(file.path(era5_dir, paste0(grid, ".Rds"))))
  stopifnot(file.exists(file.path(dsm_dir, paste0(grid, ".tiff"))))
  
  day_of_month <- stringr::str_pad(day_of_month, 2, pad = "0")
  
  # -------------------------
  # Load ERA5 and compute annual energy
  # -------------------------
  era5 <- readRDS(file.path(era5_dir, paste0(grid, ".Rds")))
  era5 <- era5[, c("timestamp", "SSRD")]
  
  # ERA5 annual energy (Wh/m²/year)
  E_ERA5_annual <- sum(era5$SSRD, na.rm = TRUE) / 3600
  
  # -------------------------
  # Load DSM
  # -------------------------
  dsm <- terra::rast(file.path(dsm_dir, paste0(grid, ".tiff")))
  
  # -------------------------
  # Compute slope and aspect
  # -------------------------
  slope_r  <- terra::terrain(dsm, v = "slope",  unit = "degrees")
  aspect_r <- terra::terrain(dsm, v = "aspect", unit = "degrees")
  
  # -------------------------
  # Init GRASS
  # -------------------------
  rgrass::initGRASS(
    gisBase = gisBase,
    home = tempdir(),
    gisDbase = file.path(tempdir(), "grassdb"),
    mapset = "PERMANENT",
    override = TRUE
  )
  
  rgrass::execGRASS("g.proj", flags = "c", epsg = 27700)
  
  # -------------------------
  # Import rasters into GRASS
  # -------------------------
  tmp_dsm    <- file.path(tempdir(), "dsm.tif")
  tmp_slope  <- file.path(tempdir(), "slope.tif")
  tmp_aspect <- file.path(tempdir(), "aspect.tif")
  
  terra::writeRaster(dsm,    tmp_dsm,    overwrite = TRUE)
  terra::writeRaster(slope_r,  tmp_slope,  overwrite = TRUE)
  terra::writeRaster(aspect_r, tmp_aspect, overwrite = TRUE)
  
  rgrass::execGRASS("r.in.gdal", flags = c("o","overwrite"),
                    input = tmp_dsm, output = "dsm")
  rgrass::execGRASS("r.in.gdal", flags = c("o","overwrite"),
                    input = tmp_slope, output = "slope")
  rgrass::execGRASS("r.in.gdal", flags = c("o","overwrite"),
                    input = tmp_aspect, output = "aspect")
  
  rgrass::execGRASS("g.region", raster = "dsm")
  
  # -------------------------
  # Representative days (mid‑month)
  # -------------------------
  months <- 1:12
  rep_days <- as.Date(
    lubridate::ymd(paste0(year, "-", months, "-", day_of_month))
  )
  doy <- as.integer(format(rep_days, "%j"))
  days_in_month <- lubridate::days_in_month(rep_days)
  
  # -------------------------
  # Run r.sun (daily integration, clear sky)
  # -------------------------
  clear_sky_monthly <- vector("list", 12)
  
  for (i in seq_along(doy)) {
    
    out_map <- paste0("sun_clear_", sprintf("%02d", i))
    
    rgrass::execGRASS(
      "r.sun",
      flags = "overwrite",
      parameters = list(
        elevation = "dsm",
        slope = "slope",
        aspect = "aspect",
        day = doy[i],
        step = 1,                 # hourly integration
        nprocs = nprocs,
        glob_rad = out_map
      )
    )
    
    # Export to terra
    out_tif <- file.path(tempdir(), paste0(out_map, ".tif"))
    rgrass::execGRASS(
      "r.out.gdal",
      flags = "overwrite",
      input = out_map,
      output = out_tif,
      format = "GTiff",
      type = "Float32",
      nodata = -9999
    )
    
    clear_sky_monthly[[i]] <- terra::rast(out_tif)
  }
  
  # -------------------------
  # Build clear‑sky annual energy (Wh/m²/year)
  # -------------------------
  E_clear_annual <- clear_sky_monthly[[1]] * days_in_month[1]
  
  for (i in 2:12) {
    E_clear_annual <- E_clear_annual +
      clear_sky_monthly[[i]] * days_in_month[i]
  }
  
  # -------------------------
  # Compute scaling factor from ERA5
  # -------------------------
  mean_clear <- terra::global(E_clear_annual, "mean", na.rm = TRUE)[1,1]
  scale_factor <- E_ERA5_annual / mean_clear
  
  # -------------------------
  # Final annual solar energy
  # -------------------------
  # kWh/m2/year - I think (Leeds is ~100 to 1700 kWh, mean is 1000 kWh)
  # Solar Atlas Gives 983 kWh/m2 for Leeds
  E_final_annual <- E_clear_annual * scale_factor / 1000
  
  # -------------------------
  # Write output (Wh/m²/year)
  # -------------------------
  out_file <- file.path(
    out_dir,
    paste0(grid, "_annual_insolation_Whm2.tif")
  )
  
  terra::writeRaster(
    E_final_annual,
    filename = out_file,
    overwrite = TRUE,
    datatype = "FLT4S",
    gdal = "COMPRESS=LZW",
    NAflag = -9999
  )
  
  # -------------------------
  # Clean up
  # -------------------------
  unlink(c(tmp_dsm, tmp_slope, tmp_aspect))
  unlink(file.path(tempdir(), "grassdb"), recursive = TRUE)
  
  invisible(out_file)
}