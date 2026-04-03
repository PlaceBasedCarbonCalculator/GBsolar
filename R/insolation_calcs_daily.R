insolation_daily <- function(grid = "SE23",
                             day_of_month = 15,
                             year = 2020,
                             dsm_dir = "", 
                             era5_dir = "sampleData/ERA5/byGrid/" ,
                             out_dir = "C:/rastTemp/solar",
                             gisBase = "C:/Program Files/GRASS GIS 8.4",
                             nprocs = 35,
                             skip = TRUE
){
  
  # Check Paths
  stopifnot(dir.exists(out_dir))
  stopifnot(file.exists(file.path(era5_dir,paste0(grid,".Rds"))))
  stopifnot(file.exists(file.path(dsm_dir,paste0(grid,".tiff"))))
  
  day_of_month <- stringr::str_pad(day_of_month,2,pad="0")
  
  # Load Weather Data
  era5 <- readRDS(file.path(era5_dir,paste0(grid,".Rds")))
  era5 <- era5[,c("timestamp","SSRD","TCC")]
  era5$era5_ghi_wm2 <- era5$SSRD / 3600
  era5$era5_tcc <- era5$TCC
  
  # Load DSM
  dsm <- terra::rast(file.path(dsm_dir,paste0(grid,".tiff")))
  
  # Location centre for daylight filter
  extent <- terra::ext(dsm)
  centre <- sf::st_point(c((extent[1] + extent[2])/2,
                           (extent[3] + extent[4])/2))
  centre <- sf::st_as_sf(sf::st_sfc(list(centre), crs = 27700))
  centre <- sf::st_transform(centre, 4326)
  
  # Selected days = 15th of each month
  selected_days <- lubridate::ymd(
    paste0(year,"-",1:12,"-",day_of_month)
  ) |> as.POSIXct()
  
  dawn <- suntools::crepuscule(centre, selected_days,
                               solarDep = 6, direction = "dawn", POSIXct.out = TRUE)
  dusk <- suntools::crepuscule(centre, selected_days,
                               solarDep = 6, direction = "dusk", POSIXct.out = TRUE)
  dawn$time <- lubridate::ceiling_date(dawn$time, "hour")
  dusk$time <- lubridate::floor_date(dusk$time, "hour")
  
  # Create list of times per day
  all_times <- unlist(lapply(1:nrow(dawn), function(i){
    seq(dawn$time[i], dusk$time[i], by="hour")
  }))
  all_times <- as.POSIXct(all_times, tz = "UTC")
  
  # Compute slope/aspect
  slope_r <- terra::terrain(dsm, v = "slope", unit = "degrees", neighbors = 8)
  aspect_r <- terra::terrain(dsm, v = "aspect", unit = "degrees", neighbors = 8)
  
  # Init GRASS
  rgrass::initGRASS(gisBase = gisBase,
                    home = tempdir(),
                    gisDbase = file.path(tempdir(),"grassdb"),
                    mapset = "PERMANENT",
                    override = TRUE)
  rgrass::execGRASS("g.proj", flags = "c", epsg = 27700)
  
  # Import DSM, slope, aspect
  tmp_dsm <- file.path(tempdir(),"dsm_for_grass.tif")
  terra::writeRaster(dsm, tmp_dsm, overwrite = TRUE)
  rgrass::execGRASS("r.in.gdal", flags=c("o","overwrite"), input=tmp_dsm, output="dsm")
  
  tmp_slope <- file.path(tempdir(),"slope_for_grass.tif")
  terra::writeRaster(slope_r, tmp_slope, overwrite=TRUE)
  rgrass::execGRASS("r.in.gdal", flags=c("o","overwrite"), input=tmp_slope, output="slope")
  
  tmp_aspect <- file.path(tempdir(),"aspect_for_grass.tif")
  terra::writeRaster(aspect_r, tmp_aspect, overwrite=TRUE)
  rgrass::execGRASS("r.in.gdal", flags=c("o","overwrite"), input=tmp_aspect, output="aspect")
  
  rgrass::execGRASS("g.region", raster="dsm")
  
  # Prepare day/hour vectors
  doy  <- as.integer(format(all_times, "%j"))
  hour <- as.numeric(format(all_times,"%H")) + 
    as.numeric(format(all_times,"%M"))/60
  
  # ---- DAILY ACCUMULATION ----
  current_day <- NULL
  daily_accum <- NULL
  
  for(i in seq_along(all_times)){
    
    t <- all_times[i]
    day_i  <- doy[i]
    hour_i <- hour[i]
    
    day_date <- as.Date(t)
    
    # If new day starts → write previous day summary
    if (!is.null(current_day) && day_date != current_day) {
      # write daily sum raster
      daily_out <- file.path(out_dir, paste0(grid,"_",current_day,"_daily_total.tif"))
      terra::writeRaster(daily_accum, daily_out,
                         overwrite=TRUE, datatype="FLT4S",
                         gdal="COMPRESS=LZW", NAflag=-9999)
      daily_accum <- NULL
    }
    
    if (is.null(current_day)) current_day <- day_date
    current_day <- day_date
    
    # run r.sun
    out_prefix <- paste0("sun_", format(t, "%Y%m%d%H"))
    rgrass::execGRASS("r.sun",
                      flags="overwrite",
                      parameters=list(
                        elevation="dsm",
                        aspect="aspect",
                        slope="slope",
                        day=day_i,
                        time=hour_i,
                        nprocs=nprocs,
                        glob_rad=paste0(out_prefix,"_glob")
                      ))
    
    # export hourly raster
    out_tif <- file.path(tempdir(), paste0(grid,"_ghi_",format(t,"%Y%m%d%H"),".tif"))
    rgrass::execGRASS("r.out.gdal", flags="overwrite",
                      input=paste0(out_prefix,"_glob"),
                      output=out_tif, format="GTiff",
                      type="Float32", nodata=-9999)
    
    ghi_r <- terra::rast(out_tif)
    unlink(out_tif)
    
    # scale with ERA5 bias correction
    era5_row <- era5[era5$timestamp == t,]
    scale <- era5_row$era5_ghi_wm2 / mean(values(ghi_r), na.rm=TRUE)
    cloud_factor <- pmax(0.1, 1 - era5_row$TCC)
    ghi_corr <- ghi_r * scale * cloud_factor
    
    # ---- accumulate into daily raster ----
    if (is.null(daily_accum)) {
      daily_accum <- ghi_corr
    } else {
      daily_accum <- daily_accum + ghi_corr
    }
  }
  
  # Write final day
  if (!is.null(daily_accum)) {
    last_out <- file.path(out_dir, paste0(grid,"_",current_day,"_daily_total.tif"))
    terra::writeRaster(daily_accum, last_out,
                       overwrite=TRUE, datatype="FLT4S",
                       gdal="COMPRESS=LZW", NAflag=-9999)
  }
  
  # Clean Up
  unlink(c(tmp_dsm, tmp_slope, tmp_aspect))
  unlink(file.path(tempdir(), "grassdb"), recursive = TRUE)
  
  return(invisible(NULL))
}