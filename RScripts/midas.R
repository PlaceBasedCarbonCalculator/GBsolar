midas_csv <- "midas_hourly_solar_2020.csv"
midas <- fread(midas_csv)

# Example column names – adjust to your file
# station_id, latitude, longitude, obs_time, global_irradiance
midas[, obs_time := as.POSIXct(obs_time, tz = "UTC")]

# Convert to W/m² if needed (if J/m² per hour)
if (attr(midas$global_irradiance, "units") == "J/m^2") {
  midas[, ghi_wm2 := global_irradiance / 3600]
} else {
  midas[, ghi_wm2 := global_irradiance]
}

#3. Match MIDAS stations to ERA5 grid and compute monthly bias
# Pick stations near your tile; simplest is nearest ERA5 grid cell.


# Build ERA5 grid cell centroids in lat/lon
era5_ll <- rast(era5_nc)          # original lat/lon grid
era5_pts <- as.points(era5_ll[[1]])  # one layer is enough
era5_sf  <- st_as_sf(era5_pts)

# Convert MIDAS station coords to sf
midas_sf <- st_as_sf(midas[, .(station_id, latitude, longitude)],
                     coords = c("longitude", "latitude"), crs = 4326)
midas_sf <- st_unique(midas_sf)

# Nearest ERA5 cell for each station
idx <- st_nearest_feature(midas_sf, era5_sf)
midas_sf$era5_cell <- idx

# Join back to MIDAS table
midas <- merge(midas, st_drop_geometry(midas_sf), by = "station_id")

# Extract ERA5 GHI time series at those ERA5 cells:

# Extract ERA5 GHI at those cells (lat/lon grid)
cell_ids <- unique(midas$era5_cell)
era5_ghi_ts <- terra::extract(ghi_era5, era5_pts[cell_ids, ], cells = TRUE)

# era5_ghi_ts: data.frame(cell, layer1, layer2, ...)
# Build long table: time × cell
era5_time <- as.POSIXct(terra::time(ghi_era5), tz = "UTC")
era5_long <- data.table::as.data.table(era5_ghi_ts)
era5_long <- melt(era5_long, id.vars = "cell",
                  variable.name = "layer", value.name = "ghi_wm2_era5")
era5_long[, time := era5_time[as.integer(gsub("lyr", "", layer))]]
era5_long[, layer := NULL]

# Join MIDAS and ERA5 by station/time:

midas_era5 <- merge(
  midas[, .(station_id, era5_cell, obs_time, ghi_wm2)],
  era5_long,
  by.x = c("era5_cell", "obs_time"),
  by.y = c("cell", "time"),
  all = FALSE
)

# Compute monthly bias ratio: MIDAS / ERA5
midas_era5[, month := format(obs_time, "%m")]
bias_month <- midas_era5[
  !is.na(ghi_wm2) & !is.na(ghi_wm2_era5) & ghi_wm2_era5 > 0,
  .(bias = median(ghi_wm2 / ghi_wm2_era5)),
  by = month
]
bias_month

#4. Apply bias correction to ERA5 rasters on the DSM grid
#We now have a monthly multiplicative bias factor for GHI. Apply it to the reprojected ERA5 GHI/DNI rasters:

# ghi_bng has one layer per hour; get its time vector
tvec <- as.POSIXct(time(ghi_bng), tz = "UTC")
months <- format(tvec, "%m")

# Build a vector of bias factors per layer
bias_vec <- bias_month$bias[match(months, bias_month$month)]

# Apply bias to each layer
for (i in seq_along(bias_vec)) {
  ghi_bng[[i]] <- ghi_bng[[i]] * bias_vec[i]
  dni_bng[[i]] <- dni_bng[[i]] * bias_vec[i]   # same factor applied to DNI
}
