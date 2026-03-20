
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