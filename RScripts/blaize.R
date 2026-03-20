#remotes::install_github("Pakillo/CityShadeMapper")
library(CityShadeMapper)
library(terra)
library(sf)
library(dplyr)
library(suncalc)
library(ecmwfr)
library(ncdf4)
library(birk)

#source("R/calc.R")
source("R/create_dates.R")
source("R/create_df_from_nc.R")
#source("R/get.R")

# # define site
# site = "Manchester"
# # site directory
# sdir = paste0("sites/",site,"/")
# #read in point for site
# pt = readRDS(paste0("sites/",site,"/DAT/site_location.RDS"))
# # read in date range
#date_df = readRDS(paste0("sites/",site,"/DAT/date_range.RDS"))[1:2,]

date_df = create_dates_df("2020-01-01","2020-12-31","Europe/London","1 hour")
date_df$date <- lubridate::with_tz(date_df$date, tzone = "Europe/London")

# read in DSM data
dsm = rast("sampleData/DSM/SE23.tiff")

# calculate the slope and aspect
slope_r = terrain(dsm, v = "slope", unit = "degrees", neighbors = 8)
aspect_r = terrain(dsm, v = "aspect", unit = "degrees", neighbors = 8)

# shade function requires name to be height_canopy
names(dsm) = "height_canopy"

# returns percentage of irradiation in local time
# started 13:57
shade = make_shademap(height.ras = dsm, date = date_df$date, hour = 0:23,omit.nights = FALSE, type = "canopy", multicore = TRUE)

# get lat lon of centroid of domain
dsm_domain = st_as_sf(as.polygons(ext(dsm_ll), crs = "EPSG:4326"))
domain_cent = st_centroid(dsm_domain)
lat = st_coordinates(domain_cent)[2]
lon = st_coordinates(domain_cent)[1]

# for periods that straddle clock shifts or other countries the time offset will be needed, not implemented yet as is GMT/UTC
# Find time zone from lon-lat
tz = lutz::tz_lookup_coords(lat = lat, lon = lon, method = "fast", warn = FALSE)

# Get time difference with UTC
time.offset = lutz::tz_offset(date_df$date, tz = tz)$utc_offset_h

# variables to import
variablez = c("surface_solar_radiation_downwards","2m_temperature","total_sky_direct_solar_radiation_at_surface")

# extract ERA data for site location
all_variables = create_df_from_nc(ncdf_dir = paste0(sdir,"MET/"), variable = variablez, sf_point = pt, start_date = min(date_df$date), end_date = max(date_df$date))

# remove cumulative for ssrd
all_variables$ssrd = c(all_variables$ssrd[1], diff(all_variables$ssrd))
all_variables$ssrd[all_variables$ssrd < 0] = 0
all_variables$ssrd[1] = 0

# create data frame of meteo values
meteo_dat = all_variables |>
  transmute(date,
            ssrd_wm2 = ssrd/3600,
            tdir_wm2 = fdir/3600,
            air_temp = t2m-273.13) |> 
  mutate(tdiff_wm2 = ssrd_wm2-tdir_wm2)

# calculate tilted wm2 for each grid cell in domain
tilted_wm2 = calc_solar_potential(meteo_dat = meteo_dat,lat = lat,lon = lon,aspect_rast = aspect_r,slope_rast = slope_r,shade_rast = shade)

# plot first day
plot(tilted_wm2[[9:19]])

# plot total
total_wm2 = app(tilted_wm2,sum)
plot(total_wm2)

# calculate estimated yield for some estimated PV parameters and using air temp top calculate cell efficiency
# output is a multi layered raster for each hour of the period
yield_kWh = calc_pv_potential(tilted_wm2 = tilted_wm2,meteo_dat = meteo_dat)

