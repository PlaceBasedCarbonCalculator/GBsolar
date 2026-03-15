#read ERA5 NetCDF and align to DSM
library(terra)
library(data.table)
library(sf)

era5_nc  <- "era5_uk_2020_solar.nc"
dsm_tif  <- "tile_10km_dsm_27700.tif"

# Load DSM (target grid)
dsm <- rast(dsm_tif)              # EPSG:27700
crs(dsm) <- "EPSG:27700"

# Load ERA5 as SpatRaster
era5 <- rast(era5_nc)             # each variable × time
era5
# Check names – adapt if different
# e.g. "surface_solar_radiation_downwards", "total_sky_direct_solar_radiation_at_surface", ...

# Subset variables of interest
ghi_era5   <- era5$surface_solar_radiation_downwards
dni_era5   <- era5$total_sky_direct_solar_radiation_at_surface
t2m_era5   <- era5$`2m_temperature`
u10_era5   <- era5$`10m_u_component_of_wind`
v10_era5   <- era5$`10m_v_component_of_wind`

# ERA5 is lat/lon; reproject to BNG and resample to DSM grid
ghi_bng <- project(ghi_era5, dsm, method = "bilinear")
dni_bng <- project(dni_era5, dsm, method = "bilinear")

# ERA5 radiation is J/m² per hour → W/m²
ghi_bng <- ghi_bng / 3600
dni_bng <- dni_bng / 3600
names(ghi_bng) <- names(dni_bng) <- names(ghi_era5)
