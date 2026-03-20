library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(data.table)

# Set the file path (use forward slashes or double backslashes on Windows)
path <- "D:/OneDrive - University of Leeds/Data/ERA5/ERA5_hourly_2020_UK.grib"
bounds <- readRDS("../build/_targets/objects/bounds_la")
bounds <- bounds[substr(bounds$LAD25CD,1,1) != "N",]
bounds <- st_simplify(bounds, TRUE, 1000)
bounds <- st_buffer(bounds, 1000)
bounds <- st_union(bounds)

# Make Grid
os_grid_10 = read_sf("../../OrdnanceSurvey/OS-British-National-Grids/os_bng_grids/os_bng_grids.gpkg", layer = "10km_grid")
os_grid_10 = os_grid_10[bounds,] #9100 squares -> 2903
os_cents = st_centroid(os_grid_10)


r <- rast(path)  # GRIB; lazy connection
os_cents_rascrs <- project(vect(os_cents), crs(r))  # project points to raster CRS


bands <- sf::gdal_utils(source = path,
                        options = c("-json"), 
                        quiet = T) %>%
  jsonlite::fromJSON()

bands_df = bands$bands
#bands_df$timestamp = unlist((bands_df$metadata[[1]]$GRIB_REF_TIME))
#bands_df$timestamp = as.POSIXct(as.integer(bands_df$timestamp), origin="1970-01-01")

bands_df$timestamp = unlist((bands_df$metadata[[1]]$GRIB_VALID_TIME))
bands_df$timestamp = as.POSIXct(as.integer(bands_df$timestamp), origin="1970-01-01")

bands_df$description = bands$bands$description

bands_df$comment = unlist((bands_df$metadata[[1]]$GRIB_COMMENT))
bands_df$element = unlist((bands_df$metadata[[1]]$GRIB_ELEMENT))
bands_df$unit = unlist((bands_df$metadata[[1]]$GRIB_UNIT))

bands_df = bands_df[,c("band","timestamp","comment","element","unit")]

bands_df$element = gsub(" of table 228 of center ECMWF","",bands_df$element)

bands_df$dups = (duplicated(bands_df[,c("timestamp2","element")]))

#bands_df = bands_df[bands_df$element %in% c("SSRC","SSRD","2T","SP","TCC","10V","10U")]

# Using the same 'bands_df' from your stars metadata (band -> timestamp/element)
types <- unique(bands_df$element)
bands_df <- bands_df |> mutate(month = format(timestamp, "%Y-%m"))

extract_batch_terra <- function(idx) {
  r_sub <- r[[idx]]                 # subset the layers
  m <- terra::extract(r_sub, os_cents_rascrs)  # data.frame: ID + one column per layer
  m$ID <- NULL
  dt <- as.data.table(m)
  setnames(dt, paste0("b", idx))
  dt[, grid := os_cents$tile_name]
  long <- melt(dt, id.vars = "grid", variable.name = "b", value.name = "value")
  long[, band := as.integer(sub("^b", "", b))]
  long[, b := NULL]
  long <- merge(long,
                bands_df[, c("band","timestamp","element")],
                by = "band",
                all.x = TRUE,
                sort = FALSE)
  long[, band := NULL]
  tibble::as_tibble(long)
}

data_all <- purrr::map_dfr(types, function(tp) {
  idx_tp <- which(bands_df$element == tp)
  idx_by_month <- split(idx_tp, bands_df$month[idx_tp])
  purrr::map_dfr(idx_by_month, extract_batch_terra)
}, .progress = TRUE)

# data_all$grid = as.factor(data_all$grid)
# data_all$element = as.factor(data_all$element)


data_wide <- data_all |>
  pivot_wider(id_cols = c("grid","timestamp"),
              values_from = "value",
              names_from = "element")


saveRDS(data_wide,"sampleData/ERA5/summary_data_all.Rds")

data_wide_lst <- data_wide |>
  group_by(grid) |>
  group_split()

for(i in 1:length(data_wide_lst)){
  sub <- data_wide_lst[[i]]
  saveRDS(sub, paste0("sampleData/ERA5/byGrid/",sub$grid[1],".Rds"))
}
