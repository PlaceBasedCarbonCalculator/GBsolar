# Merge hourly into daily to save space
library(terra)

out_dir = "F:/DTM_DSM/GB_10k/solarDaily/"

grids = list.files("F:/DTM_DSM/GB_10k/DSM")
grids = gsub(".tiff","",grids)

raster_list = list.files("F:/DTM_DSM/GB_10k/solar/")
# remove non adjusted rasters
# raster_remove = raster_list[!grepl("_era5_adj.tif", raster_list)]
# unlink(file.path("F:/DTM_DSM/GB_10k/solar/",raster_remove))

day_of_month = 15

# Loop over Grids
for(i in 1:416){
  message(Sys.time()," ",grids[i])
  sub_list = raster_list[grepl(paste0(grids[i],"_ghi_2020[0-9]{2}",day_of_month),raster_list)]
  # Loop over days
  dates = unique(substr(sub_list, 10,17))
  for(j in seq_along(dates)){
    message(dates[j])
    day_list = sub_list[grepl(paste0("_",dates[j]),sub_list)]
    if(length(day_list) == 0){
      message("No rasters, skipping")
      next
    }
    # Loop over hours
    for(k in seq_along(day_list)){
      message(day_list[k])
      rast_hourly = terra::rast(file.path("F:/DTM_DSM/GB_10k/solar",day_list[k]))
      if (!exists("accum")) {
        accum <- rast_hourly         # initialise
      } else {
        accum <- accum + rast_hourly # running sum
      }
    }
    out_tif_corr = file.path("F:/DTM_DSM/GB_10k/solarDaily",paste0(grids[i],"_",dates[j],"_dailysum_era5_adj.tif"))
    terra::writeRaster(accum, 
                       filename=out_tif_corr, 
                       overwrite=TRUE, datatype="FLT4S", 
                       gdal="COMPRESS=LZW", NAflag=-9999)
    rm(accum, rast_hourly)
    # Remove hourly rasters
    unlink(file.path("F:/DTM_DSM/GB_10k/solar",day_list))
    # End of Day Loop
  }
  # End of Grid Loop
}
