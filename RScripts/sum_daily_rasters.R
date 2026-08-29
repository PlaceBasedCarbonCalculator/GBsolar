library(terra)

path_in = "F:/DTM_DSM/GB_10k/solarDaily"
path_out = "F:/DTM_DSM/GB_10k/solarTotal"

fls = list.files(file.path(path_in), full.names = FALSE, pattern = ".tif")

grids = as.data.frame(table(substr(fls,1,4)))

if(length(unique(grids$Freq)) != 1){
  stop("Differnt number of day for each grid")
}

if(!dir.exists(path_out)){dir.create(path_out)}

# Days per month (non-leap year)
month_days = c(31,28,31,30,31,30,31,31,30,31,30,31)

overwrite = FALSE
performance_ratio = 0.85 # Fraction of Solar that becomes Electricity

for(i in 1:nrow(grids)){
  
  if(file.exists(file.path(path_out,paste0(grids$Var1[i],".tif"))) & !overwrite){
    message("Skipping ",grids$Var1[i])
    next
  }
  message(Sys.time()," Processing ",grids$Var1[i])
  fls_sub = fls[grepl(grids$Var1[i],fls)]
  
  # Ensure correct ordering (important!)
  fls_sub = sort(fls_sub)
  
  # Load rasters
  r_stack = rast(file.path(path_in,fls_sub))
  
  # Check we have 12 months
  if(nlyr(r_stack) != 12){
    stop(paste("Expected 12 rasters, got", nlyr(r_stack), "for", grid_id))
  }
  
  # Apply monthly weighting
  # Sum to annual total (kWh/m²/year)
  # Annual average per day (optional)
  # 0.85 Perfo ratio
  annual_total = sum(r_stack * month_days) * performance_ratio / 1000
  #annual_avg =  annual_total / sum(month_days) 
  
  writeRaster(annual_total, file.path(path_out,paste0(grids$Var1[i],".tif")), overwrite= overwrite)
  
}

