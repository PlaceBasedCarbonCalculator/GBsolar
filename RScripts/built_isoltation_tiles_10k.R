# Loop over every DSM tiles and produce hourly insulation maps. 
# DO this for a given month of the year
#source("R/insolation_calcs_daily.R")
source("R/insolation_calcs_intergrated.R")

grids = list.files("F:/DTM_DSM/GB_10k/DSM")
grids = gsub(".tiff","",grids)

#SJ19 has no ERA5 file
# Places with no land
grids = grids[!grids %in% c("SJ19","SS55","SX34","TF53","TF54","TM32") ]

ear5s = list.files("sampleData/ERA5/byGrid/")
ear5s = gsub(".Rds","",ear5s)
grids[!grids %in% ear5s]

# Check Already Done
#done = list.files("F:/DTM_DSM/GB_10k/solarDaily")
done = list.files("F:/DTM_DSM/GB_10k/solarAnnual")
done = unique(substr(done,1,4))
start_at = grep(done[length(done)],grids) + 1

day_of_month = 15

for(i in start_at:length(grids)){
  message(Sys.time()," ",grids[i])
  # insolation_daily(grid = grids[i],
  insolation_annual_strategy(grid = grids[i],
             day_of_month = day_of_month,
             year = 2020,
             dsm_dir = "F:/DTM_DSM/GB_10k/DSM",
             era5_dir = "sampleData/ERA5/byGrid/" ,
             #out_dir = "F:/DTM_DSM/GB_10k/solarDaily",
             out_dir = "F:/DTM_DSM/GB_10k/solarAnnual",
             nprocs = 35)
  
}
