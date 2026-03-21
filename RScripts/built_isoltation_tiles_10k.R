# Loop over every DSM tiles and produce hourly insulation maps. 
# DO this for a given month of the year
source("R/insolation_calcs.R")

grids = list.files("F:/DTM_DSM/GB_10k/DSM")
grids = gsub(".tiff","",grids)

day_of_month = 15

for(i in 1:length(grids)){
  message(Sys.time()," ",grids[i])
  insolation(grid = grids[i],
             day_of_month = day_of_month,
             year = 2020,
             dsm_dir = "F:/DTM_DSM/GB_10k/DSM", 
             ear5_dir = "sampleData/ERA5/byGrid/" ,
             out_dir = "F:/DTM_DSM/GB_10k/solar")
  
}
