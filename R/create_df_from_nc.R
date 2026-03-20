# #library(aermod)
# library(tidyverse)
# library(ecmwfr)
# library(lubridate)
# library(ncdf4)
# library(birk)


create_df_from_nc <- function(ncdf_dir,
                  variable,
                  sf_point,
                  start_date,
                  end_date){

  kr8_d8 <- format(seq(from = start_date, to = end_date, by = "month"))
  
  lat = sf::st_coordinates(sf::st_transform(sf::sf_point,4326))[2]
  lon = sf::st_coordinates(sf::st_transform(sf::sf_point,4326))[1]

  all_varz <- list()
    for (v in variable){

      all_d8s <- list()
      for (d in kr8_d8){

    m <- month(ymd(d))
    y <- year(ymd(d))

    nf <- ncdf4::nc_open(paste0(ncdf_dir, v, "_",sprintf("%02d",m),"_", y, ".nc"))

    thyme <- ncvar_get(nf, "valid_time")
    
    d8_time <- lubridate::ymd(str_sub(nf$dim$valid_time$units, 15, -1)) + lubridate::seconds(thyme)
    v_nam <- names(nf$var)[3]
    
    met_lon <- ncvar_get(nf, 'longitude')
    met_lat <- ncvar_get(nf, 'latitude')
    
    lon_ind <- which.closest(met_lon, lon)
    lat_ind <- which.closest(met_lat, lat)
    
    ##extract no2 variables for entire domain for 1 time step and altitude
    df <- data.frame(date = d8_time,
                        nam = ncvar_get(nf, v_nam, start = c(lon_ind,lat_ind,1), count = c(1,1,NROW(d8_time))))

    names(df) <- c("date", v_nam)

    ncdf4::nc_close(nf)

    all_d8s[[d]] <- df

    message("Processed ", v, " for ", d)

    }


  if (length(all_d8s) > 0) {
    all_yrz_df <- do.call(rbind, all_d8s)
    #write.csv(all_yrz_df, file = paste0(path_out, v_nam, "_", y, ".csv"), row.names = FALSE)
    all_varz[[v]] <- all_yrz_df
  }

}

all_variables <- Reduce(function(x, y) dplyr::left_join(x, y, by = "date"), all_varz)

return(all_variables)

}



