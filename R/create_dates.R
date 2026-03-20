
create_dates_df = function(start_date, end_date, time_zone, timestep){

# Define Areas ------------------------------------------------------------
start_d8 = as.POSIXct(start_date, tz="UTC")
end_d8 = as.POSIXct(end_date, tz="UTC")

d8_df = data.frame(date = seq(
  from=start_d8,
  to=end_d8,
  by=timestep
)  ) %>% 
  mutate(variable = paste0("X", sprintf("%02d", lubridate::hour(date)+1))) 


d8_df$dow = lubridate::wday(d8_df$date, label = TRUE)
d8_df$hr <- lubridate::hour(d8_df$date)+1
d8_df$doy <- lubridate::yday(d8_df$date)
d8_df$dow <- tolower(as.character(d8_df$dow))
d8_df$month <- tolower(as.character(lubridate::month(d8_df$date, label = TRUE)))

#d8_df = mutate(d8_df, date = format(date, "%d%m%Y %H%M"))

return(d8_df)

}
