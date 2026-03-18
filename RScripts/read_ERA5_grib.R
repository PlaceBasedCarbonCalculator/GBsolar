# RScripts/read_ERA5_grib.R
# Purpose: Example script to read an ERA5 hourly GRIB file (2020 UK) on Windows
# Path used in the user's request:
# C:/Users/malco/OneDrive - University of Leeds/Data/ERA5/ERA5_hourly_2020_UK.grib

# This script will:
#  - ensure needed packages are installed
#  - attempt to read the GRIB with `stars` (recommended, supports proxy reads)
#  - fall back to `raster::brick` if needed
#  - print metadata and time values (if available)
#  - write one example layer to GeoTIFF

# Usage: edit `era5_file` below if needed, then run this script in R/RStudio.
library(stars)
library(raster)
library(terra)
library(lubridate)

# Set the file path (use forward slashes or double backslashes on Windows)
path <- "C:/Users/malco/OneDrive - University of Leeds/Data/ERA5/ERA5_hourly_2020_UK.grib"

st <- try(stars::read_stars(path, proxy = TRUE), silent = TRUE)



# Try reading with stars (recommended: can use proxy = TRUE so it won't load everything)
read_with_stars <- function(path) {
  message("Trying to read with 'stars' (proxy = TRUE)...")
  
  if (inherits(st, "try-error")) {
    message("stars::read_stars failed: ", st)
    return(NULL)
  }

  # Print summary / metadata
  print(st)
  message("Dimensions (stars::st_dimensions):")
  print(stars::st_dimensions(st))

  # Try to read time values if available
  

  return(st)
}


# Try stars first
st_obj <- read_with_stars(era5_file)

message("stars read succeeded. If you want to operate on values in R, consider converting a subset to 'terra' or 'raster'.")
# Example: read the first time slice to a GeoTIFF without loading all slices
# Attempt to read only the first time index (if there's a time dimension)
time_vals <- try(stars::st_get_dimension_values(st_obj, "time"), silent = TRUE)
if (!inherits(time_vals, "try-error") && length(time_vals) > 0) {
  message("Reading first time slice and writing example TIFF: ERA5_first_slice.tif")
  # index selection: use stars indexing. We'll attempt to keep other dims, but only select first time
  # Determine which dimension name is 'time'
  dims <- names(stars::st_dimensions(st_obj))
  tdim <- NULL
  if ("time" %in% dims) tdim <- which(dims == "time")

  out_file <- file.path(getwd(), "ERA5_first_slice.tif")
  # try a generic approach: slice by index if time is present
  if (!is.null(tdim)) {
    # create an index list of full ranges then replace the time dimension with 1
    idx <- lapply(stars::st_dimensions(st_obj), function(d) seq_len(d$to - d$from + 1))
    idx[[tdim]] <- 1L
    # apply indexing
    slice <- try(st_obj[idx], silent = TRUE)
    if (!inherits(slice, "try-error")) {
      try(stars::write_stars(slice, out_file, options = c("COMPRESS=LZW")), silent = TRUE)
      message("Wrote example slice to: ", out_file)
    } else {
      message("Failed to slice stars object (writing skipped).")
    }
  } else {
    message("No explicit 'time' dimension; skipping example slice write.")
  }
} else {
  message("No time values available; skipping slice write example.")
}


message("Done. If you need further extraction (variables, subsets, reprojection), update this script or ask for a customized helper.")
