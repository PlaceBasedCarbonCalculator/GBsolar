# Annual insolation, normalised by a clear-sky index rather than by each
# tile's own mean.
#
# WHAT WAS WRONG
#
# insolation_annual_strategy() in R/insolation_calcs_intergrated.R ended with
#
#     mean_clear   <- terra::global(E_clear_annual, "mean", na.rm = TRUE)[1,1]
#     scale_factor <- E_ERA5_annual / mean_clear
#     E_final      <- E_clear_annual * scale_factor / 1000
#
# `mean_clear` is the spatial mean of modelled clear-sky energy over that one
# 10 km tile, so every tile was rescaled until its own mean equalled its own
# ERA5 total. Measured on the outputs, that holds exactly: for fourteen tiles
# checked, mean(E_final) == E_ERA5_annual / 1000 to within 0.000%.
#
# Two consequences.
#
# First, terrain cancels out of every tile mean. A tile of steep north-facing
# ground has a low mean_clear, so it gets a large multiplier and is inflated
# back up to the ERA5 mean; a flat tile is pushed down. Real differences in
# solar resource between one 10 km square and the next - which is much of the
# point of a terrain-based solar map - are erased, leaving only ERA5's own
# ~1% variation. Relative pattern *within* a tile survives; between tiles it
# does not.
#
# Second, because the multiplier is constant within a tile and jumps at the
# edge, the map shows the 10 km grid. Measured across ten boundaries: mean
# step 5.5%, maximum 14.2%, against 0.04% across an arbitrary line inside a
# tile.
#
# WHAT THIS DOES INSTEAD
#
# The ERA5 ratio is meant to convert modelled clear-sky radiation into
# realistic all-sky radiation - it is a cloud correction. That is the
# clear-sky index, and its denominator has to be clear-sky radiation on a flat
# horizontal surface, which is what a reanalysis grid cell reports:
#
#     CSI       = E_ERA5_allsky / E_clear_FLAT
#     E_final   = E_clear_terrain * CSI
#
# Dividing by the terrain-averaged mean instead conflates cloud cover with
# terrain shading, which is exactly how terrain came to cancel. With a flat
# reference the denominator is terrain-independent, so terrain survives into
# the result, and because both CSI terms vary smoothly with latitude the
# multiplier no longer jumps at tile edges.
#
# E_clear_FLAT is cheap: it needs no shadows, so it is computed on a coarse
# region (500 m cells, ~20x20 for a 10 km tile) with constant elevation and
# zero slope. That is a few seconds against the hours the terrain run takes.
#
# E_ERA5_allsky comes from R/era5_field.R, which smooths the reanalysis over
# GB first. The per-square extracts are nearest-neighbour samples of a ~31 km
# grid - 38% of east-west neighbours hold an identical value and 32% step by
# more than 1% - so used raw they would print the ERA5 cell edges onto the map.
#
# THE CLEAR-SKY RASTER IS NOW SAVED
#
# E_clear_terrain is the expensive product and it was previously discarded,
# keeping only the normalised result. That is why this correction needs a full
# re-run: the normalisation cannot be undone after the fact, because
# mean_clear was never recorded. Saving it means any future change to the
# normalisation is a few minutes of arithmetic rather than another pass over
# 2,253 tiles.

#' Annual insolation for one OS 10 km square, clear-sky-index normalised
#'
#' @param grid OS 10 km square, e.g. "SE33"
#' @param era5_field data.frame from R/era5_field.R with `grid` and
#'   `era5_smooth` (Wh/m2/year). Read once and passed in, not re-read per tile.
#' @param clear_dir where to save the terrain clear-sky raster. Set to NULL to
#'   skip saving it, which is not recommended - see above.
#' @param flat_res cell size in metres for the flat reference run.
insolation_annual_csi <- function(
    grid = "SE23",
    era5_field,
    day_of_month = 15,
    year = 2020,
    dsm_dir = "F:/DTM_DSM/GB_10k/DSM",
    out_dir = "F:/DTM_DSM/GB_10k/solarAnnualCSI",
    clear_dir = "F:/DTM_DSM/GB_10k/solarClearSky",
    gisBase = "C:/Program Files/GRASS GIS 8.4",
    nprocs = 35,
    flat_res = 500
) {

  stopifnot(dir.exists(out_dir))
  if (!is.null(clear_dir)) dir.create(clear_dir, showWarnings = FALSE, recursive = TRUE)
  stopifnot(file.exists(file.path(dsm_dir, paste0(grid, ".tiff"))))

  i_era5 <- match(grid, era5_field$grid)
  if (is.na(i_era5)) stop("No ERA5 value for grid ", grid)
  E_ERA5_annual <- era5_field$era5_smooth[i_era5]
  if (!is.finite(E_ERA5_annual) || E_ERA5_annual <= 0)
    stop("Bad ERA5 annual total for ", grid, ": ", E_ERA5_annual)

  day_of_month <- stringr::str_pad(day_of_month, 2, pad = "0")

  dsm      <- terra::rast(file.path(dsm_dir, paste0(grid, ".tiff")))
  slope_r  <- terra::terrain(dsm, v = "slope",  unit = "degrees")
  aspect_r <- terra::terrain(dsm, v = "aspect", unit = "degrees")
  mean_elev <- terra::global(dsm, "mean", na.rm = TRUE)[1, 1]
  if (!is.finite(mean_elev)) mean_elev <- 0

  rgrass::initGRASS(gisBase = gisBase, home = tempdir(),
                    gisDbase = file.path(tempdir(), "grassdb"),
                    mapset = "PERMANENT", override = TRUE)
  rgrass::execGRASS("g.proj", flags = "c", epsg = 27700)

  tmp_dsm    <- file.path(tempdir(), "dsm.tif")
  tmp_slope  <- file.path(tempdir(), "slope.tif")
  tmp_aspect <- file.path(tempdir(), "aspect.tif")
  terra::writeRaster(dsm,      tmp_dsm,    overwrite = TRUE)
  terra::writeRaster(slope_r,  tmp_slope,  overwrite = TRUE)
  terra::writeRaster(aspect_r, tmp_aspect, overwrite = TRUE)

  rgrass::execGRASS("r.in.gdal", flags = c("o", "overwrite"), input = tmp_dsm,    output = "dsm")
  rgrass::execGRASS("r.in.gdal", flags = c("o", "overwrite"), input = tmp_slope,  output = "slope")
  rgrass::execGRASS("r.in.gdal", flags = c("o", "overwrite"), input = tmp_aspect, output = "aspect")
  rgrass::execGRASS("g.region", raster = "dsm")

  months <- 1:12
  rep_days <- as.Date(lubridate::ymd(paste0(year, "-", months, "-", day_of_month)))
  doy <- as.integer(format(rep_days, "%j"))
  days_in_month <- lubridate::days_in_month(rep_days)

  annual_from_monthly <- function(maps) {
    out <- maps[[1]] * days_in_month[1]
    for (i in 2:12) out <- out + maps[[i]] * days_in_month[i]
    out
  }

  # --- flat reference, cheap, FIRST ----------------------------------------
  # Same model, same days, no terrain: constant elevation, zero slope. Run on
  # a coarse region because with no relief there is nothing to resolve.
  #
  # Deliberately ahead of the terrain pass. This is seconds of work and the
  # terrain pass is the better part of an hour, so anything wrong here should
  # surface before that time is spent, not after it.
  #
  # Every argument goes inside `parameters`: rgrass rejects a call that mixes
  # GRASS options given as R arguments with a `parameters` list.
  rgrass::execGRASS("g.region",
                    parameters = list(raster = "dsm", res = as.character(flat_res)))
  rgrass::execGRASS("r.mapcalc", flags = "overwrite",
                    parameters = list(expression = sprintf("flat_elev = %.3f", mean_elev)))
  rgrass::execGRASS("r.mapcalc", flags = "overwrite",
                    parameters = list(expression = "flat_slope = 0"))
  rgrass::execGRASS("r.mapcalc", flags = "overwrite",
                    parameters = list(expression = "flat_aspect = 0"))

  flat_monthly <- numeric(12)
  for (i in seq_along(doy)) {
    out_map <- paste0("sun_flat_", sprintf("%02d", i))
    rgrass::execGRASS("r.sun", flags = "overwrite", parameters = list(
      elevation = "flat_elev", slope = "flat_slope", aspect = "flat_aspect",
      day = doy[i], step = 1, nprocs = nprocs, glob_rad = out_map))
    tif <- file.path(tempdir(), paste0(out_map, ".tif"))
    rgrass::execGRASS("r.out.gdal", flags = "overwrite", parameters = list(
      input = out_map, output = tif, format = "GTiff", type = "Float32",
      nodata = -9999))
    flat_monthly[i] <- terra::global(terra::rast(tif), "mean", na.rm = TRUE)[1, 1]
  }
  E_clear_flat_annual <- sum(flat_monthly * days_in_month)
  if (!is.finite(E_clear_flat_annual) || E_clear_flat_annual <= 0)
    stop("Flat clear-sky reference failed for ", grid, ": ", E_clear_flat_annual)

  # --- terrain clear sky, the expensive part -------------------------------
  rgrass::execGRASS("g.region", parameters = list(raster = "dsm"))
  clear_monthly <- vector("list", 12)
  for (i in seq_along(doy)) {
    out_map <- paste0("sun_clear_", sprintf("%02d", i))
    rgrass::execGRASS("r.sun", flags = "overwrite", parameters = list(
      elevation = "dsm", slope = "slope", aspect = "aspect",
      day = doy[i], step = 1, nprocs = nprocs, glob_rad = out_map))
    tif <- file.path(tempdir(), paste0(out_map, ".tif"))
    rgrass::execGRASS("r.out.gdal", flags = "overwrite", parameters = list(
      input = out_map, output = tif, format = "GTiff", type = "Float32",
      nodata = -9999))
    clear_monthly[[i]] <- terra::rast(tif)
  }
  E_clear_annual <- annual_from_monthly(clear_monthly)

  # --- clear-sky index and result -----------------------------------------
  csi <- E_ERA5_annual / E_clear_flat_annual
  if (csi <= 0 || csi > 1.2)
    warning(grid, ": clear-sky index ", round(csi, 3),
            " is outside the plausible range - check units")

  E_final_annual <- E_clear_annual * csi / 1000   # kWh/m2/year

  if (!is.null(clear_dir)) {
    terra::writeRaster(E_clear_annual,
                       file.path(clear_dir, paste0(grid, "_clearsky_annual_Whm2.tif")),
                       overwrite = TRUE, datatype = "FLT4S",
                       gdal = "COMPRESS=LZW", NAflag = -9999)
  }
  out_file <- file.path(out_dir, paste0(grid, "_annual_insolation_Whm2.tif"))
  terra::writeRaster(E_final_annual, filename = out_file, overwrite = TRUE,
                     datatype = "FLT4S", gdal = "COMPRESS=LZW", NAflag = -9999)

  invisible(data.frame(
    grid = grid,
    era5_annual = E_ERA5_annual,
    clear_flat_annual = E_clear_flat_annual,
    clear_terrain_mean = terra::global(E_clear_annual, "mean", na.rm = TRUE)[1, 1],
    csi = csi,
    final_mean = terra::global(E_final_annual, "mean", na.rm = TRUE)[1, 1],
    stringsAsFactors = FALSE))
}
