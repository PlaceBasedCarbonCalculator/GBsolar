# Correct the existing insolation rasters in place of a full re-run.
#
# THE ALGEBRA
#
# The old code produced, for tile g,
#
#     E_current(x) = E_clear(x) * s_g,     s_g = E_ERA5_raw_g / (1000 * mean_clear_g)
#
# and what is wanted is
#
#     E_correct(x) = E_clear(x) * CSI_g / 1000,   CSI_g = E_ERA5_smooth_g / E_clear_flat_g
#
# Dividing one by the other, E_clear(x) cancels completely:
#
#     E_correct(x) = E_current(x) * T_g * R_g
#       T_g = mean_clear_g / E_clear_flat_g      the tile's mean terrain factor
#       R_g = E_ERA5_smooth_g / E_ERA5_raw_g     de-blocking the reanalysis
#
# So the defect is entirely a per-tile scalar, and the 2 m within-tile pattern
# in the existing rasters - the expensive, correct part - can be kept as it is.
# Only the level of each tile was wrong.
#
# R_g is free. T_g is the one thing that needs r.sun. But T_g is a spatial mean
# over ~25 million pixels, not a per-pixel value, and spatial means survive
# downsampling far better than the fields they come from. Computing it on a
# coarsened DSM is therefore worth testing against the full-resolution answer:
# if the two agree, the whole correction costs hours rather than the ~27 days
# a full 2 m re-run of all 2,253 tiles would take.
#
# `terrain_factor()` computes T_g at whatever resolution it is given.
# `calibrate_terrain_factor()` is what decides whether the cheap version is
# good enough - do not skip it.

#' Check an existing tile really was produced by the old normalisation
#'
#' If it was, its spatial mean equals its own raw ERA5 annual total over 1000,
#' by construction. A tile that fails this was made some other way and must not
#' be corrected by rescaling - it has to be rebuilt.
verify_old_normalisation <- function(grids, era5_field,
                                     in_dir = "F:/DTM_DSM/GB_10k/solarAnnual",
                                     tol_pct = 0.01) {
  do.call(rbind, lapply(grids, function(g) {
    f <- file.path(in_dir, paste0(g, "_annual_insolation_Whm2.tif"))
    if (!file.exists(f))
      return(data.frame(grid = g, mean_actual = NA_real_, mean_expected = NA_real_,
                        diff_pct = NA_real_, ok = FALSE, stringsAsFactors = FALSE))
    m <- terra::global(terra::rast(f), "mean", na.rm = TRUE)[1, 1]
    i <- match(g, era5_field$grid)
    expect <- if (is.na(i)) NA_real_ else era5_field$era5_raw[i] / 1000
    d <- 100 * (m - expect) / expect
    data.frame(grid = g, mean_actual = m, mean_expected = expect,
               diff_pct = d, ok = is.finite(d) && abs(d) < tol_pct,
               stringsAsFactors = FALSE)
  }))
}

#' Mean terrain factor for one tile: mean clear-sky over terrain, divided by
#' clear-sky on flat ground at the same place
#'
#' @param res cell size in metres. NULL uses the DSM's native resolution.
#' @param months which months to use as representative days. T_g is a ratio of
#'   two annual sums built the same way, so much of the seasonal detail may
#'   cancel; using fewer days is only legitimate if the ratio is shown to be
#'   unchanged, which is what calibrating against months = 1:12 is for. Each
#'   chosen month is weighted by the days of every month nearest to it, so the
#'   weights always total 365/366.
terrain_factor <- function(grid,
                           res = NULL,
                           months = 1:12,
                           day_of_month = 15,
                           year = 2020,
                           dsm_dir = "F:/DTM_DSM/GB_10k/DSM",
                           gisBase = "C:/Program Files/GRASS GIS 8.4",
                           nprocs = 35,
                           flat_res = 500) {

  f_dsm <- file.path(dsm_dir, paste0(grid, ".tiff"))
  stopifnot(file.exists(f_dsm))
  day_of_month <- stringr::str_pad(day_of_month, 2, pad = "0")

  dsm <- terra::rast(f_dsm)
  if (!is.null(res)) {
    fact <- max(1, round(res / terra::res(dsm)[1]))
    if (fact > 1) dsm <- terra::aggregate(dsm, fact = fact, fun = "mean", na.rm = TRUE)
  }
  slope_r   <- terra::terrain(dsm, v = "slope",  unit = "degrees")
  aspect_r  <- terra::terrain(dsm, v = "aspect", unit = "degrees")
  mean_elev <- terra::global(dsm, "mean", na.rm = TRUE)[1, 1]
  if (!is.finite(mean_elev)) mean_elev <- 0

  rgrass::initGRASS(gisBase = gisBase, home = tempdir(),
                    gisDbase = file.path(tempdir(), "grassdb"),
                    mapset = "PERMANENT", override = TRUE)
  rgrass::execGRASS("g.proj", flags = "c", epsg = 27700)

  wr <- function(r, nm) {
    p <- file.path(tempdir(), paste0(nm, ".tif"))
    terra::writeRaster(r, p, overwrite = TRUE)
    rgrass::execGRASS("r.in.gdal", flags = c("o", "overwrite"), input = p, output = nm)
  }
  wr(dsm, "dsm"); wr(slope_r, "slope"); wr(aspect_r, "aspect")
  rgrass::execGRASS("g.region", raster = "dsm")

  rep_days <- as.Date(lubridate::ymd(paste0(year, "-", months, "-", day_of_month)))
  doy <- as.integer(format(rep_days, "%j"))

  # Weight each representative month by the days of every month that is closer
  # to it than to any other, wrapping December round to January so winter is
  # not split. With months = 1:12 this is just days_in_month.
  all_days <- as.numeric(lubridate::days_in_month(
    as.Date(lubridate::ymd(paste0(year, "-", 1:12, "-", day_of_month)))))
  circ <- function(a, b) pmin(abs(a - b), 12 - abs(a - b))
  owner <- vapply(1:12, function(m) months[which.min(circ(m, months))], numeric(1))
  dim_ <- vapply(months, function(m) sum(all_days[owner == m]), numeric(1))
  stopifnot(abs(sum(dim_) - sum(all_days)) < 1e-9)

  run_mean <- function(elev, slope, aspect, tag) {
    v <- numeric(12)
    for (i in seq_along(doy)) {
      m <- paste0(tag, sprintf("%02d", i))
      rgrass::execGRASS("r.sun", flags = "overwrite", parameters = list(
        elevation = elev, slope = slope, aspect = aspect,
        day = doy[i], step = 1, nprocs = nprocs, glob_rad = m))
      tif <- file.path(tempdir(), paste0(m, ".tif"))
      rgrass::execGRASS("r.out.gdal", flags = "overwrite", parameters = list(
        input = m, output = tif, format = "GTiff", type = "Float32",
        nodata = -9999))
      v[i] <- terra::global(terra::rast(tif), "mean", na.rm = TRUE)[1, 1]
    }
    sum(v * dim_)
  }

  # Flat reference first: it is seconds of work and the terrain pass is not,
  # so a mistake here should not cost the terrain pass. Every argument goes
  # inside `parameters` - rgrass rejects a call that mixes GRASS options given
  # as R arguments with a `parameters` list.
  rgrass::execGRASS("g.region",
                    parameters = list(raster = "dsm", res = as.character(flat_res)))
  rgrass::execGRASS("r.mapcalc", flags = "overwrite",
                    parameters = list(expression = sprintf("flat_elev = %.3f", mean_elev)))
  rgrass::execGRASS("r.mapcalc", flags = "overwrite",
                    parameters = list(expression = "flat_slope = 0"))
  rgrass::execGRASS("r.mapcalc", flags = "overwrite",
                    parameters = list(expression = "flat_aspect = 0"))
  clear_flat <- run_mean("flat_elev", "flat_slope", "flat_aspect", "tf_flat_")

  rgrass::execGRASS("g.region", parameters = list(raster = "dsm"))
  mean_clear <- run_mean("dsm", "slope", "aspect", "tf_terr_")

  data.frame(grid = grid, res = if (is.null(res)) terra::res(dsm)[1] else res,
             mean_clear = mean_clear, clear_flat = clear_flat,
             terrain_factor = mean_clear / clear_flat, stringsAsFactors = FALSE)
}

#' Apply the per-tile correction to an existing raster
apply_correction <- function(grid, terrain_factor_g, era5_field,
                             in_dir  = "F:/DTM_DSM/GB_10k/solarAnnual",
                             out_dir = "F:/DTM_DSM/GB_10k/solarAnnualCSI") {
  i <- match(grid, era5_field$grid)
  if (is.na(i)) stop("No ERA5 value for ", grid)
  R_g <- era5_field$era5_smooth[i] / era5_field$era5_raw[i]
  mult <- terrain_factor_g * R_g

  r <- terra::rast(file.path(in_dir, paste0(grid, "_annual_insolation_Whm2.tif")))
  out <- r * mult
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  terra::writeRaster(out, file.path(out_dir, paste0(grid, "_annual_insolation_Whm2.tif")),
                     overwrite = TRUE, datatype = "FLT4S",
                     gdal = "COMPRESS=LZW", NAflag = -9999)
  invisible(data.frame(grid = grid, terrain_factor = terrain_factor_g,
                       era5_ratio = R_g, multiplier = mult,
                       stringsAsFactors = FALSE))
}
