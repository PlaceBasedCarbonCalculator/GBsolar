# Build a smooth GB-wide field of annual all-sky irradiation from ERA5.
#
# Why this exists
#
# The ERA5 extracts in sampleData/ERA5/byGrid are one time series per OS 10 km
# square, sampled nearest-neighbour from a ~31 km reanalysis grid. Neighbouring
# squares therefore share identical annual totals and then step abruptly: among
# sixteen adjacent squares around Leeds and Sheffield there are only nine
# distinct annual values. Using those values directly as a per-tile multiplier
# would print the ERA5 cell boundaries onto the map.
#
# The underlying quantity - annual surface solar radiation downwards - varies
# smoothly over GB. So the block structure is a sampling artefact, not signal,
# and is removed here by fitting a smooth surface through the square centroids
# and evaluating it per square.

#' Decode an OS 10 km grid reference to its centre in EPSG:27700
#'
#' Two letters (500 km then 100 km square, "I" omitted from the sequence) then
#' two digits giving the 10 km square within it.
os_grid_centre <- function(ref) {
  L <- c(LETTERS[1:8], LETTERS[10:26])          # A-Z without I
  l1 <- match(substr(ref, 1, 1), L) - 1
  l2 <- match(substr(ref, 2, 2), L) - 1
  if (anyNA(c(l1, l2))) stop("Not an OS grid reference: ", ref)
  e100 <- ((l1 %% 5) * 5 + (l2 %% 5)) * 1e5 - 1e6
  n100 <- ((4 - l1 %/% 5) * 5 + (4 - l2 %/% 5)) * 1e5 - 5e5
  data.frame(
    grid = ref,
    easting  = e100 + as.integer(substr(ref, 3, 3)) * 1e4 + 5e3,
    northing = n100 + as.integer(substr(ref, 4, 4)) * 1e4 + 5e3,
    stringsAsFactors = FALSE
  )
}

#' Annual all-sky irradiation per grid square, Wh/m2/year
#'
#' SSRD is an accumulated flux in J/m2 over each hourly step, so the annual
#' total is the sum divided by 3600 - the same conversion the insolation code
#' uses, kept identical on purpose.
era5_annual_by_grid <- function(era5_dir = "sampleData/ERA5/byGrid") {
  files <- list.files(era5_dir, pattern = "\\.Rds$", full.names = TRUE)
  if (length(files) == 0) stop("No ERA5 .Rds files in ", era5_dir)
  out <- vapply(files, function(f) {
    e <- readRDS(f)
    sum(e$SSRD, na.rm = TRUE) / 3600
  }, numeric(1))
  data.frame(grid = sub("\\.Rds$", "", basename(files)),
             era5_raw = as.numeric(out), stringsAsFactors = FALSE)
}

#' Smooth the annual field over GB
#'
#' A thin-plate spline through the square centroids. `k` sets how much
#' structure survives: the point is to keep the real north-south and
#' coast-inland gradient while discarding steps that are narrower than ERA5
#' can actually resolve. k = 120 leaves features of roughly 60-80 km, which is
#' comfortably coarser than the ~31 km sampling and comfortably finer than the
#' national gradient.
era5_smooth_field <- function(annual, k = 120) {
  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("mgcv is needed to smooth the ERA5 field")
  xy <- do.call(rbind, lapply(annual$grid, os_grid_centre))
  d <- merge(annual, xy, by = "grid")
  d <- d[is.finite(d$era5_raw) & d$era5_raw > 0, ]
  fit <- mgcv::gam(era5_raw ~ s(easting, northing, k = k), data = d)
  d$era5_smooth <- as.numeric(predict(fit, newdata = d))
  d$residual_pct <- 100 * (d$era5_raw - d$era5_smooth) / d$era5_smooth
  attr(d, "fit") <- fit
  d
}

#' Build and save the lookup used by the insolation step
build_era5_field <- function(era5_dir = "sampleData/ERA5/byGrid",
                             out_csv = "data/era5_annual_smooth.csv",
                             k = 120) {
  annual <- era5_annual_by_grid(era5_dir)
  d <- era5_smooth_field(annual, k = k)
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  write.csv(d[, c("grid", "easting", "northing", "era5_raw", "era5_smooth",
                  "residual_pct")],
            out_csv, row.names = FALSE)
  message("Wrote ", nrow(d), " squares to ", out_csv)
  invisible(d)
}
