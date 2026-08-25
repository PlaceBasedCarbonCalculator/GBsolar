# Generate the Turbo colour ramp used to render the GB annual-insolation raster.
#
# Produces two files in the working directory below:
#   turbo_100_1500.txt  - the gdaldem color-relief colour table
#   turbo_100_1500.csv  - the same mapping as data, for the documentation table
#
# Why this and not rio-rgbify: rio-rgbify encodes a number into RGB in pseudo
# base-256 so the original value can be decoded client-side. That is right for a
# DEM, where MapLibre reads the elevation back out. Here we want a *colour ramp*
# - the pixel's colour is the message - so the encoding step is a colour table
# applied with `gdaldem color-relief`, and the mapping below is the definition
# of what any colour on the map means.
#
# DOMAIN: 100 -> 1500 kWh/m2/year, linear.
#
# The first published version of this map used 0 -> 2000, taken from a brief
# asking that dark red be "the maximum value in any of the rasters ~2000". That
# domain no longer fits the data. Measured over the CSI-corrected mosaic
# (`gdalinfo -approx_stats`, 2026-08-25):
#
#              min     mean      sd     max
#   old       98.2   1074.2   271.9  2288.0
#   new       87.7    965.2   238.9  1743.2
#
# so on the corrected rasters a 0-2000 domain leaves the top 13% of Turbo
# unreachable - dark red would be a colour no pixel ever takes - and squeezes
# the central 90% of pixels into 41% of the ramp.
#
# 100-1500 was chosen against the measured percentiles instead. It costs very
# little clipping and buys a lot of contrast:
#
#   pixels above 1500 (clamped to dark red)   ~0.1%
#   pixels below 100  (clamped to dark blue)  ~0.009%
#   central 90% of pixels now span            58% of the ramp (was 41%)
#
# Both ends are therefore clamps rather than extremes: `#7A0403` means
# "1500 or above" and `#30123B` means "100 or below". That is a deliberate
# trade for hot-spot separation, and it is why the legend must be labelled with
# the domain rather than with the data range. Change `vmin`/`vmax` below to
# revisit it, and rename the output files to match - the domain is in the
# filename on purpose, so a tile pyramid can always be traced to the ramp that
# coloured it.

out_dir <- "F:/DTM_DSM/large_rasters/SolarCSI"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

vmin <- 100
vmax <- 1500
n    <- 65   # anchors; gdaldem interpolates linearly between them

stem <- sprintf("turbo_%g_%g", vmin, vmax)

cols <- viridisLite::turbo(n)              # exact Google Turbo
rgb  <- t(grDevices::col2rgb(cols))        # n x 3, values 0-255
vals <- round(seq(vmin, vmax, length.out = n), 3)

tbl <- data.frame(value_kWhm2 = vals,
                  r = rgb[, 1], g = rgb[, 2], b = rgb[, 3],
                  hex = substr(cols, 1, 7))

# gdaldem colour table. "nv" is the nodata entry (-9999 in these rasters) and is
# fully transparent, so sea and unmapped ground drop out rather than rendering
# as the bottom colour of the ramp.
lines <- c(sprintf("%s %d %d %d 255", format(tbl$value_kWhm2, trim = TRUE, scientific = FALSE),
                   tbl$r, tbl$g, tbl$b),
           "nv 0 0 0 0")
writeLines(lines, file.path(out_dir, paste0(stem, ".txt")))
write.csv(tbl, file.path(out_dir, paste0(stem, ".csv")), row.names = FALSE)

cat("wrote", n, "anchors spanning", vmin, "-", vmax, "kWh/m2/year to", stem, "\n")
print(utils::head(tbl, 3))
print(utils::tail(tbl, 3))
