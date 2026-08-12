# Generate the Turbo colour ramp used to render the GB annual-insolation raster.
#
# Produces two files in F:/DTM_DSM/large_rasters/Solar:
#   turbo_0_2000.txt  - the gdaldem color-relief colour table
#   turbo_0_2000.csv  - the same mapping as data, for the documentation table
#
# Why this and not rio-rgbify: rio-rgbify encodes a number into RGB in pseudo
# base-256 so the original value can be decoded client-side. That is right for a
# DEM, where MapLibre reads the elevation back out. Here we want a *colour ramp*
# - the pixel's colour is the message - so the encoding step is a colour table
# applied with `gdaldem color-relief`, and the mapping below is the definition
# of what any colour on the map means.
#
# Domain: 0 -> 2000 Wh/m2, as specified. gdaldem clamps outside the table, so
# anything above 2000 renders as the final dark red rather than being lost.
# Note the observed data minimum is about 97-130 Wh/m2, so the darkest blues at
# the very bottom of the ramp are in practice unused; the domain starts at 0
# because a fixed, round domain keeps the legend stable and comparable between
# renders.

out_dir <- "F:/DTM_DSM/large_rasters/Solar"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

vmin <- 0
vmax <- 2000
n    <- 65   # anchors; gdaldem interpolates linearly between them

cols <- viridisLite::turbo(n)              # exact Google Turbo
rgb  <- t(grDevices::col2rgb(cols))        # n x 3, values 0-255
vals <- round(seq(vmin, vmax, length.out = n), 3)

tbl <- data.frame(value_Whm2 = vals,
                  r = rgb[, 1], g = rgb[, 2], b = rgb[, 3],
                  hex = substr(cols, 1, 7))

# gdaldem colour table. "nv" is the nodata entry (-9999 in these rasters) and is
# fully transparent, so sea and unmapped ground drop out rather than rendering
# as the bottom colour of the ramp.
lines <- c(sprintf("%s %d %d %d 255", format(tbl$value_Whm2, trim = TRUE, scientific = FALSE),
                   tbl$r, tbl$g, tbl$b),
           "nv 0 0 0 0")
writeLines(lines, file.path(out_dir, "turbo_0_2000.txt"))
write.csv(tbl, file.path(out_dir, "turbo_0_2000.csv"), row.names = FALSE)

cat("wrote", n, "anchors spanning", vmin, "-", vmax, "Wh/m2\n")
print(utils::head(tbl, 3))
print(utils::tail(tbl, 3))
