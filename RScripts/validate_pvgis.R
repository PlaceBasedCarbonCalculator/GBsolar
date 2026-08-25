#!/usr/bin/env Rscript
#
# Validate the absolute level of the corrected rasters against PVGIS.
#
#   Rscript RScripts/validate_pvgis.R
#
# Writes data/validation_pvgis.csv. Takes a few minutes: it is polite to the
# PVGIS API and sleeps between calls. See METHOD.md section 6 for the results
# and what the latitude trend in the residuals means.
#
# mean(corrected_g) = T_g * E_ERA5_smooth_g, so mean/T_g recovers the modelled
# flat-ground annual GHI for the tile. PVGIS's radiation database is SARAH2
# (Meteosat-derived), independent of the ERA5 field used here, so comparing the
# two tests the level the map is published at.
#
# The query point is the centroid of the tile's LAND pixels, not of the tile.
# Sea is NA in these rasters, so the tile mean already describes land only; a
# geometric centroid falls in the water for coastal and island squares and
# PVGIS rejects it ("Location over the sea"). Where even the land centroid is
# rejected - a horseshoe-shaped coastline - it falls back to actual sampled
# land pixels, nearest the centroid first.
suppressPackageStartupMessages({library(terra); library(jsonlite)})
terraOptions(progress = 0)

DIR <- "F:/DTM_DSM/GB_10k/solarAnnualCSI"
p <- read.csv("data/csi_correction_progress.csv", stringsAsFactors = FALSE)
p <- p[p$status == "ok", ]; p <- p[!duplicated(p$grid), ]

hdr <- do.call(rbind, lapply(p$grid, function(g) {
  e <- ext(rast(file.path(DIR, paste0(g, "_annual_insolation_Whm2.tif"))))
  data.frame(grid = g, y = mean(e[3:4]))
}))
p <- merge(p, hdr, by = "grid")

set.seed(7)
p$band <- cut(p$y, breaks = 12, labels = FALSE)
sel <- do.call(rbind, lapply(split(p, p$band), function(d) d[sample(nrow(d), min(2, nrow(d))), ]))

to_ll <- function(x, y) {
  g <- geom(project(vect(data.frame(x = x, y = y), geom = c("x", "y"),
                         crs = "EPSG:27700"), "EPSG:4326"))
  data.frame(lon = g[, "x"], lat = g[, "y"])
}

pvgis <- function(lat, lon) {
  u <- sprintf(paste0("https://re.jrc.ec.europa.eu/api/v5_2/MRcalc?lat=%.5f&lon=%.5f",
                      "&horirrad=1&startyear=2020&endyear=2020&outputformat=json"), lat, lon)
  r <- tryCatch(suppressWarnings(fromJSON(u)), error = function(e) NULL)
  if (is.null(r) || is.null(r$outputs$monthly)) return(NA_real_)
  sum(r$outputs$monthly$`H(h)_m`)
}

out <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i) {
  g <- sel$grid[i]
  r <- rast(file.path(DIR, paste0(g, "_annual_insolation_Whm2.tif")))
  m <- global(r, "mean", na.rm = TRUE)[1, 1]
  flat <- m / sel$terrain_factor[i]

  s <- spatSample(r, 4000, method = "regular", na.rm = TRUE, xy = TRUE)
  # Land centroid first, then real land pixels ordered by distance from it.
  cand <- rbind(data.frame(x = mean(s$x), y = mean(s$y)),
                s[order((s$x - mean(s$x))^2 + (s$y - mean(s$y))^2), c("x", "y")][1:40, ])
  ll <- to_ll(cand$x, cand$y)

  gh <- NA_real_; k <- 0
  for (j in seq_len(nrow(ll))) {
    gh <- pvgis(ll$lat[j], ll$lon[j]); k <- j
    Sys.sleep(0.6)
    if (is.finite(gh)) break
  }
  data.frame(grid = g, lat = round(ll$lat[k], 3), lon = round(ll$lon[k], 3),
             tries = k, land_pct = round(100 * mean(!is.na(values(r))), 1),
             T_g = round(sel$terrain_factor[i], 4),
             tile_mean = round(m, 1), model_flat = round(flat, 1),
             pvgis = round(gh, 1), diff_pct = round(100 * (flat - gh) / gh, 2))
}))

out <- out[order(out$lat), ]
print(out, row.names = FALSE)
d <- out$diff_pct[is.finite(out$diff_pct)]
cat("\nn =", length(d), "of", nrow(out), "\n")
cat(sprintf("mean bias   : %+.2f%%\n", mean(d)))
cat(sprintf("median bias : %+.2f%%\n", median(d)))
cat(sprintf("sd          :  %.2f%%\n", sd(d)))
cat(sprintf("range       : %+.2f%% to %+.2f%%\n", min(d), max(d)))
ok <- out[is.finite(out$diff_pct), ]
cat(sprintf("correlation : %.4f\n", cor(ok$model_flat, ok$pvgis)))
write.csv(out, "data/validation_pvgis.csv", row.names = FALSE)
cat("\nwrote data/validation_pvgis.csv\n")
