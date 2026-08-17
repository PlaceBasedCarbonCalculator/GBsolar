#!/usr/bin/env Rscript
#
# Correct every 10 km insolation tile from the per-tile-mean normalisation to
# a clear-sky index. See R/apply_csi_correction.R for the algebra and
# R/insolation_calcs_csi.R for why the old normalisation was wrong.
#
# For each tile this computes the true 2 m terrain factor
#   T_g = mean(clear-sky over terrain) / clear-sky on flat ground
# with r.sun, then writes E_current * T_g * R_g to solarAnnualCSI.
#
# Four representative days rather than twelve. T_g is a ratio of two annual
# sums built identically, so the seasonal detail cancels; measured against the
# full twelve-day answer on three deliberately unalike tiles:
#
#   SE33  urban Leeds      0.83756 -> 0.83723   -0.039%
#   NN17  Ben Nevis        0.82549 -> 0.82539   -0.012%
#   TF20  Lincolnshire fen 0.96581 -> 0.96558   -0.024%
#
# and it cuts a tile from 20.3 to 6.7 minutes. Note those three span 0.825 to
# 0.966: that 14% spread is the real between-tile variation the old code was
# erasing, and is why this is worth doing.
#
# The DSM must stay at its native 2 m. Estimating T_g from a coarsened DSM was
# tried and fails badly - 0.838 becomes 0.966 at 10 m and 0.996 at 50 m,
# because most of Leeds' shading is buildings and averaging deletes them. That
# is a larger error than the defect being fixed.
#
# RESUMABLE. Every tile appends a row to the progress CSV as it finishes, and
# a re-run skips whatever is already recorded ok. At roughly 6.7 min a tile
# this takes about ten days, so it will almost certainly be interrupted at
# some point; just start it again. A tile that errors is recorded and skipped
# rather than stopping the run.
#
#   Rscript RScripts/run_csi_correction.R           # run, resuming
#   Rscript RScripts/run_csi_correction.R --status  # progress only, run nothing

suppressPackageStartupMessages(library(terra))
source("R/apply_csi_correction.R")

IN_DIR   <- "F:/DTM_DSM/GB_10k/solarAnnual"
OUT_DIR  <- "F:/DTM_DSM/GB_10k/solarAnnualCSI"
DSM_DIR  <- "F:/DTM_DSM/GB_10k/DSM"
PROGRESS <- "F:/DTM_DSM/GB_10k/csi_correction_progress.csv"
MONTHS   <- c(1, 4, 7, 10)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
field <- read.csv("data/era5_annual_smooth.csv", stringsAsFactors = FALSE)

grids <- sub("_annual_insolation_Whm2\\.tif$", "",
             list.files(IN_DIR, pattern = "\\.tif$"))

read_progress <- function() {
  if (!file.exists(PROGRESS))
    return(data.frame(grid = character(0), status = character(0),
                      stringsAsFactors = FALSE))
  utils::read.csv(PROGRESS, stringsAsFactors = FALSE)
}

done <- read_progress()
# "empty" is as final as "ok" - an all-sea square has nothing to correct and
# will be empty again next time. Only errors are worth retrying on a restart.
todo <- setdiff(grids, done$grid[done$status %in% c("ok", "empty")])

if ("--status" %in% commandArgs(trailingOnly = TRUE)) {
  cat("tiles total     :", length(grids), "\n")
  cat("done ok         :", sum(done$status == "ok"), "\n")
  cat("skipped (empty) :", sum(done$status == "empty"), "\n")
  cat("errored         :", sum(done$status == "error"), "\n")
  cat("remaining       :", length(todo), "\n")
  if (nrow(done) && any(done$status == "ok")) {
    m <- mean(done$minutes[done$status == "ok"], na.rm = TRUE)
    cat("mean min/tile   :", round(m, 2), "\n")
    cat("est. remaining  :", round(length(todo) * m / 60 / 24, 2), "days\n")
  }
  quit(save = "no", status = 0)
}

append_row <- function(row) {
  utils::write.table(row, PROGRESS, sep = ",", row.names = FALSE,
                     col.names = !file.exists(PROGRESS), append = file.exists(PROGRESS))
}

message("tiles to process: ", length(todo), " of ", length(grids),
        " (", length(grids) - length(todo), " already done)")
t_start <- Sys.time()

for (i in seq_along(todo)) {
  g <- todo[i]
  t0 <- Sys.time()

  row <- tryCatch({
    # All-sea squares have an entirely empty DSM - NS01 and NS05 are the two in
    # GB. There is nothing to scale and r.sun has nothing to work on.
    dsm_f <- file.path(DSM_DIR, paste0(g, ".tiff"))
    if (!file.exists(dsm_f)) stop("no DSM")
    if (global(!is.na(rast(dsm_f)), "sum", na.rm = TRUE)[1, 1] == 0) {
      data.frame(grid = g, terrain_factor = NA_real_, era5_ratio = NA_real_,
                 multiplier = NA_real_, minutes = 0, status = "empty",
                 finished = format(Sys.time()), stringsAsFactors = FALSE)
    } else {
      tf <- terrain_factor(g, res = NULL, months = MONTHS, dsm_dir = DSM_DIR)
      ap <- apply_correction(g, tf$terrain_factor, field,
                             in_dir = IN_DIR, out_dir = OUT_DIR)
      data.frame(grid = g, terrain_factor = tf$terrain_factor,
                 era5_ratio = ap$era5_ratio, multiplier = ap$multiplier,
                 minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
                 status = "ok", finished = format(Sys.time()),
                 stringsAsFactors = FALSE)
    }
  }, error = function(e) {
    message("  ERROR on ", g, ": ", conditionMessage(e))
    data.frame(grid = g, terrain_factor = NA_real_, era5_ratio = NA_real_,
               multiplier = NA_real_,
               minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
               status = "error", finished = format(Sys.time()),
               stringsAsFactors = FALSE)
  })

  append_row(row)

  el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  rate <- el / i
  message(sprintf("[%d/%d] %s  T=%s  x%s  %.1f min | elapsed %.1f h | eta %.1f days",
                  i, length(todo), g,
                  ifelse(is.na(row$terrain_factor), "-", sprintf("%.4f", row$terrain_factor)),
                  ifelse(is.na(row$multiplier), "-", sprintf("%.4f", row$multiplier)),
                  row$minutes, el / 60, rate * (length(todo) - i) / 60 / 24))
}

message("finished in ", round(as.numeric(difftime(Sys.time(), t_start, units = "hours")), 1), " hours")
