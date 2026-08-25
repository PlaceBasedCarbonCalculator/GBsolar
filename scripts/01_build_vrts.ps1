# Stage 1: mosaic -> reproject -> colour, entirely as VRTs.
#
# GBDEM materialised each of these stages as a real raster (a ~95 GB mosaic,
# then a reprojected copy, then an RGB copy) and its README notes the warp alone
# "will still take a day". None of that is necessary here. A GDAL VRT is a small
# XML description evaluated on read, so stages 1-3 cost 1.1 MB on disk and a few
# seconds, and the pixels are only ever computed once - when the tiler asks for
# them.
#
# Order matters: reproject FIRST (resampling insolation values), colour SECOND.
# Colouring first would make the warp average RGB triples instead of kWh/m2/yr,
# which is wrong wherever the ramp is non-linear in colour space.
#
# SOURCE: solarAnnualCSI, the clear-sky-index-corrected rasters, NOT solarAnnual.
# The two differ by a per-tile scalar spanning 0.62-1.05 (see
# R/apply_csi_correction.R), so a pyramid built from solarAnnual is wrong by up
# to 38% and cannot be patched after the fact - it has to be re-tiled.
#
# The work directory is deliberately separate from the old `Solar` one so the
# previously published GBsolar.pmtiles survives until this rebuild is verified.
#
# Prereq: RScripts/make_turbo_ramp.R has written $ramp into $work.

param(
    [string]$src  = "F:\DTM_DSM\GB_10k\solarAnnualCSI",
    [string]$work = "F:\DTM_DSM\large_rasters\SolarCSI",
    [string]$ramp = "turbo_100_1500.txt"
)

$gdal = "C:\OSGeo4W\bin"
New-Item -ItemType Directory -Force -Path $work | Out-Null
if (-not (Test-Path "$work\$ramp")) { throw "colour table $work\$ramp missing - run RScripts/make_turbo_ramp.R first" }

# gdalwarp refuses to write over an existing output and gdalbuildvrt would
# silently reuse a stale file list, so clear all three before rebuilding.
# They cost seconds to regenerate, so this is always the cheap option.
foreach ($v in "solar_27700.vrt", "solar_3857.vrt", "solar_rgba_3857.vrt") {
    Remove-Item "$work\$v" -ErrorAction SilentlyContinue
}

# --- 1. Mosaic the 2251 OS-grid 10 km tiles -----------------------------------
# -input_file_list because 2251 paths blow the command-line length limit.
# 2251, not 2253: NS01 and NS05 are entirely sea and have no corrected output.
Get-ChildItem $src -Filter *.tif | Select-Object -ExpandProperty FullName |
    Set-Content "$work\tif_list.txt" -Encoding ascii
& "$gdal\gdalbuildvrt.exe" -input_file_list "$work\tif_list.txt" "$work\solar_27700.vrt"

# --- 2. Reproject 27700 -> 3857 ----------------------------------------------
# -tr is the exact z14 resolution for 512 px tiles:
#     156543.03392804097 / 2^14 / 2 = 4.77731426716 projected m/px
# Matching the tile grid resolution here means the tiler does no further
# rescaling of the base tiles.
# -r average because this is a downsample (2 m native -> ~4.78 projected m);
# 'near' would alias badly on a surface this noisy.
# nodata is carried through as -9999 so the colour table's "nv" entry can make
# it transparent rather than rendering sea as the bottom of the ramp.
$res = 4.77731426716
& "$gdal\gdalwarp.exe" -of VRT -s_srs EPSG:27700 -t_srs EPSG:3857 -r average `
    -tr $res $res -srcnodata -9999 -dstnodata -9999 -multi `
    "$work\solar_27700.vrt" "$work\solar_3857.vrt"

# --- 3. Apply the Turbo colour table -----------------------------------------
# -of VRT keeps this virtual too; -alpha adds the 4th band so "nv 0 0 0 0"
# yields genuine transparency.
& "$gdal\gdaldem.exe" color-relief "$work\solar_3857.vrt" "$work\$ramp" `
    "$work\solar_rgba_3857.vrt" -of VRT -alpha

& "$gdal\gdalinfo.exe" -nofl "$work\solar_rgba_3857.vrt" |
    Select-String "Size is|Pixel Size|Band |ColorInterp"
