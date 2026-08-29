# Stage 2: XYZ tile pyramid from the RGBA colour VRT.
#
# GBDEM drove this from the QGIS processing toolbox (qgis:tilesxyzdirectory).
# gdal2tiles is used here instead: it is scriptable without a QGIS session, it
# reads the VRT chain directly, and --tilesize removes the reason GBDEM needed
# QGIS in the first place (its note that "the mbtiles version can't create
# high-resolution 512x512 images").
#
# gdal2tiles is not on PATH in this OSGeo4W install and its .bat wrapper does not
# resolve, so it is invoked as a module after sourcing o4w_env.bat.
#
# --xyz            OSM/slippy numbering, which is what MapLibre expects.
#                  Without it you get TMS and the map is upside down.
# --tilesize=512   matches GBDEM's output and halves the tile count.
# -z 5-$maxZoom    must match the resolution baked into the VRT chain by
#                  scripts/01_build_vrts.ps1's $maxZoom - see the note there.
#                  z5 shows the whole island at the low end regardless.
# -r average       for building the overview zooms.
# -w none          no HTML viewer boilerplate.
#
# gdal2tiles' base zoom is where nearly all the time goes - it is the one level
# resampled directly from the source VRT, and every coarser zoom is then built
# cheaply by 2x2-downsampling the level below. Going from a z14 to a z15 base
# roughly quadruples the base tile count (404x824 -> ~808x1648), so this step
# roughly quadruples too: the z5-14 pyramid took 10h08m on 2026-08-25/26,
# so budget on the order of a day and a half for z5-15, not "a bit longer".
#
# gdal2tiles cannot resume, so run it detached and keep the log:
#   $p = Start-Process powershell -PassThru -WindowStyle Hidden -ArgumentList `
#          '-File','scripts\02_tiles.ps1','-RedirectStandardOutput',"$work\tiling.log"

param(
    [string]$work    = "F:\DTM_DSM\large_rasters\SolarCSI",
    [int]   $maxZoom = 15
)

$vrt = "$work\solar_rgba_3857.vrt"
$out = "$work\tiles_png"
if (-not (Test-Path $vrt)) { throw "$vrt missing - run scripts/01_build_vrts.ps1 first" }

# A part-written pyramid from an interrupted run would be silently merged into
# this one, so refuse rather than guess. Delete it deliberately to restart.
if (Test-Path $out) { throw "$out already exists - remove it before re-tiling" }

# Built as one string and handed to cmd as a single argument. The previous
# version relied on `cmd /c "..." + "..."`, which PowerShell parses in argument
# mode: `+` is passed through as a literal argument rather than concatenating,
# so the command cmd actually received was not the one written here.
$cmd = 'call C:\OSGeo4W\bin\o4w_env.bat && ' +
       'python -m osgeo_utils.gdal2tiles ' +
       "--xyz --tilesize=512 -z 5-$maxZoom --processes=14 -r average --no-kml -w none " +
       "`"$vrt`" `"$out`""

Write-Output "running: $cmd"
$t0 = Get-Date
cmd /c $cmd
Write-Output ("gdal2tiles finished in {0:n2} hours" -f ((Get-Date) - $t0).TotalHours)
