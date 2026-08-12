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
# -z 5-14          z14 at 512 px is 4.78 projected m/px, just coarser than the
#                  2 m source at GB latitudes, so z14 is the natural floor.
#                  z5 shows the whole island.
# -r average       for building the overview zooms.
# -w none          no HTML viewer boilerplate.

$work = "F:\DTM_DSM\large_rasters\Solar"

cmd /c "call C:\OSGeo4W\bin\o4w_env.bat && python -m osgeo_utils.gdal2tiles " +
       "--xyz --tilesize=512 -z 5-14 --processes=14 -r average --no-kml -w none " +
       "$work\solar_rgba_3857.vrt $work\tiles_png"
