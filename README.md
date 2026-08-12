# GBsolar
A spin off of GB DEM to create a national solar map

## Status: PMTiles layer built (2026-08-09)

`F:\DTM_DSM\large_rasters\Solar\GBsolar.pmtiles` — 5.81 GB, WebP, zoom 5–14,
Turbo colour ramp linear over 0–2000 Wh/m², nodata transparent.

**See [METHOD.md](METHOD.md)** for the full method, the colour↔value table, what
was done differently from GBDEM and why, measured outcomes, and the one failure
encountered. Scripts are in [`scripts/`](scripts/) and
[`RScripts/`](RScripts/), in the order they run.

Two things to decide before publishing:

1. The measured data maximum is **2288 Wh/m²**, not ~2000, so with the specified
   domain dark red means "2000 or above" rather than the true maximum. The data
   minimum is 98, so the darkest blues never appear and the render uses roughly
   the middle half of Turbo. Both are a one-line change in
   `RScripts/make_turbo_ramp.R`.
2. Tiles are **lossy** WebP (q85), unlike GBDEM's lossless. Safe here because the
   colour is only ever looked at, never decoded, and it is the difference
   between ~5.8 GB and ~35 GB. Measured cost is ~4.5/255 of colour error.

Add it to MapLibre as `type: "raster"` — **not** `raster-dem`, which would try to
decode the colours as terrain heights.

# Instructions for Claude

This task is to only be attempted if you have completed the all the jobs of the PublicTransportAnalysis and build repos.

This repo "GBsolar" is a addtion to the Carbon and Place analysis by creating a high resolution (2m) map of solar power potential in Great Britain. Some work has already been done by running the /RScripts/built_isolation_tiles_10k.R script. Which has made a folder 10km x 10km insolation maps saved in F:\DTM_DSM\GB_10k\solarAnnual.

This map builds of the ../GBDEM repo which created the digital terrain map, which was a key input to the analysis.

These itf files are aligned with the British National Grid tiles (https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid).

Your task is to go from a folder of 10km rasters tifs to a pmtiles file of coloured raster images that show the levels of insolation across Great Britain.

Fortunately there are some guides on how to do this in the ../GBDEM/README.md

Key steps are:

1. Mosaic the 10km rasters into a single large raster of whole GB 
2. Reproject the raster from epsg:27700 to epsg:3857
3. Convert from single band to RGB encoded raster
4. Generate XYZ Tiles (directory)
5. Convert png images to webp using imagemagic
6. Convert to MBTiles
7. Convert to PMtiles

The pipeline outlines in GBDEM is a bit convoluted but was discover through trial an error. Due to the massive size of the raster some alterative methods simple crash or run forever. 

Bring these instruction from GBDEM over to GBsolar as the starting point for your work. But be prepared to encounter some problems. For example the first stage GB_10km_mosaic.R  has references to the code crashing and has code commented out. It is possible that these are the wrong instructions but the actual method successfully used has been lost. So treat GBDEM as a rough guide not perfect code.

Another challenge for you is the conversion from a single band of insolation values in Wh/m2 to colours. rio-rgbify was designed to encode arbitrary bit depth rasters in pseudo base-256 as RGB which works for a DEM raster. But for the solar raster we want a colour ramp, use Turbo from QGIS (Dark Blue - Blue - Green - Yellow - Orange - Red - Dark Red see https://research.google/blog/turbo-an-improved-rainbow-colormap-for-visualization/ for details). Where dark blue represents 0 and dark red is the maximum value in any of the rasters ~2000 Wh/m2. It is important you document the colour to value mapping. rio-rgbify may be the wrong tool for this job. 

The end goal is for the pmtiles file to be added to the Carbon and Place map a a colour raster layer clearly showing solar hot spots. 

Do as many of these tasks as you can. And write clear documentation about what you did what you tried and the outcomes. Create new code and documentation in this repo down change existing work in other repos. If you create intermediate rasters/files keep them saved in F:\DTM_DSM\large_rasters\Solar however. The F drive is a HDD and slow. If you need faster disk speed create a tempSolar folder on the C drive which is an SSD and faster. However it has limited free space ~ 200GB so you may have to move files to F drive after creation.
