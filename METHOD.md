# From 10 km insolation rasters to a GB PMTiles colour map

Method, colour definition and outcomes for turning
`F:\DTM_DSM\GB_10k\solarAnnualCSI` (2,251 OS-grid 10 km GeoTIFFs of annual
insolation, 170 GB) into a single PMTiles raster layer for carbon.place.

Written 2026-08-09, rebuilt 2026-08-26 and again from 2026-08-29. Scripts
referenced live in `scripts/` and `RScripts/`.

> **The source changed on 2026-08-25.** The first build read `solarAnnual`,
> whose per-tile level was set by a normalisation later found to be wrong (see
> [`R/apply_csi_correction.R`](R/apply_csi_correction.R)). The corrected rasters
> are in `solarAnnualCSI` and differ from the originals by a per-tile scalar
> spanning **0.62 to 1.05**. That is not a uniform shift, so the old pyramid
> could not be adjusted after the fact and the map was re-tiled from scratch.
>
> **The base zoom changed on 2026-08-29**, from 14 to 15. The 100-1500 ramp
> (previous change) made full-zoom banding more noticeable, and the fix is a
> finer base zoom rather than a different resampling method - see section 2.2.
> This is a same-source, same-colours rebuild: only the tile pyramid's depth
> changes, not what any pixel means.
>
> Everything below describes the current (z15) build; the superseded figures
> from both earlier builds are kept in section 5 for comparison rather than
> deleted.

---

## 1. What the source actually is

| property | value |
|---|---|
| tiles | 2,251 GeoTIFFs, OS National Grid 10 km squares |
| total size | 170 GB |
| mosaic dimensions | 263,000 x 568,932 px (**150 billion pixels**) |
| resolution | 2 m |
| CRS | EPSG:27700 |
| band | 1 x Float32, nodata `-9999` |
| units | kWh/m2/year, annual insolation |

2,251 and not 2,253: NS01 and NS05 are entirely sea, have no land pixel to
correct, and produce no output.

Measured value distribution over the whole mosaic (`gdalinfo -approx_stats`),
before and after the clear-sky-index correction:

| statistic | old (`solarAnnual`) | **new (`solarAnnualCSI`)** |
|---|---:|---:|
| minimum | 98.2 | 87.7 |
| maximum | 2288.0 | **1743.2** |
| mean | 1074.2 | **965.2** |
| std. dev. | 271.9 | 238.9 |

**The maximum fell from 2288 to 1743.** That is what forced the colour domain to
be redefined in section 3: on the corrected data a 0-2000 ramp would leave its
top 13% unreachable, so the reddest colour on the map would be one no pixel ever
takes.

---

## 2. The pipeline, and why it differs from GBDEM

GBDEM's README is the starting point, as instructed. Its shape is kept —
mosaic, reproject, encode to RGB, XYZ tiles, WebP, MBTiles, PMTiles — but three
stages are done differently. GBDEM warns its own instructions are a rough guide
recovered by trial and error, and two of its steps do not survive contact with
this dataset.

### 2.1 Everything before tiling stays virtual

GBDEM materialises a ~95 GB mosaic, then a reprojected copy, then an RGB copy,
and notes the warp alone "will still take a day". At Float32 this mosaic is
~600 GB uncompressed, and the RGBA render of it larger still.

None of it needs writing. A GDAL **VRT** is an XML description evaluated on
read, so the mosaic, the reprojection and the colour mapping chain together and
the pixels are computed exactly once, by the tiler that consumes them.

```
solar_27700.vrt       gdalbuildvrt over 2,253 tiles
  -> solar_3857.vrt        gdalwarp -of VRT   (reprojection, value space)
  -> solar_rgba_3857.vrt   gdaldem color-relief -of VRT -alpha
```

**Total on disk for all three: 1.1 MB.** The RGBA product is 206,492 x 421,532
in EPSG:3857 with R/G/B/Alpha bands. `gdaldem color-relief` supporting
`-of VRT` is the load-bearing detail — it is what removes the last big
intermediate.

### 2.2 Reproject before colouring, not after

Order is not free. Resampling has to happen on **insolation values**, then the
colour table applied to the result. Colouring first would make the warp average
RGB triples, and Turbo is not linear in colour space, so the average of the
colours for 800 and 1200 kWh/m2 is not the colour for 1000. Doing it in the wrong
order produces plausible-looking but wrong colours everywhere the surface is
rough — which, at 2 m over buildings, is everywhere.

Warp settings:

- `-tr` is the exact resolution of the tiler's base zoom for 512 px tiles
  (`156543.03392804097 / 2^maxZoom / 2`). Matching the tile grid means the
  tiler does no further rescaling when cutting base tiles. Base zoom was 14
  (4.77731426716 m/px) on the first two builds and is **15** (2.38865713391
  m/px) from 2026-08-29, because the 100-1500 ramp made banding at full zoom
  more visible than the 0-2000 one did, and z14 was already matched to the 2 m
  source about as tightly as it can be - the fix is a finer base zoom, not a
  different resampling method. `$maxZoom` is a parameter in
  `scripts/01_build_vrts.ps1` and must be kept in sync with the same parameter
  in `scripts/02_tiles.ps1`.
- `-r average` — this is a downsample even at z15 (2 m native to ~2.39
  projected m/px at GB latitudes). `near`, which GBDEM uses, aliases badly on a
  surface this noisy.
- `-srcnodata -9999 -dstnodata -9999` — carried through so the colour table can
  make it transparent.

### 2.3 rio-rgbify is the wrong tool, as suspected

`rio-rgbify` encodes a number into RGB in pseudo base-256 so the client can
decode it back. That is right for a DEM, where MapLibre reads elevation out of
the pixel. Here the pixel's **colour is the message**, so the encoding step is a
colour table applied with `gdaldem color-relief`. That also means the mapping is
a plain text file rather than an arithmetic convention — see below.

### 2.4 gdal2tiles instead of the QGIS toolbox

GBDEM drives tiling from `qgis:tilesxyzdirectory`, because "the mbtiles version
can't create high-resolution (512x512) images". `gdal2tiles --tilesize=512`
removes that constraint, needs no QGIS session, and reads the VRT chain
directly.

One trap: `--xyz` is required. Without it gdal2tiles emits TMS numbering and the
layer appears vertically flipped in MapLibre.

---

## 3. The colour mapping (this is the definition)

Ramp: **Turbo**, from `viridisLite::turbo()`, which is Google's published
colormap - the same one QGIS ships. Generated by
[`RScripts/make_turbo_ramp.R`](RScripts/make_turbo_ramp.R) into
`F:\DTM_DSM\large_rasters\SolarCSI	urbo_100_1500.txt` (the `gdaldem` colour
table) and `turbo_100_1500.csv` (the same data for reference). The domain is in
the filename deliberately, so any tile pyramid can be traced back to the ramp
that coloured it.

**Domain: 100 to 1500 kWh/m2/year, linear**, over 65 evenly spaced anchors with
linear interpolation between them.

| kWh/m2/year | colour | RGB |
|---:|---|---|
| 100 | `#30123B` | 48, 18, 59 |
| 187.5 | `#4040A2` | 64, 64, 162 |
| 275 | `#466BE3` | 70, 107, 227 |
| 362.5 | `#4293FF` | 66, 147, 255 |
| 450 | `#28BBEC` | 40, 187, 236 |
| 537.5 | `#18DCC3` | 24, 220, 195 |
| 625 | `#31F299` | 49, 242, 153 |
| 712.5 | `#6BFE64` | 107, 254, 100 |
| 800 | `#A2FC3C` | 162, 252, 60 |
| 887.5 | `#CCED34` | 204, 237, 52 |
| 975 | `#EDD03A` | 237, 208, 58 |
| 1062.5 | `#FDAD35` | 253, 173, 53 |
| 1150 | `#FB8022` | 251, 128, 34 |
| 1237.5 | `#EC520F` | 236, 82, 15 |
| 1325 | `#D23105` | 210, 49, 5 |
| 1412.5 | `#AC1701` | 172, 23, 1 |
| 1500 | `#7A0403` | 122, 4, 3 |

`nodata (-9999)` maps to `0 0 0 0`, fully transparent, so sea and unmapped
ground drop out instead of rendering as the bottom of the ramp.

### Why the domain changed from 0-2000

The first build used 0-2000, from a brief asking that dark red be "the maximum
value in any of the rasters ~2000". That never quite fitted - the measured
maximum was 2288, so dark red meant "2000 or above" rather than "the maximum" -
and after the clear-sky-index correction it fitted much worse: the corrected
maximum is 1743, which would leave the top 13% of the ramp unreachable, making
the reddest colour on the map one that no pixel ever takes.

Re-chosen against the measured percentiles of the corrected mosaic instead:

| domain | central 90% of pixels | clipped low | clipped high |
|---|---:|---:|---:|
| 0-2000 (old) | 41% of ramp | none | none |
| 0-1750 | 47% of ramp | none | none |
| **100-1500 (chosen)** | **58% of ramp** | **0.009%** | **0.1%** |

100-1500 buys a lot of contrast for very little clipping. The cost is that
**both ends are now clamps rather than extremes**: `#7A0403` means "1500 or
above" and `#30123B` means "100 or below". A legend must therefore be labelled
with the domain and not with the data range, and the clamped tails must not be
read as flat - roughly one pixel in a thousand is genuinely brighter than the
reddest colour on the map.

This is a one-line change: `vmin`/`vmax` at the top of `make_turbo_ramp.R`.
Rename the outputs to match if it is revisited.

### MapLibre usage

This is a plain colour raster, **not** `raster-dem`. It must be added as
`type: "raster"`; declaring it as `raster-dem` would make MapLibre try to decode
the Turbo colours as terrain heights.

```json
"sources": {
  "solar": {
    "type": "raster",
    "url": "pmtiles://https://www.carbon.place/GBsolar.pmtiles",
    "tileSize": 512,
    "minzoom": 5,
    "maxzoom": 15
  }
}
```

---

## 4. Output stages

| stage | tool | output |
|---|---|---|
| 1. mosaic + reproject + colour | `gdalbuildvrt`, `gdalwarp`, `gdaldem` | 1.1 MB of VRT |
| 2. XYZ tiles | `gdal2tiles --xyz --tilesize=512 -z 5-15` | PNG pyramid |
| 3. WebP | ImageMagick `mogrify` in WSL | WebP pyramid |
| 4. MBTiles | `mb-util --image_format=webp --scheme=xyz` | `.mbtiles` |
| 5. PMTiles | `pmtiles-convert` | `.pmtiles` |

### Lossy WebP, deliberately

GBDEM uses `webp:lossless=true`, and must: its tiles are Terrarium-encoded
elevation, so one altered byte is a wrong height. These tiles are a colour ramp,
only ever looked at, never decoded — so lossy WebP is safe here. On imagery this
noisy that is roughly a 10x difference in the size of the deliverable.
`scripts/03_webp_mbtiles_pmtiles.sh lossless` reproduces GBDEM's setting.

### Tooling notes for whoever runs this next

- GDAL is **not on PATH**; it lives in `C:\OSGeo4W\bin` (GDAL 3.4.0).
- `gdal2tiles.bat` in `C:\OSGeo4W\bin` does not resolve. Invoke the module after
  sourcing the environment:
  `cmd /c "call C:\OSGeo4W\bin\o4w_env.bat && python -m osgeo_utils.gdal2tiles ..."`
- GDAL 3.4 predates the WebP tile driver (added 3.6), which is why tiles are
  written as PNG and converted afterwards — the same reason GBDEM gives.
- WSL `Ubuntu-22.04` already had `mogrify`, `mb-util` and `tippecanoe`.
  `tippecanoe` is vector-only and irrelevant here. **`pmtiles` was missing**; I
  installed it user-level with `pip3 install --user pmtiles` (v3.7.0), giving
  `~/.local/bin/pmtiles-convert`.
- **`Ubuntu-22.04` is not the default WSL distro** - `Ubuntu-18.04` is, and it
  has `mb-util` and `pmtiles-convert` but neither `mogrify` nor `sqlite3`. A
  plain `wsl -- bash ...` therefore runs a long way before failing on the WebP
  conversion. Always pass `-d Ubuntu-22.04`, and confirm with
  `wsl --list --verbose`.
- Passing shell one-liners from PowerShell into `wsl -- bash -lc "..."` mangles
  quoting. Write a `.sh` file, `sed -i s/\r$//` it, then `wsl -- bash file`.

---

## 5. Results

### The 2026-08-29 rebuild, base zoom 15 (current)

Same source (`solarAnnualCSI`), same 100-1500 ramp, same `SolarCSI` work
directory - the only change from the previous build is `$maxZoom` in
`scripts/01_build_vrts.ps1` and `scripts/02_tiles.ps1`, 14 -> 15, made because
the higher-contrast ramp showed visible blur/banding at full zoom, and z14 was
already matched to the 2 m source about as tightly as a base zoom can be (see
section 2.2). Everything upstream of the VRT chain - the raster values, the
correction, the colour table - is untouched.

The warped mosaic is **412,983 x 843,063** in EPSG:3857 (double the z14
build's linear dimensions, off by one pixel each way from GDAL's extent
rounding), so the z15 base tile grid is
**807 x 1,647 = 1,329,129** tiles, essentially 4x the z14 base's 332,896 (ratio
3.99, the shortfall from exactly 4x being rounding at the tile-grid edges) -
gdal2tiles' cost is dominated by the base zoom, which is resampled directly
from the source VRT, so this build cost roughly 4x the z5-14 one rather than a
proportionate "one more level" amount.

| stage | time | output |
|---|---:|---|
| VRT chain (mosaic + warp + colour) | seconds | 1.12 MB |
| `gdal2tiles` z5-15 | _running_ | _pending_ |
| PNG -> WebP q85 (14 workers) + `disk_to_mbtiles` | _pending_ | _pending_ |
| `mbtiles_to_pmtiles` | _pending_ | _pending_ |

### The 2026-08-25/26 rebuild, base zoom 14 (superseded 2026-08-29)

Built from `solarAnnualCSI` with the 100-1500 ramp, into
`F:\DTM_DSM\large_rasters\SolarCSI\`. The old `Solar` directory is left
untouched until this build is verified and published, so the superseded
`GBsolar.pmtiles` remains available for rollback.

The warped mosaic is **206,492 x 421,532** in EPSG:3857, identical to the first
build, so the tile pyramid geometry and counts below are unchanged - only the
pixel values and their colours differ.

Completed 2026-08-26. Final deliverable:
**`F:\DTM_DSM\large_rasters\SolarCSI\GBsolar.pmtiles`, 6.98 GiB.**

| stage | time | output |
|---|---:|---|
| VRT chain (mosaic + warp + colour) | seconds | 1.12 MB |
| `gdal2tiles` z5-14 | 10h 08m | 444,879 PNG, 51.43 GB |
| PNG -> WebP q85 (14 workers) + `disk_to_mbtiles` | 4h 55m | 7.3 GB WebP, 7.28 GiB MBTiles |
| `mbtiles_to_pmtiles` | 27m | **6.98 GiB** |

Total 08:30 on 2026-08-25 to 00:06 on 2026-08-26, 15h 36m.

Tile counts per zoom are identical to the first build - same mosaic geometry, so
the same pyramid - and PMTiles deduplication again collapses 444,879 addressed
tiles to **122,669 entries over 117,400 distinct blobs**, because most of the
bounding box is sea and every such tile is byte-identical.

Verified by reading the archive back with
[`scripts/04_verify_pmtiles.py`](scripts/04_verify_pmtiles.py): `TileType.WEBP`,
zoom 5-14, bounds `-6.8637,49.9291,1.9979,60.2151`, all nine metadata keys
present, and valid WebP returned for London, Bristol, Edinburgh, Cardiff,
Inverness and Norwich at both z10 and z14. Sample tiles rendered from the
archive are in `SolarCSI/samples/`: the z14 London tile resolves the Thames,
individual roof pitches and the ribs of a station train shed, and reads
distinctly warmer than the z14 Edinburgh tile - which is the between-tile
variation the old normalisation was erasing, now visible on the map.

### The contrast cost about 1.2 GB

The deliverable grew from 5.81 to 6.98 GiB, and the PNG pyramid from 49.12 to
51.43 GB, even though the geometry is unchanged. That is the 100-1500 domain
being paid for: spreading the data over more of Turbo makes neighbouring pixels
differ more in colour, so the imagery carries more high-frequency detail and
lossy WebP compresses it less well (7.0x rather than 8.2x). Worth knowing before
tightening the domain any further - contrast and archive size trade against each
other directly here.

### The 2026-08-09 build (superseded)

Kept for comparison. This one was coloured from the **uncorrected** `solarAnnual`
rasters on a 0-2000 domain, so its levels are wrong by the per-tile scalar
described at the top of this file. Final deliverable was
**`F:\DTM_DSM\large_rasters\Solar\GBsolar.pmtiles`, 5.81 GB.**

| stage | time | output |
|---|---:|---|
| VRT chain (mosaic + warp + colour) | seconds | 1.1 MB |
| `gdal2tiles` z5-14 | 8h 23m | 444,879 PNG, 49.12 GB |
| PNG -> WebP q85, 14 workers | 32m | 5.99 GB (**8.2x smaller**) |
| `disk_to_mbtiles` | ~35m | 6.11 GB |
| `mbtiles_to_pmtiles` | 69s | **5.81 GB** |

Tile pyramid, which matches the expected geometry exactly
(z14 = ceil(206492/512) x ceil(421532/512) = 404 x 824 = 332,896):

| zoom | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| tiles | 4 | 12 | 28 | 91 | 338 | 1,352 | 5,408 | 21,114 | 83,636 | 332,896 |

Verified by reading the finished archive back: `TileType.WEBP`, zoom 5-14,
bounds `-6.8626,49.9346,1.5680,60.1264`, all ten metadata keys present, and
valid WebP returned for London, Bristol, Edinburgh and Cardiff at z10 and z14.
A z14 London tile extracted from the archive renders correctly — Thames orange,
north-facing roof pitches blue, station train-shed roofs resolved individually.

**PMTiles deduplication does the heavy lifting on empty sea.** 444,879 addressed
tiles collapse to 122,669 entries over 117,400 distinct blobs, because most of
the bounding box is sea or unmapped and every such tile is byte-identical. This
is why it was not worth excluding empty tiles at the gdal2tiles stage.

### How much the lossy WebP actually costs

Measured on a dense 210 KB London z14 tile, RMSE against the original PNG:

| quality | size | RMSE (0-1) | mean channel error |
|---|--:|--:|--:|
| 80 | 13 KB | 0.0192 | 4.9 / 255 |
| **85** | **17 KB** | **0.0178** | **4.5 / 255** |
| 90 | 23 KB | 0.0170 | 4.3 / 255 |
| 95 | 34 KB | 0.0165 | 4.2 / 255 |
| lossless | 140 KB | 0 | 0 |

Error barely falls from q85 to q95 while size doubles — the source is too
high-frequency for lossy WebP to ever match it — so q85 is the sensible point.
4.5/255 of colour distance is on the order of 35 kWh/m2/year read off the ramp,
against a data standard deviation of 272. Fine for a visualisation; if the layer
were ever to be read quantitatively, rebuild with `lossless`.

### One failure worth recording

`mb-util` only writes an MBTiles metadata table if the tile directory contains a
`metadata.json`, and `gdal2tiles` does not produce one. The MBTiles therefore
came out with all 444,879 tiles correct and an **empty metadata table**, and
`pmtiles-convert` died on `KeyError: 'format'`. Nothing was wrong with the
tiles. The fix is to populate `metadata` (`format`, `bounds`, `center`,
`minzoom`, `maxzoom`, `name`, `type`, `version`) directly in the SQLite file and
re-run the conversion — 69 seconds, no need to rebuild anything.

---

## 6. Validation against PVGIS

The clear-sky-index correction makes a tile mean physically meaningful for the
first time - under the old normalisation each tile was scaled to its own mean by
construction, so no external comparison could fail. It is now worth testing.

For tile *g*, `mean(corrected) = T_g x E_ERA5_smooth`, so dividing the tile mean
by its terrain factor recovers the modelled **flat-ground** annual GHI. PVGIS
answers the same question from **SARAH2**, a Meteosat-derived radiation
database independent of the ERA5 field used here. Script:
[`RScripts/validate_pvgis.R`](RScripts/validate_pvgis.R); results in
[`data/validation_pvgis.csv`](data/validation_pvgis.csv).

24 tiles, stratified into 12 latitude bands from Cornwall to Shetland, year 2020:

| statistic | value |
|---|---:|
| sites compared | 22 of 24 |
| correlation | **0.956** |
| mean bias | **+1.28%** |
| median bias | +0.95% |
| sd of differences | 4.20% |
| range | -4.85% to +8.20% |

A mean bias of about **1%** on an absolute radiation level is a good result, and
it is a genuinely independent check: nothing in this pipeline was tuned to
SARAH2.

**But the residuals are not random.** They trend with latitude at
**+1.10% per degree** (p = 0.0001, R2 = 0.56):

| | n | mean difference |
|---|---:|---:|
| south of 55N | 11 | **-1.05%** |
| north of 55N | 11 | **+3.62%** |

So relative to SARAH2 this map reads slightly low in southern England and
several percent high in northern Scotland. Residual scatter about the latitude
fit is 2.79%. The likely cause is a known disagreement between ERA5 and
satellite radiation retrievals at high latitude, where low sun angles, snow and
persistent cloud all make the retrieval harder; the ERA5 field is the input
here, so the pattern is inherited rather than introduced. Part of the scatter is
also methodological - PVGIS answers for a *point* while these are 10 km *means*,
and the two cannot agree exactly over broken terrain.

**Two sampling traps, both real.** PVGIS refuses any query it considers offshore
("Location over the sea"), and a first attempt using tile *centroids* lost five
of 24 sites: sea is nodata in these rasters, so the tile mean already describes
land only, but the geometric centre of a coastal or island square often sits in
the water. The query point must be the centroid of the tile's **land** pixels,
falling back to actual sampled land pixels for horseshoe-shaped coastlines. Two
sites (SW65, 11.6% land; HY63, Orkney) are still refused and are reported as
missing rather than quietly dropped.

## 7. Corrections to the GBsolar README

- It points at `../GBDEM/README.md`. There is no GBDEM beside this repo; it is
  at **`F:\GitHub\mem48\GBDEM`**.
- It cites `/RScripts/built_isolation_tiles_10k.R`. The file is really
  `RScripts/built_isolt` + `ation_tiles_10k.R` — i.e.
  **`built_isoltation_tiles_10k.R`**, with "isolation" misspelt "isoltation".
  Worth renaming; a search for the name in the README finds nothing.
