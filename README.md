# GBsolar

A 2 m resolution map of the solar energy reaching every square metre of ground and
roof in Great Britain, built from the surface model produced by
[GBDEM](https://github.com/PlaceBasedCarbonCalculator/GBDEM).

The finished layer is the **Solar Potential** layer in the
[Retrofit Explorer](https://www.carbon.place/retrofit/) on
[Carbon & Place](https://www.carbon.place), and it is described for a general
audience in the [manual](https://www.carbon.place/manual/#solarpotential).

Values are **kilowatt hours per square metre per year (kWh/m²/year)**. For
scale, a flat unshaded site in Great Britain receives roughly 750 to 1,100
kWh/m² a year; the mean across the whole national mosaic is 965 and the
maximum, on steep south-facing surfaces, is 1,743 (measured on the
clear-sky-index-corrected mosaic — see *Cloud correction from ERA5* below).

Because the model works from a *surface* model rather than a terrain model, it
resolves individual roof pitches: at full zoom the south-facing slope of a roof
reads orange or red while the north-facing slope of the same roof reads blue,
and the shadow a tall building casts over its neighbours appears as a cool
patch.

---

## What this is, and what it is not

This layer maps **the resource, not the yield**. It models how much solar energy
arrives at a surface over a year, given that surface's slope, its aspect, the
shading cast on it by terrain and by nearby buildings and trees, and the local
cloudiness.

It says nothing about what a photovoltaic installation would actually generate.
Converting irradiation to electricity needs assumptions about panel efficiency,
system losses, inverter behaviour, how much of a roof is usable, and whether the
roof can carry the load. None of those are applied here. A bright roof on this
map is a candidate for investigation, not an estimate of output.

It is also a **visualisation**, not a data product. The published tiles encode
values as colours (see *Colour mapping* below), and are stored as lossy WebP, so
they should be read to the nearest colour band. If you need the numbers, use the
source rasters rather than sampling the tiles.

---

## Method

### 1. Clear-sky irradiation from the surface model

`R/insolation_calcs_intergrated.R` is the function that runs in production
(`insolation_annual_strategy()`), driven tile by tile from
`RScripts/built_isoltation_tiles_10k.R`.

For each of the 2,253 OS National Grid 10 km squares:

1. Load the 2 m digital surface model tile from GBDEM and derive slope and
   aspect with `terra::terrain()`.
2. Import the DSM, slope, and aspect into GRASS GIS and run
   [`r.sun`](https://grass.osgeo.org/grass-stable/manuals/r.sun.html) in daily
   integration mode, hourly step, for one representative day per month (the 15th
   by default). `r.sun` handles the solar geometry and, importantly, casts
   shadows from the surface itself, which is where the roof-pitch and
   building-shadow detail comes from.
3. Weight each representative day by the number of days in its month and sum, to
   give a clear-sky annual total for every pixel.

Twelve representative days rather than 365 is the compromise that makes a
national run feasible. It captures the seasonal cycle of sun angle and day
length well; it will not capture day-to-day variation, which is handled at the
next step and only in the aggregate.

### 2. Cloud correction from ERA5

Clear-sky totals are far too high for Britain. Each tile is therefore scaled to
match the surface solar radiation downwards (`SSRD`) recorded for that location
in the ECMWF [ERA5](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)
reanalysis, so that the tile's *mean* matches the reanalysis while the *pattern*
within the tile comes from the surface model.

That description is of the *original* method, and it was wrong. Scaling each
tile to its own mean forced every 10 km square to the same average clear-sky
index, erasing the real differences between them - a flat fen and a Highland
glen came out equally sunny on average. `R/apply_csi_correction.R` documents the
defect and the algebra that fixes it: each tile is now normalised by a **clear-sky
index**, the ratio of the reanalysis total to clear-sky radiation on *flat*
ground, so between-tile variation survives.

Because the defect was a per-tile scalar, the expensive 2 m within-tile
pattern could be kept and only the level of each tile rescaled - a 27-day
re-run became a 10-day one, and the resulting corrections span **0.62 to 1.05**.
That is the size of the error that had been there.

**This correction is complete.** All 2,251 land tiles were reprocessed between
2026-08-12 and 2026-08-25 by
[`RScripts/run_csi_correction.R`](RScripts/run_csi_correction.R), with no
errors; the per-tile record is in
[`data/csi_correction_progress.csv`](data/csi_correction_progress.csv). The
corrected rasters are `solarAnnualCSI` and are what the published map is built
from. `solarAnnual` is superseded and should not be used.

**Note the resolution mismatch.** ERA5 is a global reanalysis on a grid of tens
of kilometres. Applying it per 10 km tile can leave visible steps at tile
boundaries, which is why the correction code de-blocks the reanalysis field
(`data/era5_annual_smooth.csv`, built by `R/era5_field.R`) rather than using raw
per-tile values.

### 3. Tiling for the web

Turning 2,251 corrected rasters (170 GB, a 150-billion-pixel mosaic) into a single
PMTiles archive is documented in full in **[METHOD.md](METHOD.md)**, which
covers the VRT chain that avoids materialising any intermediate raster, why
reprojection must happen before colouring, the colour table, the tile pyramid,
the measured cost of lossy WebP, and the failures encountered.

**Base zoom is 15** (previously 14, raised 2026-08-29). The tighter 100-1500
colour domain below made the resampling blur beyond the tileset's top zoom
more visible than the original 0-2000 domain did, and z14 was already matched
to the 2 m source about as tightly as a base zoom can be — the fix is one more
zoom level, not a different resampling method. `$maxZoom` is a parameter of
`scripts/01_build_vrts.ps1` and `scripts/02_tiles.ps1`; the two must agree, or
the tiler either rescales the base tiles again or samples a source coarser
than the zoom level implies. It roughly quadruples the tiling cost, because
gdal2tiles' time is dominated by the base zoom — see METHOD.md for the numbers.

---

## Colour mapping

The published layer encodes value as colour using the
[Turbo](https://research.google/blog/turbo-an-improved-rainbow-colormap-for-visualization/)
ramp, **linear over a fixed domain of 100 to 1500 kWh/m²/year**, generated by
`RScripts/make_turbo_ramp.R`. Nodata is fully transparent, so sea and unmapped
ground drop out rather than rendering as the bottom of the ramp.

The domain is fixed rather than fitted so that the legend means the same thing
every time the layer is rebuilt, and it moved from an original 0-2000 to
100-1500 when the clear-sky-index correction (above) shifted the mosaic's
whole distribution downward — 0-2000 would have left the top 13% of the ramp
unreachable. **Both ends are now clamps, not extremes:**

- **The top of the ramp is a clamp.** The measured maximum is 1,743 kWh/m²/year,
  so the darkest red means "1500 or above", not "the maximum" — about 0.1% of
  pixels are affected.
- **The bottom of the ramp is also a clamp**, unlike the previous domain. The
  measured minimum is 88 kWh/m²/year, so the darkest blue means "100 or below" —
  about 0.009% of pixels are affected.

That trade buys real contrast: the central 90% of pixels now span 58% of the
ramp, against 41% under the old domain. The cost is a legend that must be read
as a domain, not a data range, and roughly 1.2 GB of extra archive size, because
spreading the data across more of Turbo gives lossy WebP more high-frequency
colour detail to compress. Both ends are a one-line change (`vmin`/`vmax` in
`make_turbo_ramp.R`) — rename the output files to match if it is revisited, so
a pyramid can always be traced back to the ramp that coloured it. The full
value-to-colour table and the domain-selection reasoning are in
[METHOD.md](METHOD.md).

### Adding the layer to MapLibre

This is a plain colour raster. It must be declared as `type: "raster"` and
**not** as `raster-dem`, which would make MapLibre try to decode the Turbo
colours as terrain heights.

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

## A note on units and filenames

Everything published from this repository is in **kWh/m²/year**: the rasters, the
colour ramp domain, the figures in [METHOD.md](METHOD.md), and the legend on
Carbon & Place. This was previously documented as Wh/m² throughout, which was
wrong by a factor of a thousand and has been corrected.

Inside the pipeline the intermediate quantities really are in Wh/m²/year — both
the ERA5 annual totals (`SSRD / 3600`) and the clear-sky output of `r.sun`. The
conversion happens at the last step, where `E_clear_annual * scale / 1000`
produces the final raster. Comments in `R/` mark which is which, and the one
remaining inconsistency is the output filename, described below.

## Known issues

**Output filenames say `Whm2` and are wrong.** The rasters are named
`<grid>_annual_insolation_Whm2.tif` but hold kWh/m²/year, a legacy of the units
being mislabelled when the pipeline was written. The name is retained because
`apply_csi_correction.R`, `run_csi_correction.R` and the tiling stage all match
on it, and 2,253 files exist on disk under it; renaming would need all four
changed together. Do not infer the units from the filename. The intermediate
`<grid>_clearsky_annual_Whm2.tif` rasters, by contrast, genuinely are in
Wh/m²/year, as are the ERA5 fields — only the final product was converted.

**Twelve representative days.** The seasonal cycle is sampled monthly. Locations
whose shading changes sharply within a month, and any effect of the specific
choice of the 15th, are not captured.

**Cloud is applied as a per-tile scalar.** Within a 10 km tile, every pixel
receives the same cloudiness correction. Real local variation in cloud, for
example orographic cloud on one side of a hill, is not represented.

**The surface model is a snapshot.** It reflects the buildings and vegetation
present when the LIDAR was flown. New buildings, demolitions, and tree growth
or felling since are absent, and LIDAR coverage itself is uneven: over 95% in
England, about 70% in Wales, and about 40% in Scotland, with gaps filled from
lower-resolution data.

**Lossy tiles.** Published tiles are WebP at quality 85, measured at about
4.5/255 of mean colour error on the original 0-2000/z14 build (see METHOD.md) —
scaled to the current 100-1500 domain that is on the order of 25 kWh/m²/year
read off the ramp, against a data standard deviation of 239; not yet
re-measured on the z15 tileset. Fine for looking at; not suitable for
quantitative sampling. Rebuild with `lossless` if that is needed.

**Some inland squares have no data at all, silently.** Both
`RScripts/built_isoltation_tiles_10k.R` and `RScripts/run_csi_correction.R`
enumerate grid squares from the DSM directory, so a square with no DSM never
enters either pipeline — it does not appear as an error or a flagged gap, it is
simply absent from the published map, indistinguishable at a glance from a
mostly-sea square with real data too thin to see (e.g. NZ90, 0.5% land). A
national scan (2026-08-29) found **169** such squares nationally, **136** with
more than 10% land in GBDEM's DTM (bare earth), so genuinely missing rather
than genuinely at sea. Of those, **107** — including NZ80 in the North York
Moors — already have a `DSM_infill` file sitting unused on disk, because
neither script falls back to it; the remaining **29**, mostly Orkney, the
Western Isles, Caithness and Galloway, have no surface model at all, native or
infilled. Not yet fixed; the DSM pipeline is unchanged pending a decision on
whether to add the `DSM_infill` fallback.

**The level reads high in the far north.** Validated against PVGIS/SARAH2 the
map carries a mean bias of about +1%, but the residuals trend with latitude at
+1.1% per degree: roughly -1% south of 55N and +3.6% north of it. See METHOD.md
section 6.

---

## Repository layout

| Path | Contents |
| --- | --- |
| `R/` | Functions. `insolation_calcs_intergrated.R` holds the production routine; `apply_csi_correction.R` the ERA5 correction and its derivation; `era5_field.R` and `create_df_from_nc.R` the reanalysis handling. |
| `RScripts/` | Scripts that run the functions, in roughly the order they are used. `built_isoltation_tiles_10k.R` is the national loop (note the misspelling of "insolation" in the filename). `run_csi_correction.R` applies the clear-sky-index correction across GB and is resumable. `make_turbo_ramp.R` defines the colour table, and `validate_pvgis.R` checks the result against PVGIS. |
| `sampleData/` | Small inputs, including the per-grid ERA5 extracts. |
| `data/` | `era5_annual_smooth.csv`, the de-blocked reanalysis field; `csi_correction_progress.csv`, the per-tile record of the correction run; `validation_pvgis.csv`, the PVGIS comparison. |
| `scripts/` | The tiling pipeline, run in order: `01_build_vrts.ps1` (mosaic + reproject + colour, as VRTs), `02_tiles.ps1` (`gdal2tiles`), `03_webp_mbtiles_pmtiles.sh` (WSL; PNG → WebP → MBTiles → PMTiles), `04_verify_pmtiles.py` (reads the finished archive back and checks it). |
| `METHOD.md` | The full tiling method, colour table, results, and tooling notes. |

## Requirements

R, with `terra`, `sf`, `rgrass`, `suntools`, `lubridate`, `stringr`,
`viridisLite`, and `jsonlite` (for `validate_pvgis.R`); **GRASS GIS 8.4** for
`r.sun`; and GDAL for the tiling stages.
The national run is long: `gdal2tiles` took 10h08m for the z5-14 pyramid and
roughly four times that for the current z5-15 one, since its cost is dominated
by the base zoom (see METHOD.md); a full 2 m re-run of the insolation
calculation for all 2,253 tiles would take roughly four weeks on the machine
used.

## Related repositories

- [GBDEM](https://github.com/PlaceBasedCarbonCalculator/GBDEM) — the terrain and
  surface models this is built from.
- [build](https://github.com/PlaceBasedCarbonCalculator/build) — the main
  Carbon & Place analysis pipeline.

## Data sources

- Digital surface model derived from LIDAR published by the
  [Environment Agency](https://www.data.gov.uk/dataset/f0db0249-f17b-4036-9e65-309148c97ce4/national-lidar-programme),
  the [Welsh Government](https://datamap.gov.wales/maps/lidar-viewer/), and the
  [Scottish Government](https://remotesensingdata.gov.scot/data#/list), via GBDEM.
- [ECMWF ERA5](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)
  reanalysis, for surface solar radiation and cloud.

## Licence

Code is published under the GNU Affero General Public Licence v3.0; see
[LICENSE](LICENSE). Input data remains under the licences of its publishers.

## Citation

Morgan, M. (2026). Carbon & Place: Data and tools to understand the spatial
variation in carbon footprints. *Environment and Planning B: Urban Analytics and
City Science*, 53(3), 538–554. <https://doi.org/10.1177/23998083251401613>
