#!/bin/bash
# Stage 3: PNG -> WebP -> MBTiles -> PMTiles. Run inside WSL (Ubuntu-22.04).
#
#   wsl -d Ubuntu-22.04 -- bash #     /mnt/f/GitHub/.../scripts/03_webp_mbtiles_pmtiles.sh [lossy|lossless]
#
# `-d Ubuntu-22.04` is NOT optional. The default distro on this machine is
# Ubuntu-18.04, which has mb-util and pmtiles-convert but no mogrify and no
# sqlite3, so a plain `wsl -- bash ...` gets most of the way through and then
# fails on the WebP conversion. Check with `wsl --list --verbose`.
#
# Do not try to pass this as a one-liner. Quoting from PowerShell or Git Bash
# into `wsl -- bash -lc "..."` mangles the command and Git Bash also rewrites
# /mnt/... into a Windows path; run it as a file, as above.
#
# Lossy is the default here, and that is a deliberate departure from GBDEM.
# GBDEM's tiles are Terrarium-encoded elevation: the RGB triple IS the number
# MapLibre decodes, so a single altered byte becomes a wrong elevation and
# lossless is mandatory. These tiles are a colour ramp - the pixel is only ever
# looked at, never decoded - so lossy WebP is safe, and on imagery this noisy it
# is the difference between a ~30 GB and a ~3 GB deliverable. Pass "lossless" to
# reproduce GBDEM's setting.
#
# WORK is SolarCSI, the clear-sky-index-corrected rebuild. The old `Solar`
# directory holds the superseded pyramid and is left alone.

set -euo pipefail
MODE="${1:-lossy}"
WORK=/mnt/f/DTM_DSM/large_rasters/SolarCSI
TILES="$WORK/tiles_png"
export PATH="$PATH:$HOME/.local/bin"

if [ "$MODE" = "lossless" ]; then
  DEFINE="webp:lossless=true"
  OUT="$WORK/GBsolar_lossless"
else
  DEFINE="webp:lossless=false"
  OUT="$WORK/GBsolar"
fi

echo "=== converting PNG -> WebP ($MODE) ==="
# -quality 85 is ignored when lossless=true. Parallelised with xargs because
# a single mogrify walk over ~300k tiles is hours of one core.
find "$TILES" -type f -name '*.png' -print0 |
  xargs -0 -P 14 -n 200 mogrify -format webp -define "$DEFINE" -quality 85
find "$TILES" -type f -name '*.png' -delete

echo "=== MBTiles ==="
# mb-util only writes an MBTiles metadata table if the tile directory contains a
# metadata.json, and gdal2tiles does not produce one. Without it the archive
# comes out with every tile correct and an EMPTY metadata table, and
# pmtiles-convert then dies on `KeyError: 'format'` - which is exactly what
# happened on the first build. Writing it here costs nothing and removes the
# failure. Bounds are the WGS84 extent of solar_rgba_3857.vrt, from
# `gdalinfo -json`; regenerate them if the source tile set ever changes.
# maxzoom must match $maxZoom in scripts/01_build_vrts.ps1 and 02_tiles.ps1.
cat > "$TILES/metadata.json" <<'JSON'
{
  "name": "GBsolar",
  "description": "GB annual solar insolation, kWh/m2/year, Turbo ramp 100-1500",
  "type": "overlay",
  "version": "1",
  "format": "webp",
  "bounds": "-6.8637,49.9291,1.9979,60.2151",
  "center": "-2.4329,55.0721,10",
  "minzoom": "5",
  "maxzoom": "15"
}
JSON

rm -f "$OUT.mbtiles"
mb-util --image_format=webp --scheme=xyz "$TILES" "$OUT.mbtiles"

echo "=== PMTiles ==="
rm -f "$OUT.pmtiles"
pmtiles-convert "$OUT.mbtiles" "$OUT.pmtiles"

ls -lh "$OUT.mbtiles" "$OUT.pmtiles"

# The metadata table is the thing that broke last time, so check it rather than
# assume it.
echo "=== metadata written ==="
sqlite3 "$OUT.mbtiles" "SELECT name || ' = ' || value FROM metadata ORDER BY name;" 2>/dev/null ||
  echo "(sqlite3 not available - check metadata manually)"
