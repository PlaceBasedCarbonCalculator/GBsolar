#!/bin/bash
# Stage 3: PNG -> WebP -> MBTiles -> PMTiles. Run inside WSL (Ubuntu-22.04).
#
#   bash 03_webp_mbtiles_pmtiles.sh [lossy|lossless]
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
# mb-util wants the directory to contain only the tile pyramid.
rm -f "$OUT.mbtiles"
mb-util --image_format=webp --scheme=xyz "$TILES" "$OUT.mbtiles"

echo "=== PMTiles ==="
rm -f "$OUT.pmtiles"
pmtiles-convert "$OUT.mbtiles" "$OUT.pmtiles"

ls -lh "$OUT.mbtiles" "$OUT.pmtiles"
