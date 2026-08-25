#!/usr/bin/env python3
"""Stage 4: read the finished PMTiles archive back and check it.

    wsl -- python3 scripts/04_verify_pmtiles.py [path]

The first build was verified by hand. This does the same checks reproducibly,
because every one of them has a plausible failure behind it:

  tile type     lossy WebP is the whole size argument; a PNG archive means the
                conversion silently did not happen.
  zoom range    gdal2tiles -z is easy to mistype and the loss is invisible
                until someone zooms in.
  metadata      an empty metadata table is what killed pmtiles-convert on the
                first build, so it is checked rather than assumed.
  real tiles    an archive can be structurally perfect and still hold nothing
                over land, so known city coordinates are fetched and their
                bytes checked for a WebP signature.
"""
import json
import math
import sys

from pmtiles.reader import MmapSource, Reader

PATH = sys.argv[1] if len(sys.argv) > 1 else \
    "/mnt/f/DTM_DSM/large_rasters/SolarCSI/GBsolar.pmtiles"

CITIES = {
    "London":    (51.5074, -0.1278),
    "Bristol":   (51.4545, -2.5879),
    "Edinburgh": (55.9533, -3.1883),
    "Cardiff":   (51.4816, -3.1791),
    "Inverness": (57.4778, -4.2247),
    "Norwich":   (52.6309,  1.2974),
}


def lonlat_to_tile(lat, lon, z):
    n = 2 ** z
    x = int((lon + 180.0) / 360.0 * n)
    r = math.radians(lat)
    y = int((1.0 - math.asinh(math.tan(r)) / math.pi) / 2.0 * n)
    return x, y


def main():
    failures = []
    with open(PATH, "rb") as f:
        reader = Reader(MmapSource(f))
        hdr = reader.header()
        meta = reader.metadata()

        print(f"archive : {PATH}")
        tt = hdr["tile_type"]
        print(f"tile type : {tt}")
        if "WEBP" not in str(tt).upper():
            failures.append(f"tile type is {tt}, expected WEBP")

        zmin, zmax = hdr["min_zoom"], hdr["max_zoom"]
        print(f"zooms   : {zmin}-{zmax}")
        if (zmin, zmax) != (5, 14):
            failures.append(f"zoom range {zmin}-{zmax}, expected 5-14")

        bounds = (hdr["min_lon_e7"] / 1e7, hdr["min_lat_e7"] / 1e7,
                  hdr["max_lon_e7"] / 1e7, hdr["max_lat_e7"] / 1e7)
        print("bounds  : %.4f,%.4f,%.4f,%.4f" % bounds)
        if not (-9 < bounds[0] < -5 and 48 < bounds[1] < 51
                and 0 < bounds[2] < 3 and 59 < bounds[3] < 62):
            failures.append(f"bounds {bounds} do not look like Great Britain")

        print(f"entries : {hdr['addressed_tiles_count']:,} addressed, "
              f"{hdr['tile_entries_count']:,} entries, "
              f"{hdr['tile_contents_count']:,} distinct blobs")

        want = {"format", "bounds", "center", "minzoom", "maxzoom",
                "name", "type", "version"}
        have = set(meta) if isinstance(meta, dict) else set()
        print(f"metadata: {len(have)} keys -> {json.dumps(meta)[:160]}")
        missing = want - have
        if missing:
            failures.append(f"metadata missing keys: {sorted(missing)}")

        print("\nsampling real tiles:")
        for z in (10, 14):
            for name, (lat, lon) in CITIES.items():
                x, y = lonlat_to_tile(lat, lon, z)
                data = reader.get(z, x, y)
                if not data:
                    failures.append(f"{name} z{z} ({x},{y}) is empty")
                    print(f"  z{z:<2} {name:<10} ({x},{y})  EMPTY")
                    continue
                sig = data[:4] == b"RIFF" and data[8:12] == b"WEBP"
                if not sig:
                    failures.append(f"{name} z{z} is not WebP")
                print(f"  z{z:<2} {name:<10} ({x},{y})  "
                      f"{len(data):>7,} bytes  {'WebP' if sig else 'NOT WEBP'}")

    print()
    if failures:
        print(f"FAILED ({len(failures)}):")
        for f_ in failures:
            print(f"  - {f_}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
