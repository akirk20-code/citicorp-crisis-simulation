# How to Extract 3D Buildings from NYC 3D Building Model

## Quick Start

```bash
# Extract buildings within 300m of a lat/lon
python extract_buildings.py --lat 40.7580 --lon -73.9690 --radius 300

# Named output, STL only
python extract_buildings.py --lat 40.7484 --lon -73.9857 --radius 200 --name empire_state --format stl
```

Output: `buildings.stl` (binary STL) + `buildings.gml` (CityGML fragments)

## What You Get

Each building has LOD2 geometry:
- **WallSurface** — vertical faces with setbacks
- **RoofSurface** — actual roof shapes (gabled, hipped, flat, etc.)
- **GroundSurface** — footprint at ground elevation
- **Attributes** — BIN, year built, floors, landmark status, assessed value

## Data Source

**NYC 3D Building Model** (v20v5)
- Photogrammetric survey from August 2014 aerial imagery
- ~1 million buildings across NYC's 5 boroughs
- LOD1 for most buildings, LOD2 for ~100 iconic buildings
- Coordinate system: EPSG:2263 (NY State Plane Long Island, US Survey Feet)
- Source: [georocket/new-york-city-model-enhanced](https://github.com/georocket/new-york-city-model-enhanced)

## How It Works

1. **Lat/lon to State Plane** — Converts your WGS84 coordinates to EPSG:2263 feet
2. **DA district detection** — NYC is divided into 20 "DA" districts; the script auto-detects which file(s) to search
3. **Streaming extraction** — Scans the large GML file (~0.5-1.2 GB each) line by line, extracting buildings whose coordinates fall within your radius
4. **STL conversion** — Triangulates all polygon surfaces and writes binary STL in local meters (centered at your query point, Z=0 at ground)

## Requirements

- Python 3.8+ (stdlib only, no pip installs)
- DA district GML files in `_citygml/` directory (auto-downloaded if missing)

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--lat` | required | Center latitude (WGS84) |
| `--lon` | required | Center longitude (WGS84) |
| `--radius` | 300 | Search radius in meters |
| `--name` | buildings | Output filename prefix |
| `--format` | both | `gml`, `stl`, or `both` |
| `--data-dir` | `_citygml/` | Where DA district files live |
| `--output-dir` | current dir | Where to write output |

## Example: Citicorp Center (100m radius)

```
python extract_buildings.py --lat 40.7580 --lon -73.9690 --radius 100

Found 31 buildings
  BIN 1036474: 284m, 60 floors, built 1978 ** INDIVIDUAL LANDMARK
  BIN 1036467: 202m, 47 floors, built 1985
  BIN 1038577: 191m, 46 floors, built 1970
  ...
  Height stats: min=8m, max=284m, mean=66m, median=34m
  Wrote buildings.stl: 4392 triangles
```

## Tips

- **Finding lat/lon**: Right-click any location in Google Maps and copy coordinates
- **Radius**: 100m = 1 block, 300m = ~3 blocks, 500m = ~5 blocks
- **Performance**: Each DA file takes 10-20 seconds to scan. The script auto-detects which district(s) to search
- **Large extractions**: For radius > 500m, expect 100+ buildings and 20k+ triangles
- **OpenFOAM use**: The binary STL output is directly usable in snappyHexMesh
