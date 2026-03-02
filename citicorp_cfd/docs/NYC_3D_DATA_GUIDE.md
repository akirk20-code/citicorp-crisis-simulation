# NYC 3D Building Data Access Guide

This document describes how to access NYC building data for CFD simulations and which tool to use for different purposes.

---

## Which Tool Should I Use?

| Tool | Best For | Data Source |
|------|----------|-------------|
| **`generators/generate_stl_nyc3d.py`** | Production CFD runs (recommended) | NYC Building Footprints API |
| **`generators/generate_stl_hybrid.py`** | Historical analysis, other locations | CityGML + Socrata API fallback |
| **`generators/generate_stl.py`** | Quick prototyping, offline use | Hardcoded 18 buildings |
| **`tools/extract_buildings.py`** | Research-grade LOD2 geometry | TUM CityGML files |

**For this project, use `generators/generate_stl_nyc3d.py`** — it provides surveyed heights (±2 ft accuracy), correct 45-degree tower rotation, and proper stilt positions.

---

## Data Sources

### 1. NYC Open Data Building Footprints API (Primary)

**Dataset ID:** `5zhs-2jue`
**Format:** GeoJSON via Socrata API
**URL:** https://data.cityofnewyork.us/City-Government/BUILDING/5zhs-2jue

**Key Attributes:**
- `height_roof`: Roof height above ground (feet)
- `groundelev`: Ground elevation (feet)
- `bin`: Building Identification Number
- `cnstrct_yr`: Construction year
- `the_geom`: Footprint geometry (MultiPolygon/Polygon)

**Example Query (Python):**
```python
import requests
params = {
    "$where": "within_box(the_geom, 40.7620, -73.9750, 40.7540, -73.9630)",
    "$limit": 5000,
    "$order": "height_roof DESC"
}
response = requests.get(
    "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson",
    params=params, timeout=60
)
```

**Authentication:** None required (public dataset)
**Accuracy:** ±2 feet for buildings >60 feet (2014 aerial survey)

### 2. NYC 3D Building Model by Community District

**Dataset ID:** `u5j4-zxpn`
**Format:** 3DM (Rhinoceros), CityGML
**URL:** https://data.cityofnewyork.us/City-Government/NYC-3D-Model-by-Community-District/u5j4-zxpn

Full 3D building massing model based on 2014 aerial survey. Hybrid LOD1/LOD2:
- LOD1: Simple block models with flat roofs
- LOD2: Detailed roof geometry for some buildings

### 3. TUM Enhanced NYC 3D Building Model

**URL:** https://github.com/georocket/new-york-city-model-enhanced
**Format:** CityGML with ~90 semantic attributes per building

Used by `tools/extract_buildings.py` for research-grade LOD2 geometry. Large dataset (~20 GB for full city).

---

## Extracting Buildings with extract_buildings.py

### Quick Start

```bash
# Extract buildings within 300m of a lat/lon
python tools/extract_buildings.py --lat 40.7580 --lon -73.9690 --radius 300

# Named output, STL only
python tools/extract_buildings.py --lat 40.7484 --lon -73.9857 --radius 200 --name empire_state --format stl
```

### What You Get

Each building has LOD2 geometry:
- **WallSurface** — vertical faces with setbacks
- **RoofSurface** — actual roof shapes (gabled, hipped, flat)
- **GroundSurface** — footprint at ground elevation
- **Attributes** — BIN, year built, floors, landmark status

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--lat` | required | Center latitude (WGS84) |
| `--lon` | required | Center longitude (WGS84) |
| `--radius` | 300 | Search radius in meters |
| `--name` | buildings | Output filename prefix |
| `--format` | both | `gml`, `stl`, or `both` |
| `--data-dir` | `_citygml/` | Where DA district files live |

### How It Works

1. **Lat/lon to State Plane** — Converts WGS84 to EPSG:2263 (US Survey Feet)
2. **DA district detection** — NYC is divided into 20 "DA" districts; auto-detects which file(s) to search
3. **Streaming extraction** — Scans large GML files (~0.5-1.2 GB each) line by line
4. **STL conversion** — Triangulates surfaces and writes binary STL in local meters

**Requirements:** Python 3.8+ (stdlib only). DA district GML files auto-downloaded if missing.

---

## Data Processing Pipeline

### Coordinate Transformation

```
WGS84 (lat, lon) --> Local Cartesian (x, y, z)
Origin: Citicorp Center (40.7579 N, 73.9690 W)
X-axis: East, Y-axis: North, Z-axis: Up
```

### Building Footprint Processing

1. Parse GeoJSON MultiPolygon/Polygon
2. Transform coordinates: WGS84 to local meters
3. Ensure CCW winding (correct normals)
4. Validate geometry (area > 10 m^2, vertices >= 3)

### Roof Type Heuristics (generate_stl_nyc3d.py)

- Height > 100m: Flat roof
- Construction year < 1950: Gabled roof
- Residential + height < 20m: Gabled roof
- Default: Flat roof

---

## Comparison: NYC Open Data vs OpenStreetMap

| Feature | NYC Open Data | OpenStreetMap |
|---------|--------------|---------------|
| Height Accuracy | ±2 feet (surveyed) | Variable (user-tagged) |
| Coverage | Complete NYC | Global, variable |
| Roof Geometry | Limited (flat tops) | Rich tags (roof:shape) |
| Authentication | None | None |
| Use Case | NYC-specific, high accuracy | Global coverage |

---

## References

- [NYC Open Data Metadata PDF](https://www.nyc.gov/assets/planning/download/pdf/data-maps/open-data/nyc-3d-model-metadata.pdf)
- [NYC Building Footprints Metadata (GitHub)](https://github.com/CityOfNewYork/nyc-geo-metadata/blob/main/Metadata/Metadata_BuildingFootprints.md)
- [TUM GIS 3D City Model Project](https://www.asg.ed.tum.de/en/gis/projects/3d-city-model-of-new-york-city/)
- [CityGML Standard (OGC)](https://www.ogc.org/standards/citygml)
- [Socrata API Documentation](https://dev.socrata.com/docs/queries/)
