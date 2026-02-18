# NYC 3D Building Model Data Access Guide

## Overview

This document describes how to access NYC DOITT (Department of Information Technology and Telecommunications) 3D Building Model data for CFD simulations. The `generate_stl_nyc3d.py` script leverages multiple NYC open data sources to create accurate building geometry for the Citicorp Center CFD case.

---

## Data Sources

### 1. NYC Open Data Building Footprints API (Primary Source)

**Dataset**: BUILDING
**Dataset ID**: `5zhs-2jue`
**Format**: GeoJSON via Socrata API
**URL**: https://data.cityofnewyork.us/City-Government/BUILDING/5zhs-2jue

#### Key Attributes

- `height_roof`: Roof height above ground (feet)
- `groundelev`: Ground elevation at building base (feet)
- `bin`: Building Identification Number (unique per building)
- `cnstrct_yr`: Construction year
- `lstmoddate`: Last modified date
- `lststatype`: Last status type (construction status)
- `the_geom`: Building footprint geometry (MultiPolygon or Polygon)

#### API Access

**Endpoint**: `https://data.cityofnewyork.us/resource/5zhs-2jue.geojson`

**Example Query** (Python with requests):
```python
import requests

params = {
    "$where": "within_box(the_geom, 40.7620, -73.9750, 40.7540, -73.9630)",
    "$limit": 5000,
    "$order": "height_roof DESC"
}

response = requests.get(
    "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson",
    params=params,
    timeout=60
)
data = response.json()
```

**Spatial Query Syntax**:
- `within_box(the_geom, north_lat, west_lon, south_lat, east_lon)`: Buildings within bounding box
- `within_circle(the_geom, lat, lon, radius_meters)`: Buildings within radius

**Authentication**: None required (public dataset)

**Rate Limits**: Socrata API allows ~1000 requests/hour without authentication. For higher limits, register for an app token at https://dev.socrata.com/

#### Data Quality

- **Completeness**: ~1 million buildings across NYC
- **Height Accuracy**: ±2 feet for buildings >60 feet (photogrammetric survey)
- **Update Frequency**: Daily updates by OTI staff
- **Last Survey**: 2014 aerial survey (base dataset)
- **Data Source**: NYC DoITT aerial imagery + NYC Department of Buildings records

---

### 2. NYC 3D Building Model by Community District

**Dataset**: NYC 3D Model by Community District
**Dataset ID**: `u5j4-zxpn`
**Format**: 3DM (Rhinoceros), CityGML (via third-party conversions)
**URL**: https://data.cityofnewyork.us/City-Government/NYC-3D-Model-by-Community-District/u5j4-zxpn

#### Description

Full 3D building massing model based on DOITT's 2014 aerial survey. Uses hybrid CityGML LOD1/LOD2 specification:
- **LOD1**: Simple block models with flat roofs
- **LOD2**: Detailed roof geometry for some buildings

#### Access Methods

1. **Direct Download**: 3DM files by community district (59 districts)
2. **Bulk Download**: Full city model (~4 GB)
3. **API Access**: Limited (primarily for metadata)

#### File Formats

- **3DM**: Native Rhinoceros format, importable to SketchUp, AutoCAD, Blender
- **CityGML**: XML-based 3D city model format (requires conversion)

**Note**: Our script uses the Building Footprints API instead of downloading 3DM files for easier automation.

---

### 3. Enhanced NYC 3D Building Model (Third-Party)

**Project**: TUM GIS Chair + GeorOcket Enhanced Model
**URL**: https://github.com/georocket/new-york-city-model-enhanced
**Format**: CityGML with PLUTO integration

#### Features

- CityGML LOD2 with ~90 semantic attributes per building
- Integrated with NYC PLUTO (Primary Land Use Tax Lot Output)
- Includes land parcels, roads, parks, terrain, water bodies
- Buildings, terrain, and infrastructure with 3D geometries

#### Access

```bash
# Clone the repository
git clone https://github.com/georocket/new-york-city-model-enhanced.git

# Files organized by borough
cd new-york-city-model-enhanced/CityGML/
# Manhattan, Brooklyn, Queens, Bronx, Staten Island subfolders
```

#### Use Cases

- High-fidelity urban simulations
- Semantic queries (e.g., "all residential buildings > 100m")
- Integration with GIS workflows

**Note**: Large dataset (~20 GB for full city), not used in our script for simplicity.

---

## Script Usage

### Basic Usage

```bash
# Generate STL files using NYC Open Data API
python generate_stl_nyc3d.py

# Offline mode (hand-crafted Citicorp only, no surrounding buildings)
python generate_stl_nyc3d.py --offline
```

### Output Files

Generated in `constant/triSurface/`:

1. **citicorp_tower.stl**: Hand-crafted Citicorp tower
   - 47.85m × 47.85m square footprint
   - 45° rotation from cardinal axes
   - 45° slanted roof (244.15m → 278.9m)
   - 22 triangles

2. **citicorp_stilts.stl**: Four stilt columns
   - 7.32m × 7.32m cross-section
   - Height: 0 → 34.75m (10 stories)
   - Positions: Face midpoints at 45° rotation
   - 48 triangles

3. **surroundings_nyc3d.stl**: Surrounding buildings from NYC Open Data
   - Variable geometry based on API response
   - Intelligent roof type inference
   - Excludes Citicorp (uses hand-crafted model)

### Configuration

Edit constants in `generate_stl_nyc3d.py`:

```python
# CFD domain bounds (meters relative to Citicorp center)
DOMAIN_X_MIN = -200
DOMAIN_X_MAX = 520
DOMAIN_Y_MIN = -360
DOMAIN_Y_MAX = 360

# Filtering
MIN_HEIGHT_M = 5.0  # Skip buildings shorter than 5m
CITICORP_PROXIMITY_M = 40  # Radius to identify Citicorp
```

---

## Data Processing Pipeline

### 1. Coordinate Transformation

```
WGS84 (lat, lon) → Local Cartesian (x, y, z)
Origin: Citicorp Center (40.7579°N, 73.9690°W)
X-axis: East
Y-axis: North
Z-axis: Up
```

**Conversion**:
```python
M_PER_DEG_LAT = 111320.0
M_PER_DEG_LON = 111320.0 * cos(CITICORP_LAT)

x = (lon - CITICORP_LON) * M_PER_DEG_LON
y = (lat - CITICORP_LAT) * M_PER_DEG_LAT
```

### 2. Building Footprint Processing

1. **Extract footprint**: Parse GeoJSON MultiPolygon/Polygon
2. **Transform coordinates**: WGS84 → local meters
3. **Ensure CCW winding**: Required for correct normals
4. **Validate geometry**: Check area > 10 m², vertices ≥ 3

### 3. Extrusion with Roof Inference

**Roof Type Heuristics**:
- Height > 100m → Flat roof
- Construction year < 1950 → Gabled roof
- Residential + height < 20m → Gabled roof
- Default → Flat roof

**Geometry Generation**:
- **Flat**: Simple extrusion to `height_roof`
- **Gabled**: Ridge along longest edge (future enhancement)
- **Hipped**: Pyramid roof (future enhancement)

### 4. STL Export

- **Format**: ASCII STL (human-readable, ParaView/MeshLab compatible)
- **Normals**: Computed via cross product (right-hand rule)
- **Watertight**: Guaranteed closed surfaces for meshing

---

## Metadata References

### Official Documentation

1. **NYC Open Data Metadata PDF**:
   https://www.nyc.gov/assets/planning/download/pdf/data-maps/open-data/nyc-3d-model-metadata.pdf

2. **GitHub Metadata Repository**:
   https://github.com/CityOfNewYork/nyc-geo-metadata/blob/main/Metadata/Metadata_BuildingFootprints.md

3. **NYC Planning 3D Model Page**:
   https://www.nyc.gov/site/planning/data-maps/open-data/dwn-nyc-3d-model-download.page

4. **Socrata API Documentation**:
   https://dev.socrata.com/docs/queries/

### Academic References

1. **TUM GIS 3D City Model Project**:
   https://www.asg.ed.tum.de/en/gis/projects/3d-city-model-of-new-york-city/

2. **CityGML Standard** (OGC):
   https://www.ogc.org/standards/citygml

3. **GeorOcket Enhanced Model** (GitHub):
   https://github.com/georocket/new-york-city-model-enhanced

---

## Comparison with OpenStreetMap (OSM)

| Feature | NYC Open Data | OpenStreetMap |
|---------|--------------|---------------|
| **Height Accuracy** | ±2 feet (surveyed) | Variable (user-tagged) |
| **Coverage** | Complete NYC | Global, variable completeness |
| **Update Frequency** | Daily | Continuous (crowdsourced) |
| **Roof Geometry** | Limited (flat tops) | Rich tags (`roof:shape`, `roof:height`) |
| **Authentication** | None | None |
| **API Complexity** | Moderate (Socrata) | High (Overpass QL) |
| **Use Case** | NYC-specific, high accuracy | Global, community-driven |

**Recommendation**: Use NYC Open Data for CFD simulations requiring accurate building heights in NYC. Use OSM for global coverage or detailed roof geometry when available.

---

## Troubleshooting

### API Query Fails

**Symptom**: `API query failed: HTTPError 500`

**Solutions**:
1. Check network connectivity
2. Verify bounding box is valid (NYC coordinates)
3. Reduce `$limit` parameter (try 1000 instead of 5000)
4. Retry after delay (rate limit may be hit)
5. Use `--offline` flag to skip API and use hand-crafted geometry only

### Missing Buildings

**Symptom**: Fewer buildings than expected in `surroundings_nyc3d.stl`

**Causes**:
1. Buildings below `MIN_HEIGHT_M` threshold (default: 5m)
2. Buildings outside CFD domain bounds
3. API result limit reached (increase `$limit` in query)
4. Missing `height_roof` attribute (filtered out)

**Solutions**:
- Lower `MIN_HEIGHT_M` to include shorter buildings
- Expand domain bounds (`DOMAIN_X_MIN`, etc.)
- Increase API `$limit` parameter

### Citicorp Not Identified

**Symptom**: `Warning: Citicorp not auto-identified in API data`

**Impact**: None (hand-crafted model used regardless)

**Cause**: Citicorp outside proximity threshold or BIN mismatch

**Solution**: Update `CITICORP_BIN_LIST` or `CITICORP_PROXIMITY_M`

---

## Future Enhancements

1. **Gabled Roof Geometry**: Implement ridge detection and triangulation
2. **CityGML Import**: Direct parsing of NYC 3D Model CityGML files
3. **PLUTO Integration**: Semantic attributes (land use, year built, etc.)
4. **Binary STL Export**: Faster parsing for large meshes
5. **Parallel API Queries**: Batch processing for large domains
6. **Roof Material Properties**: Assign boundary conditions based on roof type

---

## License and Attribution

### NYC Open Data

- **License**: Public Domain (NYC Open Data Policy)
- **Attribution**: "Data provided by NYC Department of Information Technology and Telecommunications (DoITT)"
- **Liability**: Data provided "as-is" without warranty

### This Script

- **License**: MIT License (same as repository)
- **Author**: Generated for OR 750 Reliability Analysis, GMU PhD Program
- **Date**: 2026-02-16

---

## Support and Contact

**NYC Open Data Support**:
https://opendata.cityofnewyork.us/contact/

**Socrata API Support**:
https://support.socrata.com/

**Repository Issues**:
https://github.com/akirk20-code/citicorp-crisis-simulation/issues

---

**Last Updated**: 2026-02-16
**Script Version**: 1.0
**NYC Data Vintage**: 2014 aerial survey + daily updates
