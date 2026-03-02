# STL Generator Comparison

This document compares the three STL generation scripts available for the Citicorp Center CFD simulation.

---

## Quick Comparison Table

| Feature | `generate_stl.py` | `generate_stl_osm.py` | `generate_stl_nyc3d.py` |
|---------|-------------------|----------------------|------------------------|
| **Primary Data Source** | NYC Open Data Footprints | OpenStreetMap Overpass | NYC Open Data Footprints |
| **Height Data Quality** | Surveyed (±2 ft) | User-tagged (variable) | Surveyed (±2 ft) |
| **Roof Geometry** | Flat only | Flat only | Flat + roof inference |
| **Citicorp Tower** | Hand-crafted | Hand-crafted (simplified) | Hand-crafted (detailed) |
| **Citicorp Rotation** | No (0°) | Yes (45°) | Yes (45°) |
| **Citicorp Roof** | Flat | Simplified slant | 45° slant (accurate) |
| **Stilt Positions** | Face midpoints | Corners | Face midpoints (accurate) |
| **Dependencies** | requests | requests (optional) | requests (optional) |
| **Fallback Mode** | 18 hardcoded buildings | None (OSM only) | None (API or offline) |
| **Output Files** | `surroundings.stl` | `surroundings_osm.stl` | `surroundings_nyc3d.stl` |
| **Best For** | Quick testing | Global applicability | Production CFD runs |

---

## Detailed Comparison

### 1. `generate_stl.py` (Original)

**File**: `c:\Users\kirka\OneDrive\Documents\GMU PhD\OR 750 Reliability, Safety, and Risk\Citi Sim\citicorp_cfd\generate_stl.py`

**Data Source**: NYC Open Data Building Footprints API (Socrata)

**Citicorp Geometry**:
- Tower: Simple box (no roof slant, no rotation)
- Dimensions: 47.85m × 47.85m × 278.9m
- Stilts: 4 columns at face midpoints (correct per Morgenstern 1995)
- Rotation: 0° (not rotated)

**Surrounding Buildings**:
- Fetches from NYC Open Data API
- Extrudes to `height_roof` (flat roofs only)
- Fallback: 18 hardcoded buildings if API fails

**Strengths**:
- Reliable fallback mode
- Fast execution
- Proven geometry for meshing

**Weaknesses**:
- Citicorp tower lacks slanted roof
- No 45° rotation (doesn't match real building orientation)
- Simple flat roofs for all buildings

**Use Case**: Quick CFD testing, prototyping mesh parameters

---

### 2. `generate_stl_osm.py` (OpenStreetMap)

**File**: `c:\Users\kirka\OneDrive\Documents\GMU PhD\OR 750 Reliability, Safety, and Risk\Citi Sim\citicorp_cfd\generate_stl_osm.py`

**Data Source**: OpenStreetMap Overpass API

**Citicorp Geometry**:
- Tower: 60m × 60m (approximate, not exact)
- Roof: Simplified slant (single plane)
- Stilts: 4 columns at corners (8m × 8m)
- Rotation: 45° (correct orientation)

**Surrounding Buildings**:
- Fetches from OpenStreetMap
- Height estimation from OSM tags:
  - `height`: Direct meter value
  - `building:levels`: Levels × 3.5m
  - `building:height`: Alternative height tag
  - Default: 15m (4-5 stories)
- Extrudes to estimated height (flat roofs)

**Strengths**:
- Global coverage (works anywhere, not just NYC)
- Rich OSM tagging (building types, names, etc.)
- 45° tower rotation
- No API rate limits (public Overpass instances)

**Weaknesses**:
- Height data quality varies (user-tagged)
- Citicorp dimensions approximate (60m vs 47.85m)
- Stilt positions incorrect (corners vs face midpoints)
- No fallback if Overpass API fails

**Use Case**: Global CFD studies, exploratory analysis, OSM data validation

---

### 3. `generate_stl_nyc3d.py` (NYC 3D Model) **[RECOMMENDED]**

**File**: `c:\Users\kirka\OneDrive\Documents\GMU PhD\OR 750 Reliability, Safety, and Risk\Citi Sim\citicorp_cfd\generate_stl_nyc3d.py`

**Data Source**: NYC Open Data Building Footprints API (Socrata)

**Citicorp Geometry**:
- Tower: Exact dimensions (47.85m × 47.85m)
- Roof: Accurate 45° slant (244.15m → 278.9m)
- Stilts: 4 columns at face midpoints (7.32m × 7.32m)
- Rotation: 45° from cardinal axes (correct per Morgenstern 1995)
- Total height: 278.9m (915 ft)

**Surrounding Buildings**:
- Fetches from NYC Open Data API
- Uses surveyed `height_roof` data (±2 ft accuracy)
- Intelligent roof type inference:
  - Tall buildings (>100m): Flat
  - Old buildings (<1950): Gabled (geometry TBD)
  - Low-rise residential: Gabled
  - Default: Flat
- Excludes Citicorp by BIN or proximity

**Strengths**:
- **Most accurate Citicorp geometry** (matches Morgenstern 1995)
- Surveyed height data (±2 ft accuracy)
- Intelligent roof inference
- 45° tower rotation (correct orientation)
- Correct stilt positions (face midpoints)
- Comprehensive documentation

**Weaknesses**:
- NYC-only (not globally applicable)
- No offline fallback for surrounding buildings
- Roof inference heuristics (not perfect)

**Use Case**: **Production CFD simulations**, research publications, high-fidelity analysis

---

## Geometry Verification

### Citicorp Tower Dimensions

| Parameter | Real Building | `generate_stl.py` | `generate_stl_osm.py` | `generate_stl_nyc3d.py` |
|-----------|---------------|-------------------|-----------------------|------------------------|
| Footprint | 157 ft (47.85m) | 47.85m ✓ | ~60m ✗ | 47.85m ✓ |
| Total Height | 915 ft (278.9m) | 278.9m ✓ | 279m ✓ | 278.9m ✓ |
| Roof Slant | 45° (S→N) | None ✗ | Simplified | 45° (244.15→278.9m) ✓ |
| Rotation | 45° | 0° ✗ | 45° ✓ | 45° ✓ |
| Stilt Height | 114 ft (34.75m) | 34.75m ✓ | 35m ✓ | 34.75m ✓ |
| Stilt Width | 24 ft (7.32m) | 7.32m ✓ | 8m ✗ | 7.32m ✓ |
| Stilt Positions | Face midpoints | Correct ✓ | Corners ✗ | Correct ✓ |

### Morgenstern (1995) Specifications

Per the Morgenstern structural analysis paper:
- **Tower**: 157 ft × 157 ft square plan, 915 ft total height
- **Orientation**: 45° rotation from cardinal axes (diamond orientation)
- **Roof**: 45° slant from south (801 ft) to north (915 ft)
- **Stilts**: 10 stories (114 ft), 24 ft × 24 ft, at **face midpoints** (not corners)

**Winner**: `generate_stl_nyc3d.py` matches all specifications ✓

---

## Data Source Comparison

### NYC Open Data (used by `generate_stl.py` and `generate_stl_nyc3d.py`)

**Source**: NYC Department of Information Technology and Telecommunications (DoITT)

**Data Collection**:
- 2014 aerial survey (photogrammetric)
- NYC Department of Buildings records
- Daily updates from construction permits

**Accuracy**:
- Height: ±2 feet for buildings >60 feet
- Footprint: ±1 meter (sub-pixel accuracy)
- Coverage: ~1 million buildings (complete NYC)

**Attributes**:
- `height_roof`: Surveyed roof height (feet)
- `groundelev`: Ground elevation (feet)
- `bin`: Building Identification Number
- `cnstrct_yr`: Construction year
- `lststatype`: Construction status

**Advantages**:
- High accuracy (professional survey)
- Complete coverage (all NYC buildings)
- Authoritative source (city government)
- Daily updates

**Disadvantages**:
- NYC-only (no global coverage)
- Limited roof geometry (mostly flat tops)
- API rate limits (1000 req/hr without token)

---

### OpenStreetMap (used by `generate_stl_osm.py`)

**Source**: OpenStreetMap community (crowdsourced)

**Data Collection**:
- User contributions (GPS traces, digitization)
- Import from authoritative sources (varies by region)
- Continuous updates (anyone can edit)

**Accuracy**:
- Variable (depends on contributor knowledge)
- Height: Often estimated or missing
- Footprint: Generally good in urban areas
- Coverage: Global, but variable completeness

**Attributes**:
- `height`: Direct height value (meters)
- `building:levels`: Number of stories
- `building:height`: Alternative height tag
- `roof:shape`: Roof type (flat, gabled, hipped, etc.)
- `roof:height`: Roof height (meters)
- `name`: Building name

**Advantages**:
- Global coverage
- Rich semantic tags
- No API authentication
- Community-maintained (always improving)

**Disadvantages**:
- Variable accuracy (user-tagged)
- Missing height data (common)
- Inconsistent tagging schemes
- No official authority

---

## Recommendation

### For Production CFD Simulations → `generate_stl_nyc3d.py`

**Reasons**:
1. Most accurate Citicorp geometry (matches Morgenstern 1995)
2. Surveyed building heights (±2 ft accuracy)
3. Correct 45° tower rotation
4. Accurate 45° slanted roof
5. Correct stilt positions (face midpoints)
6. Intelligent roof type inference

**Command**:
```bash
python generate_stl_nyc3d.py
```

**Output**: `citicorp_tower.stl`, `citicorp_stilts.stl`, `surroundings_nyc3d.stl`

---

### For Quick Testing → `generate_stl.py`

**Reasons**:
1. Reliable fallback mode (18 hardcoded buildings)
2. Fast execution
3. Proven geometry for meshing

**Command**:
```bash
python generate_stl.py
```

**Output**: `citicorp_tower.stl`, `citicorp_stilts.stl`, `surroundings.stl`

---

### For Global Studies → `generate_stl_osm.py`

**Reasons**:
1. Works anywhere (not NYC-specific)
2. Rich OSM tagging
3. No API authentication

**Command**:
```bash
python generate_stl_osm.py
```

**Output**: `citicorp_tower.stl`, `citicorp_stilts.stl`, `surroundings_osm.stl`

---

## Migration Guide

If you're currently using `generate_stl.py` or `generate_stl_osm.py`, here's how to migrate to `generate_stl_nyc3d.py`:

### Step 1: Backup Existing Files

```bash
cd "c:\Users\kirka\OneDrive\Documents\GMU PhD\OR 750 Reliability, Safety, and Risk\Citi Sim\citicorp_cfd"
cp constant/triSurface/citicorp_tower.stl constant/triSurface/citicorp_tower.stl.bak
cp constant/triSurface/citicorp_stilts.stl constant/triSurface/citicorp_stilts.stl.bak
cp constant/triSurface/surroundings.stl constant/triSurface/surroundings.stl.bak
```

### Step 2: Generate New STL Files

```bash
python generate_stl_nyc3d.py
```

### Step 3: Update `snappyHexMeshDict`

If you're using `surroundings.stl`, update to `surroundings_nyc3d.stl`:

```cpp
// Before
surroundings.stl
{
    type triSurfaceMesh;
    name surroundings;
}

// After
surroundings_nyc3d.stl
{
    type triSurfaceMesh;
    name surroundings;
}
```

### Step 4: Re-run Meshing

```bash
# In WSL2 Ubuntu
cd /mnt/c/Users/kirka/OneDrive/Documents/GMU\ PhD/OR\ 750\ Reliability,\ Safety,\ and\ Risk/Citi\ Sim/citicorp_cfd
source /opt/openfoam13/etc/bashrc
blockMesh
snappyHexMesh -overwrite
```

### Step 5: Verify Results

Compare mesh quality metrics:
- Check surface coverage with `checkMesh`
- Verify tower rotation in ParaView
- Compare flow patterns (if available)

---

## Summary

**Choose `generate_stl_nyc3d.py`** for production CFD work requiring:
- Accurate Citicorp geometry (Morgenstern 1995 specs)
- High-quality building height data (NYC surveyed)
- 45° tower rotation (correct orientation)
- Detailed slanted roof (45° slope)

**Use `generate_stl.py`** for quick prototyping with reliable fallback mode.

**Use `generate_stl_osm.py`** for global CFD studies outside NYC.

---

**Last Updated**: 2026-02-16
**Recommended**: `generate_stl_nyc3d.py` (v1.0)
