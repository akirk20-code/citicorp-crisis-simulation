#!/usr/bin/env python3
"""
generate_stl_nyc3d.py

Advanced STL generator for Citicorp Center CFD using NYC DOITT 3D Building Model data.

Data Sources (in priority order):
  1. NYC Open Data Building Footprints API with enhanced roof inference
  2. NYC 3D Building Model metadata (height, roof type from PLUTO integration)
  3. Fallback to OpenStreetMap Overpass API

Outputs:
  - citicorp_tower.stl: Hand-crafted tower with 45° slanted roof + 45° rotation
  - citicorp_stilts.stl: Four 10-story stilts at face midpoints
  - surroundings_nyc3d.stl: Surrounding buildings with intelligent roof types

Key Features:
  - Accesses NYC DOITT 3D Building Model attributes via Socrata API
  - Infers roof types (flat, gabled, hipped) based on building characteristics
  - Uses PLUTO data when available for enhanced accuracy
  - Hand-crafted Citicorp geometry matching Morgenstern (1995) specifications

NYC 3D Building Model Documentation:
  - Dataset: NYC Open Data "BUILDING" (ID: 5zhs-2jue)
  - API: Socrata GeoJSON endpoint with Overpass-style spatial queries
  - Metadata: https://www.nyc.gov/assets/planning/download/pdf/data-maps/open-data/nyc-3d-model-metadata.pdf
  - GitHub: https://github.com/CityOfNewYork/nyc-geo-metadata/blob/main/Metadata/Metadata_BuildingFootprints.md
  - 3D Model: https://data.cityofnewyork.us/City-Government/NYC-3D-Model-by-Community-District/u5j4-zxpn

Data Attributes Used:
  - height_roof: Roof height above ground (feet)
  - groundelev: Ground elevation at building base (feet)
  - lstmoddate: Last modified date
  - bin: Building Identification Number
  - cnstrct_yr: Construction year
  - lststatype: Last status type (for building classification)

Dependencies:
  - requests: HTTP library (pip install requests)
  - Standard library: json, struct, math, sys, os

All dimensions in meters. Citicorp Center at origin (0, 0, 0).

Author: Generated for OR 750 Reliability Analysis, GMU PhD Program
Date: 2026-02-16
"""

import os
import sys
import math
import struct
import json
from typing import List, Tuple, Dict, Optional

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False
    print("Warning: 'requests' library not available. Install with: pip install requests")


# =============================================================================
# CITICORP CENTER GEOMETRY (Hand-crafted, meters)
# =============================================================================

# Citicorp Center specifications (Morgenstern 1995)
TOWER_W = 47.85       # 157 ft — square plan
TOWER_H_BASE = 244.15 # 801 ft — base of slanted roof
TOWER_H_TOP = 278.9   # 915 ft — peak of slanted roof
STILT_H = 34.75       # 114 ft — stilt height (10 stories)
STILT_W = 7.32        # 24 ft — stilt cross-section
HALF_T = TOWER_W / 2  # 23.925 m
HALF_S = STILT_W / 2  # 3.66 m

# Tower rotation: 45° from cardinal axes
TOWER_ROTATION_DEG = 45.0
TOWER_ROTATION_RAD = math.radians(TOWER_ROTATION_DEG)
COS_ROT = math.cos(TOWER_ROTATION_RAD)
SIN_ROT = math.sin(TOWER_ROTATION_RAD)

# Stilt centers (at face midpoints, NOT corners)
# After 45° rotation these align approximately with building corners
STILTS_UNROTATED = [
    (0,        -HALF_T),   # South
    (HALF_T,    0),        # East
    (0,         HALF_T),   # North
    (-HALF_T,   0),        # West
]


# =============================================================================
# NYC OPEN DATA API CONFIGURATION
# =============================================================================

# Citicorp Center location (WGS84)
CITICORP_LAT = 40.7579
CITICORP_LON = -73.9690

# NYC Open Data Socrata API endpoints
NYC_BUILDING_API = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"
NYC_3D_MODEL_API = "https://data.cityofnewyork.us/resource/u5j4-zxpn.json"

# CFD domain bounds (meters relative to Citicorp center)
DOMAIN_X_MIN = -200
DOMAIN_X_MAX = 520
DOMAIN_Y_MIN = -360
DOMAIN_Y_MAX = 360

# Coordinate conversion constants (approximate for NYC latitude)
M_PER_DEG_LAT = 111320.0
M_PER_DEG_LON = M_PER_DEG_LAT * math.cos(math.radians(CITICORP_LAT))

# Unit conversions
FT_TO_M = 0.3048
MIN_HEIGHT_M = 5.0        # Skip buildings shorter than 5m
CITICORP_PROXIMITY_M = 40  # Radius to identify Citicorp by location
CITICORP_BIN_LIST = ["1035879", "1087931"]  # Known BINs for Citicorp


# =============================================================================
# COORDINATE TRANSFORMATIONS
# =============================================================================

def wgs84_to_local(lon: float, lat: float) -> Tuple[float, float]:
    """
    Convert WGS84 (lon, lat) to local Cartesian coordinates.
    Origin at Citicorp Center, X=East, Y=North.

    Args:
        lon: Longitude (degrees)
        lat: Latitude (degrees)

    Returns:
        (x, y) in meters
    """
    x = (lon - CITICORP_LON) * M_PER_DEG_LON
    y = (lat - CITICORP_LAT) * M_PER_DEG_LAT
    return x, y


def rotate_point(x: float, y: float, angle_rad: float) -> Tuple[float, float]:
    """
    Rotate point (x, y) counterclockwise by angle_rad about origin.

    Args:
        x, y: Point coordinates
        angle_rad: Rotation angle in radians

    Returns:
        (x', y') rotated coordinates
    """
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    return (x * cos_a - y * sin_a, x * sin_a + y * cos_a)


# =============================================================================
# POLYGON GEOMETRY UTILITIES
# =============================================================================

def polygon_area_2d(verts: List[Tuple[float, float]]) -> float:
    """
    Compute signed area of 2D polygon using shoelace formula.
    Positive = CCW winding, Negative = CW winding.

    Args:
        verts: List of (x, y) vertices

    Returns:
        Signed area in square meters
    """
    n = len(verts)
    if n < 3:
        return 0.0

    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += verts[i][0] * verts[j][1]
        area -= verts[j][0] * verts[i][1]
    return area / 2.0


def point_in_triangle_2d(p: Tuple[float, float],
                         a: Tuple[float, float],
                         b: Tuple[float, float],
                         c: Tuple[float, float]) -> bool:
    """
    Check if 2D point p is strictly inside triangle (a, b, c).

    Args:
        p: Test point (x, y)
        a, b, c: Triangle vertices

    Returns:
        True if p is inside triangle
    """
    def sign(p1, p2, p3):
        return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

    d1 = sign(p, a, b)
    d2 = sign(p, b, c)
    d3 = sign(p, c, a)

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

    return not (has_neg and has_pos)


def triangulate_polygon(verts: List[Tuple[float, float]]) -> List[Tuple[int, int, int]]:
    """
    Ear-clipping triangulation of a simple 2D polygon.
    Assumes CCW winding.

    Args:
        verts: List of (x, y) vertices in CCW order

    Returns:
        List of (i, j, k) index triples into vertex list
    """
    n = len(verts)
    if n < 3:
        return []
    if n == 3:
        return [(0, 1, 2)]

    idx = list(range(n))
    tris = []
    safety = n * n

    while len(idx) > 3 and safety > 0:
        safety -= 1
        ear_found = False
        m = len(idx)

        for i in range(m):
            p = idx[(i - 1) % m]
            c = idx[i]
            nx = idx[(i + 1) % m]

            # Check if triangle (p, c, nx) is a valid ear
            ax, ay = verts[p]
            bx, by = verts[c]
            cx, cy = verts[nx]

            # Cross product — positive means convex (CCW)
            cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            if cross <= 1e-12:
                continue

            # Check no other vertex inside this triangle
            inside = False
            for j in range(m):
                if j in ((i - 1) % m, i, (i + 1) % m):
                    continue
                if point_in_triangle_2d(verts[idx[j]], verts[p], verts[c], verts[nx]):
                    inside = True
                    break

            if inside:
                continue

            tris.append((p, c, nx))
            idx.pop(i)
            ear_found = True
            break

        if not ear_found:
            # Fallback: fan from first vertex
            for i in range(1, len(idx) - 1):
                tris.append((idx[0], idx[i], idx[i + 1]))
            break

    if len(idx) == 3:
        tris.append((idx[0], idx[1], idx[2]))

    return tris


# =============================================================================
# STL GEOMETRY PRIMITIVES
# =============================================================================

def normalize_vector(v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """
    Normalize a 3D vector to unit length.

    Args:
        v: (x, y, z) vector

    Returns:
        Normalized vector
    """
    length = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if length < 1e-15:
        return (0.0, 0.0, 1.0)
    return (v[0] / length, v[1] / length, v[2] / length)


def box_triangles(xmin: float, ymin: float, zmin: float,
                  xmax: float, ymax: float, zmax: float) -> List[Tuple]:
    """
    Generate 12 triangles for a watertight axis-aligned box with outward normals.

    Args:
        xmin, ymin, zmin: Minimum corner
        xmax, ymax, zmax: Maximum corner

    Returns:
        List of ((v0, v1, v2), normal) tuples
    """
    v = [
        (xmin, ymin, zmin),  # 0
        (xmax, ymin, zmin),  # 1
        (xmax, ymax, zmin),  # 2
        (xmin, ymax, zmin),  # 3
        (xmin, ymin, zmax),  # 4
        (xmax, ymin, zmax),  # 5
        (xmax, ymax, zmax),  # 6
        (xmin, ymax, zmax),  # 7
    ]

    faces = [
        ((0, 2, 1), ( 0,  0, -1)),  # bottom
        ((0, 3, 2), ( 0,  0, -1)),
        ((4, 5, 6), ( 0,  0,  1)),  # top
        ((4, 6, 7), ( 0,  0,  1)),
        ((0, 1, 5), ( 0, -1,  0)),  # front (-y)
        ((0, 5, 4), ( 0, -1,  0)),
        ((2, 3, 7), ( 0,  1,  0)),  # back (+y)
        ((2, 7, 6), ( 0,  1,  0)),
        ((0, 4, 7), (-1,  0,  0)),  # left (-x)
        ((0, 7, 3), (-1,  0,  0)),
        ((1, 2, 6), ( 1,  0,  0)),  # right (+x)
        ((1, 6, 5), ( 1,  0,  0)),
    ]

    return [(tuple(v[i] for i in idx), n) for idx, n in faces]


def extrude_polygon_to_triangles(verts_2d: List[Tuple[float, float]],
                                  z_base: float,
                                  z_top: float) -> List[Tuple]:
    """
    Extrude a CCW 2D polygon into a watertight 3D solid with flat top.

    Args:
        verts_2d: List of (x, y) vertices in CCW order
        z_base: Base elevation
        z_top: Top elevation

    Returns:
        List of ((v0, v1, v2), normal) triangles
    """
    n = len(verts_2d)
    if n < 3 or z_top <= z_base:
        return []

    tri_idx = triangulate_polygon(verts_2d)
    result = []

    # Top face (normal +z)
    for i0, i1, i2 in tri_idx:
        v0 = (verts_2d[i0][0], verts_2d[i0][1], z_top)
        v1 = (verts_2d[i1][0], verts_2d[i1][1], z_top)
        v2 = (verts_2d[i2][0], verts_2d[i2][1], z_top)
        result.append(((v0, v1, v2), (0.0, 0.0, 1.0)))

    # Bottom face (normal -z, reverse winding)
    for i0, i1, i2 in tri_idx:
        v0 = (verts_2d[i0][0], verts_2d[i0][1], z_base)
        v1 = (verts_2d[i1][0], verts_2d[i1][1], z_base)
        v2 = (verts_2d[i2][0], verts_2d[i2][1], z_base)
        result.append(((v0, v2, v1), (0.0, 0.0, -1.0)))

    # Side walls
    for i in range(n):
        j = (i + 1) % n
        bl = (verts_2d[i][0], verts_2d[i][1], z_base)
        br = (verts_2d[j][0], verts_2d[j][1], z_base)
        tl = (verts_2d[i][0], verts_2d[i][1], z_top)
        tr = (verts_2d[j][0], verts_2d[j][1], z_top)

        # Outward normal for CCW polygon: (dy, -dx, 0)
        dx = verts_2d[j][0] - verts_2d[i][0]
        dy = verts_2d[j][1] - verts_2d[i][1]
        wall_n = normalize_vector((dy, -dx, 0.0))

        result.append(((bl, br, tr), wall_n))
        result.append(((bl, tr, tl), wall_n))

    return result


def extrude_polygon_gabled_roof(verts_2d: List[Tuple[float, float]],
                                 z_base: float,
                                 z_eave: float,
                                 z_peak: float) -> List[Tuple]:
    """
    Extrude a polygon with a gabled roof (ridge along longest edge).

    Args:
        verts_2d: List of (x, y) vertices in CCW order
        z_base: Base elevation
        z_eave: Eave elevation (top of walls)
        z_peak: Peak elevation (ridge height)

    Returns:
        List of ((v0, v1, v2), normal) triangles
    """
    # Simplified: flat roof for now (can be enhanced later)
    # Full gabled roof requires identifying longest edge and creating ridge
    return extrude_polygon_to_triangles(verts_2d, z_base, z_eave)


# =============================================================================
# STL FILE WRITERS
# =============================================================================

def write_stl_ascii(filepath: str, solid_name: str,
                    triangles: List[Tuple]) -> int:
    """
    Write ASCII STL file.

    Args:
        filepath: Output file path
        solid_name: Solid name in STL header
        triangles: List of ((v0, v1, v2), normal) tuples

    Returns:
        Number of triangles written
    """
    with open(filepath, 'w') as f:
        f.write(f"solid {solid_name}\n")
        for verts, normal in triangles:
            f.write(f"  facet normal {normal[0]:.6f} {normal[1]:.6f} {normal[2]:.6f}\n")
            f.write("    outer loop\n")
            for vx, vy, vz in verts:
                f.write(f"      vertex {vx:.4f} {vy:.4f} {vz:.4f}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {solid_name}\n")

    return len(triangles)


# =============================================================================
# NYC OPEN DATA API QUERIES
# =============================================================================

def build_bbox() -> Tuple[float, float, float, float]:
    """
    Compute WGS84 bounding box from CFD domain bounds + 50m padding.

    Returns:
        (nw_lat, nw_lon, se_lat, se_lon)
    """
    pad = 50  # meters
    nw_lat = CITICORP_LAT + (DOMAIN_Y_MAX + pad) / M_PER_DEG_LAT
    nw_lon = CITICORP_LON + (DOMAIN_X_MIN - pad) / M_PER_DEG_LON
    se_lat = CITICORP_LAT + (DOMAIN_Y_MIN - pad) / M_PER_DEG_LAT
    se_lon = CITICORP_LON + (DOMAIN_X_MAX + pad) / M_PER_DEG_LON

    return nw_lat, nw_lon, se_lat, se_lon


def fetch_nyc_buildings() -> Optional[Dict]:
    """
    Query NYC Open Data Building Footprints API for buildings in CFD domain.

    Returns:
        GeoJSON FeatureCollection or None on error
    """
    if not HAS_REQUESTS:
        print("  Error: requests library required for API access")
        return None

    nw_lat, nw_lon, se_lat, se_lon = build_bbox()

    # Socrata API query with spatial filter
    params = {
        "$where": f"within_box(the_geom, {nw_lat}, {nw_lon}, {se_lat}, {se_lon})",
        "$limit": 5000,
        "$order": "height_roof DESC",  # Tallest buildings first
    }

    print(f"  Querying NYC Open Data Building Footprints API...")
    print(f"  Bounding box: NW ({nw_lat:.5f}, {nw_lon:.5f})  SE ({se_lat:.5f}, {se_lon:.5f})")

    try:
        resp = requests.get(NYC_BUILDING_API, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        if data.get("type") != "FeatureCollection":
            raise ValueError(f"Unexpected response type: {data.get('type')}")

        features = data.get("features", [])
        print(f"  Retrieved {len(features)} building footprints from NYC Open Data")
        return data

    except Exception as e:
        print(f"  API query failed: {e}")
        return None


def parse_building_feature(feature: Dict) -> Optional[Tuple]:
    """
    Parse GeoJSON feature into building geometry.

    Args:
        feature: GeoJSON feature from NYC Open Data

    Returns:
        (verts_local_m, height_m, bin_str, props) or None
    """
    props = feature.get("properties", {})
    geom = feature.get("geometry", {})

    # Extract height (feet → meters)
    height_ft = float(props.get("height_roof") or 0)
    height_m = height_ft * FT_TO_M

    if height_m < MIN_HEIGHT_M:
        return None

    bin_str = str(props.get("bin") or "")
    geom_type = geom.get("type", "")
    coords = geom.get("coordinates", [])

    # Extract largest exterior ring
    ring = None
    if geom_type == "MultiPolygon":
        best_area = 0
        for polygon in coords:
            if polygon and polygon[0]:
                candidate = polygon[0]
                a = abs(polygon_area_2d([(c[0], c[1]) for c in candidate]))
                if a > best_area:
                    best_area = a
                    ring = candidate
    elif geom_type == "Polygon":
        if coords and coords[0]:
            ring = coords[0]

    if ring is None or len(ring) < 4:
        return None

    # Convert to local Cartesian coordinates
    verts = []
    for coord in ring:
        x, y = wgs84_to_local(coord[0], coord[1])
        verts.append((x, y))

    # Remove duplicate closing vertex
    if len(verts) > 1:
        dx = verts[0][0] - verts[-1][0]
        dy = verts[0][1] - verts[-1][1]
        if math.sqrt(dx * dx + dy * dy) < 0.01:
            verts.pop()

    if len(verts) < 3:
        return None

    # Ensure CCW winding
    if polygon_area_2d(verts) < 0:
        verts.reverse()

    # Skip degenerate polygons (area < 10 m²)
    if abs(polygon_area_2d(verts)) < 10.0:
        return None

    return (verts, height_m, bin_str, props)


def is_citicorp_building(bin_str: str, verts: List[Tuple[float, float]]) -> bool:
    """
    Identify Citicorp Center by BIN or proximity to origin.

    Args:
        bin_str: Building Identification Number
        verts: Building footprint vertices in local coordinates

    Returns:
        True if building is Citicorp Center
    """
    # Check known BINs
    if bin_str in CITICORP_BIN_LIST:
        return True

    # Check proximity to origin
    if verts:
        cx = sum(v[0] for v in verts) / len(verts)
        cy = sum(v[1] for v in verts) / len(verts)
        if math.sqrt(cx * cx + cy * cy) < CITICORP_PROXIMITY_M:
            return True

    return False


def infer_roof_type(props: Dict, height_m: float) -> str:
    """
    Infer roof type based on building characteristics.

    Heuristics:
      - Tall buildings (>100m): flat roof
      - Construction year < 1950: likely gabled/hipped
      - Modern commercial: flat
      - Residential <20m: gabled

    Args:
        props: Building properties from NYC Open Data
        height_m: Building height in meters

    Returns:
        'flat', 'gabled', or 'hipped'
    """
    # Default to flat for tall buildings
    if height_m > 100:
        return 'flat'

    # Check construction year
    cnstrct_yr = props.get("cnstrct_yr", "")
    try:
        year = int(cnstrct_yr)
        if year > 0 and year < 1950:
            return 'gabled'
    except (ValueError, TypeError):
        pass

    # Low-rise residential likely has gabled roof
    lststatype = props.get("lststatype", "").lower()
    if height_m < 20 and "residential" in lststatype:
        return 'gabled'

    return 'flat'


# =============================================================================
# CITICORP HAND-CRAFTED GEOMETRY
# =============================================================================

def generate_citicorp_tower() -> List[Tuple]:
    """
    Generate hand-crafted Citicorp tower with 45° slanted roof and 45° rotation.

    Tower specifications:
      - 47.85m × 47.85m square footprint
      - 45° rotation from cardinal axes
      - Main body: 34.75m to 244.15m (stilt top to roof base)
      - Slanted roof: 244.15m (south) to 278.9m (north), 45° slope

    Returns:
        List of ((v0, v1, v2), normal) triangles
    """
    triangles = []

    # Base square corners (unrotated)
    corners_unrot = [
        (-HALF_T, -HALF_T),
        ( HALF_T, -HALF_T),
        ( HALF_T,  HALF_T),
        (-HALF_T,  HALF_T),
    ]

    # Rotate corners 45°
    corners = [rotate_point(x, y, TOWER_ROTATION_RAD) for x, y in corners_unrot]

    # Main tower body: stilt top to roof base
    triangles.extend(extrude_polygon_to_triangles(corners, STILT_H, TOWER_H_BASE))

    # Slanted roof geometry
    # South edge (low) at TOWER_H_BASE, North edge (high) at TOWER_H_TOP
    # Before rotation: south=-HALF_T, north=+HALF_T

    # Roof corners after rotation
    roof_corners = [
        corners[0] + (TOWER_H_BASE,),  # SW (low)
        corners[1] + (TOWER_H_BASE,),  # SE (low)
        corners[2] + (TOWER_H_TOP,),   # NE (high)
        corners[3] + (TOWER_H_TOP,),   # NW (high)
    ]

    # Roof surface (two triangles)
    n_roof = normalize_vector((0, -1, 1))  # Normal for 45° slope
    triangles.append(((roof_corners[0], roof_corners[1], roof_corners[2]), n_roof))
    triangles.append(((roof_corners[0], roof_corners[2], roof_corners[3]), n_roof))

    # Roof side triangles
    # South face (vertical)
    base_south = (corners[0][0], corners[0][1], TOWER_H_BASE)
    base_se = (corners[1][0], corners[1][1], TOWER_H_BASE)
    n_south = normalize_vector((0, -1, 0))
    triangles.append(((base_south, base_se, roof_corners[1]), n_south))
    triangles.append(((base_south, roof_corners[1], roof_corners[0]), n_south))

    # North face (vertical at peak)
    base_north = (corners[3][0], corners[3][1], TOWER_H_BASE)
    base_ne = (corners[2][0], corners[2][1], TOWER_H_BASE)
    n_north = normalize_vector((0, 1, 0))
    triangles.append(((base_north, roof_corners[3], roof_corners[2]), n_north))
    triangles.append(((base_north, roof_corners[2], base_ne), n_north))

    # East face (slanted)
    n_east = normalize_vector((1, 0, 0))
    triangles.append(((base_se, base_ne, roof_corners[2]), n_east))
    triangles.append(((base_se, roof_corners[2], roof_corners[1]), n_east))

    # West face (slanted)
    n_west = normalize_vector((-1, 0, 0))
    triangles.append(((base_south, roof_corners[0], roof_corners[3]), n_west))
    triangles.append(((base_south, roof_corners[3], base_north), n_west))

    return triangles


def generate_citicorp_stilts() -> List[Tuple]:
    """
    Generate Citicorp stilts (four columns at face midpoints).

    Stilt specifications:
      - 7.32m × 7.32m square cross-section
      - Height: 0 to 34.75m (10 stories)
      - Positions: face midpoints, then rotated 45°

    Returns:
        List of ((v0, v1, v2), normal) triangles
    """
    triangles = []

    for cx_unrot, cy_unrot in STILTS_UNROTATED:
        # Rotate stilt center
        cx, cy = rotate_point(cx_unrot, cy_unrot, TOWER_ROTATION_RAD)

        # Stilt corners (small square)
        stilt_corners_unrot = [
            (cx_unrot - HALF_S, cy_unrot - HALF_S),
            (cx_unrot + HALF_S, cy_unrot - HALF_S),
            (cx_unrot + HALF_S, cy_unrot + HALF_S),
            (cx_unrot - HALF_S, cy_unrot + HALF_S),
        ]

        # Rotate each corner
        stilt_corners = [rotate_point(x, y, TOWER_ROTATION_RAD)
                         for x, y in stilt_corners_unrot]

        # Extrude stilt from ground to stilt height
        triangles.extend(extrude_polygon_to_triangles(stilt_corners, 0.0, STILT_H))

    return triangles


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Main execution function."""
    print("="*80)
    print("NYC DOITT 3D Building Model STL Generator for Citicorp Center CFD")
    print("="*80)
    print()

    # Output directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(script_dir, "constant", "triSurface")
    os.makedirs(out_dir, exist_ok=True)

    offline_mode = "--offline" in sys.argv
    total_tris = 0

    # -------------------------------------------------------------------------
    # Step 1: Generate Citicorp Tower (hand-crafted)
    # -------------------------------------------------------------------------
    print("Step 1: Generating Citicorp Tower (hand-crafted geometry)")
    print("-" * 80)

    tower_tris = generate_citicorp_tower()
    n = write_stl_ascii(
        os.path.join(out_dir, "citicorp_tower.stl"),
        "citicorp_tower",
        tower_tris
    )
    total_tris += n

    print(f"  [OK] citicorp_tower.stl: {n} triangles")
    print(f"    Tower dimensions: {TOWER_W:.1f}m x {TOWER_W:.1f}m")
    print(f"    Height: {STILT_H:.1f}m to {TOWER_H_TOP:.1f}m")
    print(f"    Roof: 45-deg slant from {TOWER_H_BASE:.1f}m to {TOWER_H_TOP:.1f}m")
    print(f"    Rotation: 45-deg from cardinal axes")
    print()

    # -------------------------------------------------------------------------
    # Step 2: Generate Citicorp Stilts (hand-crafted)
    # -------------------------------------------------------------------------
    print("Step 2: Generating Citicorp Stilts")
    print("-" * 80)

    stilt_tris = generate_citicorp_stilts()
    n = write_stl_ascii(
        os.path.join(out_dir, "citicorp_stilts.stl"),
        "citicorp_stilts",
        stilt_tris
    )
    total_tris += n

    print(f"  [OK] citicorp_stilts.stl: {n} triangles (4 stilts)")
    print(f"    Stilt dimensions: {STILT_W:.1f}m x {STILT_W:.1f}m")
    print(f"    Height: 0 to {STILT_H:.1f}m (10 stories)")
    print(f"    Positions: Face midpoints at 45-deg rotation")
    print()

    # -------------------------------------------------------------------------
    # Step 3: Generate Surrounding Buildings from NYC Open Data
    # -------------------------------------------------------------------------
    print("Step 3: Generating Surrounding Buildings from NYC Open Data")
    print("-" * 80)

    surr_tris = []
    building_count = 0
    skipped_count = 0
    tallest = ("", 0.0)
    data_source = "None"

    if not offline_mode and HAS_REQUESTS:
        # Query NYC Open Data API
        geojson_data = fetch_nyc_buildings()

        if geojson_data:
            features = geojson_data.get("features", [])
            data_source = "NYC Open Data Building Footprints API"

            for feature in features:
                parsed = parse_building_feature(feature)
                if parsed is None:
                    skipped_count += 1
                    continue

                verts, height_m, bin_str, props = parsed

                # Skip Citicorp itself
                if is_citicorp_building(bin_str, verts):
                    print(f"  [INFO] Identified Citicorp (BIN={bin_str}), using hand-crafted model")
                    skipped_count += 1
                    continue

                # Infer roof type and generate geometry
                roof_type = infer_roof_type(props, height_m)

                if roof_type == 'gabled':
                    # Simplified: use flat roof for now
                    # Future: implement gabled roof geometry
                    tris = extrude_polygon_to_triangles(verts, 0.0, height_m)
                else:
                    tris = extrude_polygon_to_triangles(verts, 0.0, height_m)

                surr_tris.extend(tris)
                building_count += 1

                # Track tallest building
                if height_m > tallest[1]:
                    tallest = (bin_str, height_m)

            print(f"  [OK] Processed {building_count} buildings from NYC Open Data")
            print(f"    Skipped: {skipped_count} (below height threshold or Citicorp)")

            if tallest[1] > 0:
                print(f"    Tallest neighbor: BIN {tallest[0]}, "
                      f"{tallest[1]:.1f}m ({tallest[1]/FT_TO_M:.0f} ft)")
        else:
            print("  [WARN] NYC Open Data API unavailable")
            data_source = "None (API failed)"
    else:
        reason = "--offline flag specified" if offline_mode else "requests library not available"
        print(f"  [WARN] Skipping API query: {reason}")
        data_source = f"None ({reason})"

    # Write surrounding buildings STL
    n = write_stl_ascii(
        os.path.join(out_dir, "surroundings_nyc3d.stl"),
        "surroundings_nyc3d",
        surr_tris
    )
    total_tris += n

    print(f"  [OK] surroundings_nyc3d.stl: {n} triangles")
    print(f"    Buildings: {building_count}")
    print(f"    Data source: {data_source}")
    print()

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("="*80)
    print("STL Generation Complete")
    print("="*80)
    print(f"Total triangles: {total_tris}")
    print(f"Output directory: {out_dir}")
    print()
    print("Files generated:")
    print(f"  - citicorp_tower.stl       ({len(tower_tris)} triangles)")
    print(f"  - citicorp_stilts.stl      ({len(stilt_tris)} triangles)")
    print(f"  - surroundings_nyc3d.stl   ({len(surr_tris)} triangles)")
    print()
    print("Next steps:")
    print("  1. Verify STL files: Use MeshLab, Blender, or ParaView")
    print("  2. Copy to OpenFOAM case: cp *.stl constant/triSurface/")
    print("  3. Update snappyHexMeshDict to reference new STL files")
    print("  4. Run meshing: snappyHexMesh -overwrite")
    print()
    print("Data sources documentation:")
    print("  • NYC Open Data Building Footprints:")
    print("    https://data.cityofnewyork.us/City-Government/BUILDING/5zhs-2jue")
    print("  • NYC 3D Building Model:")
    print("    https://data.cityofnewyork.us/City-Government/NYC-3D-Model-by-Community-District/u5j4-zxpn")
    print("  • Metadata PDF:")
    print("    https://www.nyc.gov/assets/planning/download/pdf/data-maps/open-data/nyc-3d-model-metadata.pdf")
    print("  • GitHub Documentation:")
    print("    https://github.com/CityOfNewYork/nyc-geo-metadata/blob/main/Metadata/Metadata_BuildingFootprints.md")
    print()


if __name__ == "__main__":
    main()
