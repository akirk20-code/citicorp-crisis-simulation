#!/usr/bin/env python3
"""
generate_stl_lod2.py
====================

Comprehensive LOD2 (Level of Detail 2) STL generator for Citicorp Center CFD.

This script generates three STL files for OpenFOAM CFD meshing:
  - citicorp_tower.stl    : Hand-crafted 45-deg slanted roof, 45-deg rotated tower
  - citicorp_stilts.stl   : Four 10-story stilts at face midpoints
  - surroundings_lod2.stl : Surrounding buildings with LOD2 roof geometry

Data source priority:
  1. Local CityJSON file (--cityjson <file>) - TUM NYC LOD2 or 3D BAG
  2. NYC Open Data Building Footprints API (LOD1 fallback)
  3. Hand-crafted hardcoded geometry (--offline)

===========================================================================
TUM NYC LOD2 CityGML Dataset
===========================================================================

Primary: TUM Chair of Geoinformatics + GeorOcket Enhanced NYC Model
  GitHub:   https://github.com/georocket/new-york-city-model-enhanced
  Format:   CityGML 2.0, with PLUTO semantic attributes
  Size:     ~20 GB full city, ~3 GB Manhattan borough only
  Features: LOD2.2 roof geometry, ~90 semantic attributes per building
            Roof types: flat, gabled, hipped, mansard, shed, complex

  Download steps (Git LFS required for large files):
    git clone https://github.com/georocket/new-york-city-model-enhanced.git
    cd new-york-city-model-enhanced/CityGML/Manhattan/

  Alternative - Convert CityGML → CityJSON on command line:
    pip install cjio
    cjio manhattan.gml upgrade save manhattan.city.json

Secondary: TU Delft 3D BAG
  URL:      https://3d.bk.tudelft.nl/opendata/3dbag/
  Tiles:    https://3d.bk.tudelft.nl/opendata/3dbag/v2/tiles/
  Format:   CityJSON 1.1 (preferred) or CityGML 2.0
  Note:     Coverage is Netherlands-centric; NYC coverage may be limited.

Tertiary: NYC Open Data 3D Building Model (DOITT)
  Dataset:  NYC 3D Model by Community District
  URL:      https://data.cityofnewyork.us/City-Government/NYC-3D-Model-by-Community-District/u5j4-zxpn
  Format:   3DM (Rhinoceros); CityGML via third-party conversion

===========================================================================
CityJSON Schema Structure
===========================================================================

Root object:
  {
    "type": "CityJSON",
    "version": "1.1",
    "transform": {                    # Optional scaling/translation
      "scale": [sx, sy, sz],
      "translate": [tx, ty, tz]
    },
    "vertices": [[x, y, z], ...],     # Shared vertex pool (integer indices)
    "CityObjects": {
      "building_id": {
        "type": "Building",           # or "BuildingPart"
        "attributes": {...},          # BIN, height, year_built, etc.
        "geometry": [
          {
            "type": "MultiSurface",   # or "Solid", "CompositeSurface"
            "lod": "2",
            "boundaries": [           # List of surface rings
              [[v0, v1, v2, ...]],    # Outer ring (vertex indices)
              [[v3, v4, v5, ...]]     # Inner ring (holes) - optional
            ],
            "semantics": {
              "surfaces": [
                {"type": "WallSurface"},
                {"type": "RoofSurface"},
                {"type": "GroundSurface"}
              ],
              "values": [0, 1, 2, ...] # Index into surfaces[] for each boundary
            }
          }
        ]
      }
    }
  }

Surface types (OGC CityGML standard):
  - GroundSurface  : Bottom of building (z = ground level)
  - WallSurface    : Vertical facades
  - RoofSurface    : Sloped or flat roof planes
  - OuterFloorSurface : Balconies, terraces
  - ClosureSurface : Virtual faces closing open models

===========================================================================
Coordinate System Transformations
===========================================================================

WGS84 (EPSG:4326)          → Local Cartesian
  Geographic lat/lon          X = East, Y = North (meters)
  Used by NYC Open Data        Origin = Citicorp (40.7579°N, 73.969°W)

  x = (lon - LON_CITICORP) * 111320 * cos(LAT_CITICORP)
  y = (lat - LAT_CITICORP) * 111320

NAD83 / New York State Plane (EPSG:2263, US Survey Feet) → Local Cartesian
  Used by NYC DOITT 3D model       Convert to meters first:
  x_m = x_ft * 0.3048              Then apply WGS84 transform via pyproj

RD New Dutch Grid (EPSG:28992)  → WGS84 → Local Cartesian
  Used by TU Delft 3D BAG          Requires pyproj:
  from pyproj import Transformer
  t = Transformer.from_crs("EPSG:28992", "EPSG:4326", always_xy=True)
  lon, lat = t.transform(x_rd, y_rd)

For production: use pyproj with full EPSG lookups rather than the
approximate planar formulas in this script.

===========================================================================
Dependencies
===========================================================================

Required (standard library only - always available):
  json, struct, math, sys, os, re, pathlib

Optional (install for enhanced functionality):
  requests   - HTTP API access     : pip install requests
  pyproj     - Coordinate systems  : pip install pyproj
  shapely    - Polygon operations  : pip install shapely
  trimesh    - Mesh validation     : pip install trimesh
  cjio       - CityJSON tools      : pip install cjio

===========================================================================
Usage
===========================================================================

  # Online mode: download from NYC Open Data
  python generate_stl_lod2.py

  # With local CityJSON file (LOD2 mode)
  python generate_stl_lod2.py --cityjson citicorp_area.city.json

  # Offline mode: Citicorp only (no API)
  python generate_stl_lod2.py --offline

  # Validate output STLs (requires trimesh)
  python generate_stl_lod2.py --validate

  # Output to custom directory
  python generate_stl_lod2.py --outdir ./my_stls

Author: Generated for OR 750 Reliability Analysis, GMU PhD Program
Date: 2026-02-16
Version: 2.0
"""

import os
import sys
import math
import struct
import json
import re
import time
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any, Iterator

# ============================================================================
# Optional dependencies with graceful fallbacks
# ============================================================================

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

try:
    from pyproj import Transformer
    HAS_PYPROJ = True
    _TRANSFORMER_RDNEW_TO_WGS84 = None  # Lazy init

    def _get_rdnew_transformer():
        """Lazily initialize RD New → WGS84 transformer."""
        global _TRANSFORMER_RDNEW_TO_WGS84
        if _TRANSFORMER_RDNEW_TO_WGS84 is None:
            _TRANSFORMER_RDNEW_TO_WGS84 = Transformer.from_crs(
                "EPSG:28992", "EPSG:4326", always_xy=True
            )
        return _TRANSFORMER_RDNEW_TO_WGS84

except ImportError:
    HAS_PYPROJ = False
    _TRANSFORMER_RDNEW_TO_WGS84 = None

    def _get_rdnew_transformer():
        return None

try:
    from shapely.geometry import Polygon as ShapelyPolygon
    from shapely.ops import unary_union
    HAS_SHAPELY = True
except ImportError:
    HAS_SHAPELY = False

try:
    import trimesh
    HAS_TRIMESH = True
except ImportError:
    HAS_TRIMESH = False


# ============================================================================
# CITICORP CENTER SPECIFICATIONS (Morgenstern 1995)
# All dimensions in meters. Origin at Citicorp Center ground level.
# ============================================================================

TOWER_W = 47.85       # 157 ft  — square plan side length
TOWER_H_BASE = 244.15 # 801 ft  — top of vertical facade (start of slanted roof)
TOWER_H_TOP = 278.9   # 915 ft  — peak of slanted roof
STILT_H = 34.75       # 114 ft  — stilt top (10 stories)
STILT_W = 7.32        # 24 ft   — stilt cross-section square
HALF_T = TOWER_W / 2  # 23.925 m
HALF_S = STILT_W / 2  # 3.66 m

# Tower orientation: 45° clockwise from North (aligned with midtown Manhattan grid)
TOWER_ROTATION_DEG = 45.0
TOWER_ROTATION_RAD = math.radians(TOWER_ROTATION_DEG)

# Stilt positions at midpoints of unrotated square faces
# These become the structural column positions after 45° rotation
STILTS_UNROTATED: List[Tuple[float, float]] = [
    (0,       -HALF_T),   # South face midpoint
    (HALF_T,   0),        # East face midpoint
    (0,        HALF_T),   # North face midpoint
    (-HALF_T,  0),        # West face midpoint
]


# ============================================================================
# NYC COORDINATE REFERENCE
# ============================================================================

# Citicorp Center location (WGS84)
CITICORP_LAT: float = 40.7579   # degrees North
CITICORP_LON: float = -73.9690  # degrees East (negative = West)

# Approximate planar projection constants (valid within ~50km of NYC)
M_PER_DEG_LAT: float = 111320.0
M_PER_DEG_LON: float = M_PER_DEG_LAT * math.cos(math.radians(CITICORP_LAT))
# = 111320.0 * cos(40.7579°) = 84304.7 m/deg at this latitude

# NYC Open Data Building Footprints (Socrata API)
NYC_BUILDING_API: str = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"

# CFD domain bounds relative to Citicorp center (meters)
DOMAIN_BOUNDS: Dict[str, float] = {
    "xmin": -200, "xmax": 520,
    "ymin": -360, "ymax": 360,
}

# Unit conversions
FT_TO_M: float = 0.3048
IN_TO_M: float = 0.0254

# Geometry filtering thresholds
MIN_BLDG_HEIGHT_M: float = 5.0   # Skip buildings below this height
MIN_BLDG_AREA_M2: float = 10.0   # Skip polygons below this area

# Citicorp identification
CITICORP_BIN_SET: set = {"1035879", "1087931"}
CITICORP_PROXIMITY_M: float = 40.0  # Centroid within this radius → Citicorp


# ============================================================================
# COORDINATE TRANSFORMS
# ============================================================================

def wgs84_to_local(lon: float, lat: float) -> Tuple[float, float]:
    """
    Convert WGS84 geographic coordinates to local Cartesian meters.

    Uses equirectangular (plate carree) approximation valid for small domains
    (< 100 km radius). For production, use pyproj with EPSG:4326 → EPSG:32618.

    Args:
        lon: Longitude in decimal degrees (negative = West)
        lat: Latitude in decimal degrees (positive = North)

    Returns:
        (x, y): Local east/north offsets in meters. Origin at Citicorp Center.
    """
    x = (lon - CITICORP_LON) * M_PER_DEG_LON
    y = (lat - CITICORP_LAT) * M_PER_DEG_LAT
    return x, y


def local_to_wgs84(x: float, y: float) -> Tuple[float, float]:
    """
    Inverse of wgs84_to_local. Local meters → WGS84.

    Args:
        x, y: Local Cartesian coordinates (meters)

    Returns:
        (lon, lat): WGS84 coordinates
    """
    lon = x / M_PER_DEG_LON + CITICORP_LON
    lat = y / M_PER_DEG_LAT + CITICORP_LAT
    return lon, lat


def rdnew_to_local(x_rd: float, y_rd: float) -> Optional[Tuple[float, float]]:
    """
    Convert RD New (Dutch National Grid, EPSG:28992) to local Cartesian.

    The TU Delft 3D BAG dataset uses RD New coordinates for the Netherlands.
    Requires pyproj for accurate transform.

    Args:
        x_rd: RD New easting (meters)
        y_rd: RD New northing (meters)

    Returns:
        (x, y) local Cartesian, or None if pyproj unavailable
    """
    if not HAS_PYPROJ:
        # Approximate fallback not possible without proper CRS definition
        return None

    transformer = _get_rdnew_transformer()
    if transformer is None:
        return None

    lon, lat = transformer.transform(x_rd, y_rd)
    return wgs84_to_local(lon, lat)


def nys_plane_to_local(x_ft: float, y_ft: float) -> Tuple[float, float]:
    """
    Convert NYC state plane coordinates (NAD83/NY-LI, EPSG:2263, US Survey Feet)
    to local Cartesian meters.

    Used by NYC DOITT datasets (groundelev, height_roof in EPSG:2263).
    Provides approximate planar conversion; for accuracy use pyproj.

    Args:
        x_ft: NY State Plane easting (US Survey Feet)
        y_ft: NY State Plane northing (US Survey Feet)

    Returns:
        (x, y) approximate local Cartesian meters
    """
    # NY State Plane Long Island zone origin approx (WGS84):
    # EPSG:2263 false easting = 300000 ft, central meridian ~73.5°W
    # This is an approximation — for accuracy use pyproj
    x_m = x_ft * FT_TO_M
    y_m = y_ft * FT_TO_M

    # Approximate origin in EPSG:2263 feet (Citicorp approximate)
    # Citicorp is at approximately (988700, 214600) in EPSG:2263 feet
    CITICORP_NYSPLANE_X = 988700.0  # approx
    CITICORP_NYSPLANE_Y = 214600.0  # approx

    local_x = (x_ft - CITICORP_NYSPLANE_X) * FT_TO_M
    local_y = (y_ft - CITICORP_NYSPLANE_Y) * FT_TO_M
    return local_x, local_y


def rotate_point_2d(x: float, y: float, angle_rad: float) -> Tuple[float, float]:
    """
    Rotate a 2D point counterclockwise about the origin.

    Args:
        x, y: Point coordinates
        angle_rad: Rotation angle (radians, CCW positive)

    Returns:
        (x_rot, y_rot): Rotated coordinates
    """
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    return (x * cos_a - y * sin_a, x * sin_a + y * cos_a)


# ============================================================================
# POLYGON UTILITIES
# ============================================================================

def polygon_signed_area_2d(verts: List[Tuple[float, float]]) -> float:
    """
    Signed area of a 2D polygon (shoelace formula).

    Returns:
        Area (positive = CCW, negative = CW)
    """
    n = len(verts)
    if n < 3:
        return 0.0
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += verts[i][0] * verts[j][1]
        area -= verts[j][0] * verts[i][1]
    return area * 0.5


def polygon_centroid_2d(verts: List[Tuple[float, float]]) -> Tuple[float, float]:
    """
    Centroid of a 2D polygon.

    Args:
        verts: Polygon vertices

    Returns:
        (cx, cy) centroid coordinates
    """
    if not verts:
        return (0.0, 0.0)
    n = len(verts)
    area = polygon_signed_area_2d(verts)
    if abs(area) < 1e-15:
        # Degenerate — use mean
        cx = sum(v[0] for v in verts) / n
        cy = sum(v[1] for v in verts) / n
        return cx, cy

    cx = 0.0
    cy = 0.0
    for i in range(n):
        j = (i + 1) % n
        factor = verts[i][0] * verts[j][1] - verts[j][0] * verts[i][1]
        cx += (verts[i][0] + verts[j][0]) * factor
        cy += (verts[i][1] + verts[j][1]) * factor
    cx /= (6.0 * area)
    cy /= (6.0 * area)
    return cx, cy


def polygon_longest_edge(verts: List[Tuple[float, float]]) -> Tuple[int, int, float]:
    """
    Find the longest edge of a polygon.

    Returns:
        (i, j, length): Start index, end index, and edge length
    """
    best = (0, 1, 0.0)
    n = len(verts)
    for i in range(n):
        j = (i + 1) % n
        dx = verts[j][0] - verts[i][0]
        dy = verts[j][1] - verts[i][1]
        length = math.sqrt(dx * dx + dy * dy)
        if length > best[2]:
            best = (i, j, length)
    return best


def ensure_ccw(verts: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """
    Ensure polygon vertices are in CCW order (positive signed area).

    Args:
        verts: Polygon vertices

    Returns:
        Vertices in CCW order (reversed if necessary)
    """
    if polygon_signed_area_2d(verts) < 0:
        return list(reversed(verts))
    return list(verts)


def remove_duplicate_close_vertex(
        verts: List[Tuple[float, float]],
        tol: float = 0.01) -> List[Tuple[float, float]]:
    """
    Remove closing duplicate vertex (GeoJSON polygons repeat first vertex).

    Args:
        verts: Polygon vertex list
        tol: Distance tolerance in meters

    Returns:
        Cleaned vertex list
    """
    if len(verts) < 2:
        return verts
    dx = verts[0][0] - verts[-1][0]
    dy = verts[0][1] - verts[-1][1]
    if math.sqrt(dx * dx + dy * dy) < tol:
        return verts[:-1]
    return verts


def point_in_triangle_2d(
        p: Tuple[float, float],
        a: Tuple[float, float],
        b: Tuple[float, float],
        c: Tuple[float, float]) -> bool:
    """
    Test if point p is inside (or on boundary of) triangle (a, b, c).

    Uses barycentric sign method.
    """
    def sign3(p1, p2, p3):
        return ((p1[0] - p3[0]) * (p2[1] - p3[1])
                - (p2[0] - p3[0]) * (p1[1] - p3[1]))

    d1 = sign3(p, a, b)
    d2 = sign3(p, b, c)
    d3 = sign3(p, c, a)

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (has_neg and has_pos)


def ear_clip_triangulate(verts: List[Tuple[float, float]]) -> List[Tuple[int, int, int]]:
    """
    Ear-clipping triangulation of a simple (non-self-intersecting) 2D polygon.

    Assumes CCW winding. Falls back to fan triangulation for degenerate cases.

    Complexity: O(n^2) — adequate for building footprints (n < 100 vertices).

    Args:
        verts: List of (x, y) vertices in CCW order

    Returns:
        List of (i, j, k) index triples referencing the original vertex list
    """
    n = len(verts)
    if n < 3:
        return []
    if n == 3:
        return [(0, 1, 2)]

    active = list(range(n))  # Indices of remaining (non-clipped) vertices
    tris: List[Tuple[int, int, int]] = []
    max_iters = n * n  # Prevent infinite loop on degenerate inputs

    while len(active) > 3 and max_iters > 0:
        max_iters -= 1
        m = len(active)
        ear_found = False

        for i in range(m):
            prev_idx = active[(i - 1) % m]
            curr_idx = active[i]
            next_idx = active[(i + 1) % m]

            ax, ay = verts[prev_idx]
            bx, by = verts[curr_idx]
            cx, cy = verts[next_idx]

            # Cross product z-component: positive → convex vertex (for CCW polygon)
            cross_z = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            if cross_z <= 1e-12:
                continue  # Concave vertex — cannot be an ear

            # Verify no other active vertex lies inside this triangle
            is_ear = True
            for j in range(m):
                if j in ((i - 1) % m, i, (i + 1) % m):
                    continue
                if point_in_triangle_2d(
                        verts[active[j]],
                        verts[prev_idx], verts[curr_idx], verts[next_idx]):
                    is_ear = False
                    break

            if is_ear:
                tris.append((prev_idx, curr_idx, next_idx))
                active.pop(i)
                ear_found = True
                break

        if not ear_found:
            # Fallback: fan triangulation from active[0]
            for i in range(1, len(active) - 1):
                tris.append((active[0], active[i], active[i + 1]))
            return tris

    if len(active) == 3:
        tris.append((active[0], active[1], active[2]))

    return tris


# ============================================================================
# 3D VECTOR MATH
# ============================================================================

def vec3_sub(a: Tuple[float, float, float],
             b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vec3_cross(u: Tuple[float, float, float],
               v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def vec3_normalize(v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    length = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if length < 1e-15:
        return (0.0, 0.0, 1.0)
    return (v[0] / length, v[1] / length, v[2] / length)


def face_normal(v0: Tuple[float, float, float],
                v1: Tuple[float, float, float],
                v2: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """
    Compute outward unit normal for triangle (v0, v1, v2) by right-hand rule.
    """
    return vec3_normalize(vec3_cross(vec3_sub(v1, v0), vec3_sub(v2, v0)))


# ============================================================================
# STL GEOMETRY BUILDERS
# ============================================================================

Triangle3D = Tuple[Tuple[Tuple[float, float, float],
                          Tuple[float, float, float],
                          Tuple[float, float, float]],
                    Tuple[float, float, float]]


def box_triangles(xmin: float, ymin: float, zmin: float,
                  xmax: float, ymax: float, zmax: float) -> List[Triangle3D]:
    """
    12 outward-normal triangles for an axis-aligned box.

    Args:
        xmin/ymin/zmin: Minimum corner
        xmax/ymax/zmax: Maximum corner

    Returns:
        List of ((v0, v1, v2), normal) triangle tuples
    """
    v = [
        (xmin, ymin, zmin),  # 0: ---
        (xmax, ymin, zmin),  # 1: +--
        (xmax, ymax, zmin),  # 2: ++-
        (xmin, ymax, zmin),  # 3: -+-
        (xmin, ymin, zmax),  # 4: --+
        (xmax, ymin, zmax),  # 5: +-+
        (xmax, ymax, zmax),  # 6: +++
        (xmin, ymax, zmax),  # 7: -++
    ]
    faces = [
        ((0, 2, 1), ( 0,  0, -1)),  # bottom (z-)
        ((0, 3, 2), ( 0,  0, -1)),
        ((4, 5, 6), ( 0,  0,  1)),  # top (z+)
        ((4, 6, 7), ( 0,  0,  1)),
        ((0, 1, 5), ( 0, -1,  0)),  # south (y-)
        ((0, 5, 4), ( 0, -1,  0)),
        ((2, 3, 7), ( 0,  1,  0)),  # north (y+)
        ((2, 7, 6), ( 0,  1,  0)),
        ((0, 4, 7), (-1,  0,  0)),  # west (x-)
        ((0, 7, 3), (-1,  0,  0)),
        ((1, 2, 6), ( 1,  0,  0)),  # east (x+)
        ((1, 6, 5), ( 1,  0,  0)),
    ]
    return [(tuple(v[i] for i in idx), n) for idx, n in faces]


def extrude_flat_roof(
        footprint: List[Tuple[float, float]],
        z_ground: float,
        z_top: float) -> List[Triangle3D]:
    """
    Extrude a 2D polygon footprint into a watertight 3D prism with flat top.

    Generates:
      - Top face (normal +Z)
      - Bottom face (normal -Z)
      - Side walls (outward normals from CCW winding)

    Args:
        footprint: CCW list of (x, y) vertices
        z_ground: Ground elevation
        z_top: Roof elevation

    Returns:
        List of ((v0, v1, v2), normal) triangles
    """
    n = len(footprint)
    if n < 3 or z_top <= z_ground:
        return []

    tris_2d = ear_clip_triangulate(footprint)
    result: List[Triangle3D] = []

    # Top cap
    for i0, i1, i2 in tris_2d:
        v0 = (footprint[i0][0], footprint[i0][1], z_top)
        v1 = (footprint[i1][0], footprint[i1][1], z_top)
        v2 = (footprint[i2][0], footprint[i2][1], z_top)
        result.append(((v0, v1, v2), (0.0, 0.0, 1.0)))

    # Bottom cap (reversed winding for downward normal)
    for i0, i1, i2 in tris_2d:
        v0 = (footprint[i0][0], footprint[i0][1], z_ground)
        v1 = (footprint[i1][0], footprint[i1][1], z_ground)
        v2 = (footprint[i2][0], footprint[i2][1], z_ground)
        result.append(((v0, v2, v1), (0.0, 0.0, -1.0)))

    # Vertical walls (one quad per footprint edge = 2 triangles)
    for i in range(n):
        j = (i + 1) % n
        bl = (footprint[i][0], footprint[i][1], z_ground)
        br = (footprint[j][0], footprint[j][1], z_ground)
        tl = (footprint[i][0], footprint[i][1], z_top)
        tr = (footprint[j][0], footprint[j][1], z_top)

        # Outward normal for CCW polygon edge: rotate edge vector 90° right
        dx = footprint[j][0] - footprint[i][0]
        dy = footprint[j][1] - footprint[i][1]
        wall_n = vec3_normalize((dy, -dx, 0.0))

        result.append(((bl, br, tr), wall_n))
        result.append(((bl, tr, tl), wall_n))

    return result


def extrude_gabled_roof(
        footprint: List[Tuple[float, float]],
        z_ground: float,
        z_eave: float,
        z_ridge: float) -> List[Triangle3D]:
    """
    Extrude a polygon with a gabled (peaked) roof.

    Ridge runs parallel to the longest edge of the footprint.
    Each side of the gable is triangulated independently.

    Args:
        footprint: CCW list of (x, y) vertices
        z_ground: Ground elevation
        z_eave: Wall top / eave height
        z_ridge: Ridge peak height

    Returns:
        List of ((v0, v1, v2), normal) triangles
    """
    n = len(footprint)
    if n < 3 or z_eave <= z_ground:
        return extrude_flat_roof(footprint, z_ground, z_eave)

    result: List[Triangle3D] = []

    # Build walls (ground to eave)
    for i in range(n):
        j = (i + 1) % n
        bl = (footprint[i][0], footprint[i][1], z_ground)
        br = (footprint[j][0], footprint[j][1], z_ground)
        tl = (footprint[i][0], footprint[i][1], z_eave)
        tr = (footprint[j][0], footprint[j][1], z_eave)

        dx = footprint[j][0] - footprint[i][0]
        dy = footprint[j][1] - footprint[i][1]
        wall_n = vec3_normalize((dy, -dx, 0.0))

        result.append(((bl, br, tr), wall_n))
        result.append(((bl, tr, tl), wall_n))

    # Bottom cap
    tris_2d = ear_clip_triangulate(footprint)
    for i0, i1, i2 in tris_2d:
        v0 = (footprint[i0][0], footprint[i0][1], z_ground)
        v1 = (footprint[i1][0], footprint[i1][1], z_ground)
        v2 = (footprint[i2][0], footprint[i2][1], z_ground)
        result.append(((v0, v2, v1), (0.0, 0.0, -1.0)))

    # Gabled roof: find longest edge as ridge direction
    i_ridge, j_ridge, _ = polygon_longest_edge(footprint)

    # Ridge vertices (elevated above eave)
    rp0 = (footprint[i_ridge][0], footprint[i_ridge][1], z_ridge)
    rp1 = (footprint[j_ridge][0], footprint[j_ridge][1], z_ridge)

    # Project all footprint vertices onto ridge axis to decide which side
    ridge_dx = footprint[j_ridge][0] - footprint[i_ridge][0]
    ridge_dy = footprint[j_ridge][1] - footprint[i_ridge][1]
    ridge_len = math.sqrt(ridge_dx * ridge_dx + ridge_dy * ridge_dy)

    if ridge_len < 1e-10:
        # Degenerate — fallback to flat
        return extrude_flat_roof(footprint, z_ground, z_eave)

    # Unit ridge direction and perpendicular
    ux = ridge_dx / ridge_len
    uy = ridge_dy / ridge_len
    px = -uy   # Perpendicular to ridge (90° CCW)
    py = ux

    # Split vertices into two sides based on signed perpendicular distance
    origin_x = footprint[i_ridge][0]
    origin_y = footprint[i_ridge][1]

    side_a: List[int] = []  # vertices with positive perp distance
    side_b: List[int] = []  # vertices with negative perp distance

    for k in range(n):
        dx = footprint[k][0] - origin_x
        dy = footprint[k][1] - origin_y
        perp_dist = dx * px + dy * py
        if perp_dist >= 0:
            side_a.append(k)
        else:
            side_b.append(k)

    # Create roof panels: each vertex fans to the ridge
    for side_indices in (side_a, side_b):
        if not side_indices:
            continue
        # Include ridge endpoints in eave-ring for this side
        eave_verts = [(footprint[k][0], footprint[k][1], z_eave)
                      for k in side_indices]
        # Connect to ridge via triangles (simplified fan)
        for k in range(len(eave_verts) - 1):
            v0 = eave_verts[k]
            v1 = eave_verts[k + 1]
            # Assign nearest ridge endpoint
            mid_x = (v0[0] + v1[0]) / 2
            mid_y = (v0[1] + v1[1]) / 2
            # Project mid-eave onto ridge to find interpolated ridge height
            t = ((mid_x - origin_x) * ux + (mid_y - origin_y) * uy) / ridge_len
            t = max(0.0, min(1.0, t))
            ridge_pt = (
                footprint[i_ridge][0] + t * ridge_dx,
                footprint[i_ridge][1] + t * ridge_dy,
                z_ridge
            )
            n_face = face_normal(v0, v1, ridge_pt)
            result.append(((v0, v1, ridge_pt), n_face))

    return result


def extrude_hipped_roof(
        footprint: List[Tuple[float, float]],
        z_ground: float,
        z_eave: float,
        z_peak: float) -> List[Triangle3D]:
    """
    Extrude a polygon with a hipped (pyramid-like) roof.

    The roof converges to a single apex above the footprint centroid.

    Args:
        footprint: CCW list of (x, y) vertices
        z_ground: Ground elevation
        z_eave: Eave height
        z_peak: Apex height

    Returns:
        List of ((v0, v1, v2), normal) triangles
    """
    n = len(footprint)
    if n < 3 or z_eave <= z_ground:
        return extrude_flat_roof(footprint, z_ground, z_eave)

    result: List[Triangle3D] = []

    # Walls (ground to eave)
    for i in range(n):
        j = (i + 1) % n
        bl = (footprint[i][0], footprint[i][1], z_ground)
        br = (footprint[j][0], footprint[j][1], z_ground)
        tl = (footprint[i][0], footprint[i][1], z_eave)
        tr = (footprint[j][0], footprint[j][1], z_eave)

        dx = footprint[j][0] - footprint[i][0]
        dy = footprint[j][1] - footprint[i][1]
        wall_n = vec3_normalize((dy, -dx, 0.0))

        result.append(((bl, br, tr), wall_n))
        result.append(((bl, tr, tl), wall_n))

    # Bottom cap
    tris_2d = ear_clip_triangulate(footprint)
    for i0, i1, i2 in tris_2d:
        v0 = (footprint[i0][0], footprint[i0][1], z_ground)
        v1 = (footprint[i1][0], footprint[i1][1], z_ground)
        v2 = (footprint[i2][0], footprint[i2][1], z_ground)
        result.append(((v0, v2, v1), (0.0, 0.0, -1.0)))

    # Hipped roof: apex above centroid
    cx, cy = polygon_centroid_2d(footprint)
    apex = (cx, cy, z_peak)

    for i in range(n):
        j = (i + 1) % n
        eave_v0 = (footprint[i][0], footprint[i][1], z_eave)
        eave_v1 = (footprint[j][0], footprint[j][1], z_eave)

        n_face = face_normal(eave_v0, eave_v1, apex)
        result.append(((eave_v0, eave_v1, apex), n_face))

    return result


def extrude_mansard_roof(
        footprint: List[Tuple[float, float]],
        z_ground: float,
        z_eave: float,
        z_setback: float,
        setback_dist: float) -> List[Triangle3D]:
    """
    Extrude a polygon with a mansard roof (inward-sloping facades + flat top).

    Mansard roofs have:
      - Vertical or sloped lower section (eave to setback)
      - Flat top with inset footprint (offset inward by setback_dist)

    Args:
        footprint: CCW list of (x, y) vertices
        z_ground: Ground elevation
        z_eave: Eave height (base of roof)
        z_setback: Top of sloped section (flat top height)
        setback_dist: Inward offset for flat top footprint (meters)

    Returns:
        List of ((v0, v1, v2), normal) triangles
    """
    n = len(footprint)
    if n < 3:
        return []

    result: List[Triangle3D] = []

    # Compute inset footprint for flat top
    cx, cy = polygon_centroid_2d(footprint)
    inner_footprint = []
    for vx, vy in footprint:
        # Move vertex toward centroid by setback_dist
        dx = cx - vx
        dy = cy - vy
        dist = math.sqrt(dx * dx + dy * dy)
        if dist > setback_dist:
            factor = setback_dist / dist
            inner_footprint.append((vx + dx * factor, vy + dy * factor))
        else:
            inner_footprint.append((cx, cy))

    # Walls (ground to eave)
    for i in range(n):
        j = (i + 1) % n
        bl = (footprint[i][0], footprint[i][1], z_ground)
        br = (footprint[j][0], footprint[j][1], z_ground)
        tl = (footprint[i][0], footprint[i][1], z_eave)
        tr = (footprint[j][0], footprint[j][1], z_eave)

        dx = footprint[j][0] - footprint[i][0]
        dy = footprint[j][1] - footprint[i][1]
        wall_n = vec3_normalize((dy, -dx, 0.0))

        result.append(((bl, br, tr), wall_n))
        result.append(((bl, tr, tl), wall_n))

    # Bottom cap
    tris_2d = ear_clip_triangulate(footprint)
    for i0, i1, i2 in tris_2d:
        v0 = (footprint[i0][0], footprint[i0][1], z_ground)
        v1 = (footprint[i1][0], footprint[i1][1], z_ground)
        v2 = (footprint[i2][0], footprint[i2][1], z_ground)
        result.append(((v0, v2, v1), (0.0, 0.0, -1.0)))

    # Mansard slopes (eave to setback, connecting outer to inner)
    for i in range(n):
        j = (i + 1) % n
        bl = (footprint[i][0], footprint[i][1], z_eave)
        br = (footprint[j][0], footprint[j][1], z_eave)
        tl = (inner_footprint[i][0], inner_footprint[i][1], z_setback)
        tr = (inner_footprint[j][0], inner_footprint[j][1], z_setback)

        n_face = face_normal(bl, br, tr)
        result.append(((bl, br, tr), n_face))
        result.append(((bl, tr, tl), n_face))

    # Flat top (inner footprint at z_setback)
    if len(inner_footprint) >= 3:
        inner_ccw = ensure_ccw(inner_footprint)
        inner_tris = ear_clip_triangulate(inner_ccw)
        for i0, i1, i2 in inner_tris:
            v0 = (inner_ccw[i0][0], inner_ccw[i0][1], z_setback)
            v1 = (inner_ccw[i1][0], inner_ccw[i1][1], z_setback)
            v2 = (inner_ccw[i2][0], inner_ccw[i2][1], z_setback)
            result.append(((v0, v1, v2), (0.0, 0.0, 1.0)))

    return result


# ============================================================================
# CITYJSON PARSER
# ============================================================================

class CityJSONParser:
    """
    Parser for CityJSON 1.0 and 1.1 format files.

    CityJSON is the JSON encoding of the OGC CityGML standard.
    Reference: https://www.cityjson.org/specs/1.1.3/

    Key schema features:
      - Shared vertex pool (indexed geometry)
      - Optional transform (scale + translate for integer coordinates)
      - Semantic surface types (WallSurface, RoofSurface, GroundSurface)
      - Multiple geometry types (MultiSurface, Solid, CompositeSolid)
      - Hierarchical city objects (Building → BuildingPart)

    Example usage:
      parser = CityJSONParser(filepath)
      parser.load()
      buildings = parser.iter_buildings(domain_bbox=(-200, 520, -360, 360))
      for building in buildings:
          print(building['id'], building['roof_type'], building['height_m'])
    """

    BUILDING_TYPES = {"Building", "BuildingPart", "BuildingInstallation"}

    WALL_SURFACE_TYPES = {
        "WallSurface", "OuterWallSurface", "ClosureSurface"
    }
    ROOF_SURFACE_TYPES = {
        "RoofSurface"
    }
    GROUND_SURFACE_TYPES = {
        "GroundSurface"
    }

    def __init__(self, filepath: str):
        """
        Initialize parser.

        Args:
            filepath: Path to .city.json file
        """
        self.filepath = filepath
        self.data: Optional[Dict] = None
        self.vertices: List[Tuple[float, float, float]] = []
        self.scale: Tuple[float, float, float] = (1.0, 1.0, 1.0)
        self.translate: Tuple[float, float, float] = (0.0, 0.0, 0.0)
        self.version: str = "unknown"
        self.epsg: Optional[int] = None
        self._loaded: bool = False

    def load(self) -> bool:
        """
        Load and parse CityJSON file.

        Returns:
            True on success, False on error
        """
        print(f"  Loading CityJSON: {self.filepath}")
        try:
            with open(self.filepath, 'r', encoding='utf-8') as f:
                self.data = json.load(f)
        except FileNotFoundError:
            print(f"  Error: File not found: {self.filepath}")
            return False
        except json.JSONDecodeError as e:
            print(f"  Error: JSON parse error at {e.lineno}:{e.colno}: {e.msg}")
            return False
        except MemoryError:
            print("  Error: File too large for available memory")
            return False
        except Exception as e:
            print(f"  Error: {e}")
            return False

        # Validate CityJSON
        if self.data.get("type") != "CityJSON":
            print(f"  Error: Not a CityJSON file (type={self.data.get('type')})")
            return False

        self.version = self.data.get("version", "unknown")

        # Extract coordinate transform (optional)
        xform = self.data.get("transform", {})
        s = xform.get("scale", [1.0, 1.0, 1.0])
        t = xform.get("translate", [0.0, 0.0, 0.0])
        self.scale = (float(s[0]), float(s[1]), float(s[2]))
        self.translate = (float(t[0]), float(t[1]), float(t[2]))

        # Extract CRS/EPSG if available
        meta = self.data.get("metadata", {})
        crs = meta.get("referenceSystem", "")
        m = re.search(r'(\d{4,5})$', crs)
        if m:
            self.epsg = int(m.group(1))

        # Pre-compute scaled vertices
        raw_verts = self.data.get("vertices", [])
        self.vertices = [
            (
                v[0] * self.scale[0] + self.translate[0],
                v[1] * self.scale[1] + self.translate[1],
                v[2] * self.scale[2] + self.translate[2],
            )
            for v in raw_verts
        ]

        n_objs = len(self.data.get("CityObjects", {}))
        print(f"  [OK] CityJSON v{self.version}, EPSG:{self.epsg}")
        print(f"    {len(self.vertices)} vertices, {n_objs} city objects")
        self._loaded = True
        return True

    def _get_vertex(self, idx: int) -> Optional[Tuple[float, float, float]]:
        """Get vertex by index (with bounds check)."""
        if 0 <= idx < len(self.vertices):
            return self.vertices[idx]
        return None

    def _decode_ring_3d(self, ring: List[int]) -> List[Tuple[float, float, float]]:
        """
        Decode a vertex-index ring into 3D world coordinates.

        Args:
            ring: List of integer vertex indices

        Returns:
            List of (x, y, z) world coordinates
        """
        verts = []
        for idx in ring:
            v = self._get_vertex(idx)
            if v is not None:
                verts.append(v)
        return verts

    def _ring_to_local_xy(self,
                           ring_3d: List[Tuple[float, float, float]],
                           coord_type: str = "wgs84"
                           ) -> Tuple[List[Tuple[float, float]], float]:
        """
        Project a 3D surface ring to 2D local Cartesian coordinates.

        Args:
            ring_3d: List of (x, y, z) world coordinates
            coord_type: One of 'wgs84', 'rdnew', 'nysplane', 'local_m'

        Returns:
            (verts_2d, z_mean): 2D footprint and mean Z elevation
        """
        verts_2d = []
        z_vals = []

        for wx, wy, wz in ring_3d:
            if coord_type == "wgs84":
                lx, ly = wgs84_to_local(wx, wy)
            elif coord_type == "rdnew" and HAS_PYPROJ:
                result = rdnew_to_local(wx, wy)
                if result is None:
                    continue
                lx, ly = result
            elif coord_type == "nysplane":
                lx, ly = nys_plane_to_local(wx, wy)
            else:
                # Assume already in local Cartesian meters
                lx, ly = wx, wy

            verts_2d.append((lx, ly))
            z_vals.append(wz)

        z_mean = sum(z_vals) / len(z_vals) if z_vals else 0.0
        return verts_2d, z_mean

    def _get_coord_type(self) -> str:
        """
        Determine coordinate type from EPSG code.

        Returns:
            Coordinate type string for _ring_to_local_xy
        """
        if self.epsg is None:
            return "wgs84"  # Default assumption
        if self.epsg == 4326:
            return "wgs84"
        if self.epsg == 28992:
            return "rdnew"
        if self.epsg == 2263:
            return "nysplane"
        if self.epsg in (32618, 32619):
            return "local_m"  # UTM already in meters
        return "wgs84"  # Conservative fallback

    def _extract_geometry_surfaces(self,
                                    geom: Dict,
                                    coord_type: str
                                    ) -> Dict[str, List[List[Tuple[float, float, float]]]]:
        """
        Extract named surface groups from a CityJSON geometry object.

        Handles MultiSurface and Solid geometry types.
        Groups surfaces by semantic type (WallSurface, RoofSurface, etc.)

        Args:
            geom: CityJSON geometry dict
            coord_type: Coordinate system type

        Returns:
            Dict mapping surface type → list of 3D rings
        """
        surfaces: Dict[str, List] = {
            "WallSurface": [],
            "RoofSurface": [],
            "GroundSurface": [],
            "Other": [],
        }

        geom_type = geom.get("type", "")
        boundaries = geom.get("boundaries", [])
        semantics = geom.get("semantics", {})
        sem_surfaces = semantics.get("surfaces", [])
        sem_values = semantics.get("values", [])

        def classify_surface(surface_idx: Optional[int]) -> str:
            """Return surface type label for a given semantic index."""
            if surface_idx is None:
                return "Other"
            if 0 <= surface_idx < len(sem_surfaces):
                stype = sem_surfaces[surface_idx].get("type", "Other")
                if stype in self.WALL_SURFACE_TYPES:
                    return "WallSurface"
                if stype in self.ROOF_SURFACE_TYPES:
                    return "RoofSurface"
                if stype in self.GROUND_SURFACE_TYPES:
                    return "GroundSurface"
            return "Other"

        if geom_type == "MultiSurface":
            # boundaries = [surface, ...]
            # Each surface = [outer_ring, inner_ring, ...]
            # outer_ring = [v0_idx, v1_idx, ...]
            for surf_i, surface in enumerate(boundaries):
                sem_idx = sem_values[surf_i] if surf_i < len(sem_values) else None
                stype = classify_surface(sem_idx)

                if surface and surface[0]:
                    ring_3d = self._decode_ring_3d(surface[0])  # Outer ring only
                    surfaces[stype].append(ring_3d)

        elif geom_type == "Solid":
            # boundaries = [shell, ...]
            # Each shell = [surface, ...]
            for shell in boundaries:
                for surf_i, surface in enumerate(shell):
                    sem_idx = sem_values[surf_i] if surf_i < len(sem_values) else None
                    stype = classify_surface(sem_idx)

                    if surface and surface[0]:
                        ring_3d = self._decode_ring_3d(surface[0])
                        surfaces[stype].append(ring_3d)

        elif geom_type == "CompositeSurface":
            for surf_i, surface in enumerate(boundaries):
                sem_idx = sem_values[surf_i] if surf_i < len(sem_values) else None
                stype = classify_surface(sem_idx)
                if surface and surface[0]:
                    ring_3d = self._decode_ring_3d(surface[0])
                    surfaces[stype].append(ring_3d)

        return surfaces

    def _surfaces_to_triangles(self,
                                surfaces: Dict[str, List],
                                coord_type: str) -> List[Triangle3D]:
        """
        Convert extracted 3D surface rings to STL triangles.

        Args:
            surfaces: Dict of surface type → list of 3D rings
            coord_type: Coordinate system type

        Returns:
            List of ((v0, v1, v2), normal) triangles
        """
        result: List[Triangle3D] = []

        for stype, ring_list in surfaces.items():
            for ring_3d in ring_list:
                if len(ring_3d) < 3:
                    continue

                # Project to 2D for triangulation
                verts_2d, z_mean = self._ring_to_local_xy(ring_3d, coord_type)

                if len(verts_2d) < 3:
                    continue

                # Remove duplicate closing vertex
                verts_2d = remove_duplicate_close_vertex(verts_2d)

                if len(verts_2d) < 3:
                    continue

                # Ensure CCW winding
                verts_2d = ensure_ccw(verts_2d)

                # Skip degenerate polygons
                if abs(polygon_signed_area_2d(verts_2d)) < 0.5:
                    continue

                # Triangulate with correct Z
                tris_2d = ear_clip_triangulate(verts_2d)

                # Convert ring to 3D with original Z values
                # Map 2D projected vertices back to 3D
                ring_local = []
                for wx, wy, wz in ring_3d:
                    if coord_type == "wgs84":
                        lx, ly = wgs84_to_local(wx, wy)
                    elif coord_type == "rdnew" and HAS_PYPROJ:
                        r = rdnew_to_local(wx, wy)
                        if r is None:
                            continue
                        lx, ly = r
                    else:
                        lx, ly = wx, wy
                    ring_local.append((lx, ly, wz))

                if len(ring_local) < 3:
                    continue

                # Use original 3D ring for triangulation (preserves Z gradients)
                for i0, i1, i2 in tris_2d:
                    if i0 >= len(ring_local) or i1 >= len(ring_local) or i2 >= len(ring_local):
                        continue
                    v0 = ring_local[i0]
                    v1 = ring_local[i1]
                    v2 = ring_local[i2]
                    n = face_normal(v0, v1, v2)
                    result.append(((v0, v1, v2), n))

        return result

    def iter_buildings(
            self,
            domain_bbox: Optional[Tuple[float, float, float, float]] = None,
            lod_preference: str = "2"
    ) -> Iterator[Dict[str, Any]]:
        """
        Iterate over buildings, yielding parsed building data.

        Args:
            domain_bbox: Optional (xmin, xmax, ymin, ymax) in local meters
            lod_preference: Preferred LOD ('2', '1', etc.)

        Yields:
            Dict with keys:
              'id': building ID string
              'triangles': List of ((v0,v1,v2), normal) STL triangles
              'footprint': 2D CCW footprint vertices (local meters)
              'height_m': Roof height in meters
              'z_ground': Ground elevation in meters
              'roof_type': 'flat', 'gabled', 'hipped', 'mansard', 'unknown'
              'attributes': raw attribute dict from CityJSON
              'lod': LOD used ('1', '2', etc.)
        """
        if not self._loaded:
            return

        coord_type = self._get_coord_type()
        city_objects = self.data.get("CityObjects", {})

        for obj_id, obj_data in city_objects.items():
            obj_type = obj_data.get("type", "")
            if obj_type not in self.BUILDING_TYPES:
                continue

            attributes = obj_data.get("attributes", {})
            geometry = obj_data.get("geometry", [])

            if not geometry:
                continue

            # Select best available geometry (prefer requested LOD)
            geom_selected = None
            for geom in geometry:
                lod = str(geom.get("lod", "0"))
                if lod == lod_preference:
                    geom_selected = geom
                    break

            if geom_selected is None:
                # Use first available geometry
                geom_selected = geometry[0]

            lod_used = str(geom_selected.get("lod", "0"))

            # Extract surfaces
            surfaces = self._extract_geometry_surfaces(geom_selected, coord_type)

            # Build 3D triangles from surfaces
            triangles = self._surfaces_to_triangles(surfaces, coord_type)

            if not triangles:
                continue

            # Extract footprint from ground surface (or estimate from walls)
            footprint_2d: List[Tuple[float, float]] = []
            z_ground = 0.0
            z_roof = 0.0

            ground_surfs = surfaces.get("GroundSurface", [])
            roof_surfs = surfaces.get("RoofSurface", [])

            if ground_surfs:
                ring_3d = ground_surfs[0]
                fp, z_g = self._ring_to_local_xy(ring_3d, coord_type)
                footprint_2d = ensure_ccw(remove_duplicate_close_vertex(fp))
                z_ground = z_g

            if roof_surfs:
                all_z = []
                for ring_3d in roof_surfs:
                    for _, _, wz in ring_3d:
                        all_z.append(wz)
                if all_z:
                    z_roof = max(all_z)

            height_m = z_roof - z_ground if z_roof > z_ground else 0.0

            # Override from attributes if available
            for key in ("measuredHeight", "height", "h_dak_max"):
                if key in attributes:
                    try:
                        height_m = float(attributes[key])
                    except (TypeError, ValueError):
                        pass
                    break

            # Infer roof type from attributes
            roof_type_str = "unknown"
            for key in ("roof_type", "dak_type", "roofType"):
                if key in attributes:
                    raw = str(attributes[key]).lower()
                    if "flat" in raw or "plat" in raw:
                        roof_type_str = "flat"
                    elif "gable" in raw or "zadel" in raw:
                        roof_type_str = "gabled"
                    elif "hip" in raw or "schilddak" in raw:
                        roof_type_str = "hipped"
                    elif "mansard" in raw:
                        roof_type_str = "mansard"
                    else:
                        roof_type_str = "unknown"
                    break

            # Domain filter
            if domain_bbox and footprint_2d:
                xmin_d, xmax_d, ymin_d, ymax_d = domain_bbox
                cx, cy = polygon_centroid_2d(footprint_2d)
                if not (xmin_d <= cx <= xmax_d and ymin_d <= cy <= ymax_d):
                    continue

            yield {
                "id": obj_id,
                "triangles": triangles,
                "footprint": footprint_2d,
                "height_m": height_m,
                "z_ground": z_ground,
                "roof_type": roof_type_str,
                "attributes": attributes,
                "lod": lod_used,
            }


# ============================================================================
# NYC OPEN DATA API (LOD1 Fallback)
# ============================================================================

def _bbox_from_domain(pad_m: float = 50.0) -> Tuple[float, float, float, float]:
    """
    Compute WGS84 bounding box covering the CFD domain plus padding.

    Returns:
        (nw_lat, nw_lon, se_lat, se_lon)
    """
    nw_lat = CITICORP_LAT + (DOMAIN_BOUNDS["ymax"] + pad_m) / M_PER_DEG_LAT
    nw_lon = CITICORP_LON + (DOMAIN_BOUNDS["xmin"] - pad_m) / M_PER_DEG_LON
    se_lat = CITICORP_LAT + (DOMAIN_BOUNDS["ymin"] - pad_m) / M_PER_DEG_LAT
    se_lon = CITICORP_LON + (DOMAIN_BOUNDS["xmax"] + pad_m) / M_PER_DEG_LON
    return nw_lat, nw_lon, se_lat, se_lon


def fetch_nyc_open_data(
        limit: int = 5000,
        timeout: int = 60) -> Optional[Dict]:
    """
    Query NYC Open Data Building Footprints (Socrata GeoJSON endpoint).

    Endpoint: https://data.cityofnewyork.us/resource/5zhs-2jue.geojson
    Dataset: NYC Building Footprints (DOITT)
    Auth: None required

    Rate limits: ~1000 req/hour without app token. Register at:
      https://dev.socrata.com/foundry/data.cityofnewyork.us/5zhs-2jue

    Args:
        limit: Maximum number of features to request
        timeout: Request timeout in seconds

    Returns:
        GeoJSON FeatureCollection dict, or None on error
    """
    if not HAS_REQUESTS:
        print("  Error: 'requests' library required. Install: pip install requests")
        return None

    nw_lat, nw_lon, se_lat, se_lon = _bbox_from_domain()
    params = {
        "$where": f"within_box(the_geom, {nw_lat}, {nw_lon}, {se_lat}, {se_lon})",
        "$limit": limit,
        "$order": "height_roof DESC",
    }

    print(f"  Querying NYC Open Data API (limit={limit})...")
    print(f"  Bbox: NW({nw_lat:.5f}, {nw_lon:.5f}) SE({se_lat:.5f}, {se_lon:.5f})")

    try:
        t0 = time.time()
        resp = requests.get(NYC_BUILDING_API, params=params, timeout=timeout)
        resp.raise_for_status()
        data = resp.json()
        elapsed = time.time() - t0

        n_features = len(data.get("features", []))
        print(f"  [OK] {n_features} buildings retrieved in {elapsed:.1f}s")
        return data

    except Exception as e:
        print(f"  [FAIL] API query failed: {e}")
        return None


def parse_geojson_building(
        feature: Dict) -> Optional[Tuple[List[Tuple[float, float]], float, str, Dict]]:
    """
    Parse one GeoJSON building feature from NYC Open Data.

    Args:
        feature: GeoJSON Feature dict

    Returns:
        (footprint_2d, height_m, bin_str, props) or None if invalid/too short
    """
    props = feature.get("properties", {})
    geom = feature.get("geometry", {})

    # Height extraction (field is in feet)
    raw_height = props.get("height_roof") or props.get("HEIGHTROOF") or 0
    try:
        height_m = float(raw_height) * FT_TO_M
    except (TypeError, ValueError):
        height_m = 0.0

    if height_m < MIN_BLDG_HEIGHT_M:
        return None

    bin_str = str(props.get("bin") or props.get("BIN") or "")
    geom_type = geom.get("type", "")
    coords = geom.get("coordinates", [])

    # Extract largest exterior ring from Polygon or MultiPolygon
    ring = None
    best_area = 0.0

    if geom_type == "MultiPolygon":
        for polygon in coords:
            if polygon and polygon[0] and len(polygon[0]) >= 4:
                candidate = polygon[0]
                pts = [(c[0], c[1]) for c in candidate]
                a = abs(polygon_signed_area_2d(pts))
                if a > best_area:
                    best_area = a
                    ring = candidate

    elif geom_type == "Polygon":
        if coords and coords[0] and len(coords[0]) >= 4:
            ring = coords[0]

    if ring is None or len(ring) < 4:
        return None

    # Convert geographic to local meters
    verts = [(c[0], c[1]) for c in ring]  # (lon, lat) pairs
    verts_local = [wgs84_to_local(lon, lat) for lon, lat in verts]
    verts_local = remove_duplicate_close_vertex(verts_local)

    if len(verts_local) < 3:
        return None

    # Enforce CCW winding and minimum area
    verts_local = ensure_ccw(verts_local)
    area = abs(polygon_signed_area_2d(verts_local))
    if area < MIN_BLDG_AREA_M2:
        return None

    return (verts_local, height_m, bin_str, props)


def _is_citicorp(bin_str: str, footprint: List[Tuple[float, float]]) -> bool:
    """
    Detect if a building feature represents Citicorp Center.

    Uses BIN number (authoritative) or centroid proximity (fallback).

    Args:
        bin_str: Building Identification Number
        footprint: Building footprint in local coordinates

    Returns:
        True if this is Citicorp Center
    """
    if bin_str in CITICORP_BIN_SET:
        return True

    if footprint:
        cx, cy = polygon_centroid_2d(footprint)
        if math.sqrt(cx * cx + cy * cy) < CITICORP_PROXIMITY_M:
            return True

    return False


def _infer_roof_type_from_props(props: Dict, height_m: float) -> str:
    """
    Heuristic roof type inference from building properties.

    Rules (NYC building stock):
      - High-rise (> 80m): flat parapet roof (commercial norm)
      - Pre-war (< 1940): likely gabled or hipped if < 20m
      - Residential < 15m: gabled (common in outer boroughs)
      - Residential 15-30m: hipped (brownstones/rowhouses)
      - Default: flat

    Args:
        props: Building properties dict
        height_m: Building height in meters

    Returns:
        'flat', 'gabled', 'hipped', or 'mansard'
    """
    if height_m > 80:
        return "flat"

    year_raw = props.get("cnstrct_yr") or props.get("CNSTRCT_YR")
    year = 0
    if year_raw:
        try:
            year = int(year_raw)
        except (TypeError, ValueError):
            pass

    lststatype = str(props.get("lststatype") or "").lower()
    is_residential = "residential" in lststatype or "res" in lststatype

    if height_m < 15 and is_residential and (year == 0 or year < 1960):
        return "gabled"

    if 15 <= height_m <= 35 and is_residential and (year == 0 or year < 1950):
        return "hipped"

    if 35 < height_m <= 80 and (year == 0 or year < 1900):
        return "mansard"  # Manhattan pre-1900 mid-rises often had mansard roofs

    return "flat"


# ============================================================================
# HARDCODED FALLBACK GEOMETRY
# ============================================================================

# Known buildings within the CFD domain (approximate box parameters)
# Format: (name, cx_m, cy_m, wx_m, wy_m, height_m)
HARDCODED_SURROUNDINGS = [
    # Major midtown buildings near Citicorp Center (approximate local coords)
    ("399_Park_Avenue",         -180,  -40,  55, 40, 162),   # N/W of Citicorp
    ("280_Park_Avenue",         -150,   55,  50, 35, 135),
    ("780_Third_Ave",            185,  -20,  48, 48, 175),   # E of Citicorp
    ("919_Third_Ave",            185,   65,  52, 42, 182),
    ("885_Third_Ave_Lipstick",   190,  180,  38, 38, 138),
    ("599_Lexington",            -90,   85,  43, 38, 199),
    ("731_Lexington",            -95, -110,  48, 42, 180),   # Bloomberg Tower
    ("135_E_54th_St",             10,  105,  38, 32, 150),
    ("153_E_53rd_St",             75,  -75,  33, 28, 120),
    ("Park_Avenue_Tower_S",     -390,  -45,  52, 38, 150),
    ("Park_Avenue_Tower_N",     -390,   65,  48, 32, 120),
    ("Second_Ave_Bldg_S",        395,  -55,  48, 38, 130),
    ("Second_Ave_Bldg_N",        400,   75,  42, 32, 110),
    ("E_55th_at_Lex",           -105,  175,  38, 28, 100),
    ("E_52nd_at_Lex",           -108, -185,  42, 32,  90),
    ("E_55th_at_Third",          108,  175,  32, 28,  85),
    ("E_52nd_at_Third",          112, -168,  38, 32,  95),
    ("E_53rd_mid",               -45, -105,  28, 22,  70),
    ("Greenacre_Park_Bldg",      155,  145,  45, 35,  80),
    ("St_Peters_Church",          35,   40,  30, 20,  20),  # Short church
]


def make_hardcoded_surroundings() -> List[Triangle3D]:
    """
    Generate approximate box geometry for known surrounding buildings.

    Returns:
        Combined list of STL triangles
    """
    triangles: List[Triangle3D] = []
    for _name, cx, cy, wx, wy, h in HARDCODED_SURROUNDINGS:
        triangles.extend(box_triangles(
            cx - wx / 2, cy - wy / 2, 0.0,
            cx + wx / 2, cy + wy / 2, float(h)
        ))
    return triangles


# ============================================================================
# CITICORP CENTER DETAILED GEOMETRY
# ============================================================================

def generate_citicorp_tower() -> List[Triangle3D]:
    """
    Generate the Citicorp Center tower with hand-crafted LOD2 geometry.

    Geometry specifications (Morgenstern 1995):
      Footprint:  47.85m x 47.85m (157 ft square)
      Rotation:   45 deg CCW from East (matches Manhattan grid / cardinal alignment)
      Tower body: 34.75m to 244.15m (stilt top to roof eave)
      Slanted roof: Planar surface from south vertex (244.15m) to north vertex (278.9m)
      Roof slope: 45 deg -- iconic architectural feature

    After 45 deg rotation of the unrotated square, the four corners land at:
      Corner 0 (was SW): (  0,      -R)   R = HALF_T * sqrt(2) = 33.835m
      Corner 1 (was SE): ( +R,       0)
      Corner 2 (was NE): (  0,      +R)
      Corner 3 (was NW): ( -R,       0)

    The slanted roof plane interpolates Z linearly with the Y coordinate:
      y = -R (southernmost, corner 0) -> z = TOWER_H_BASE (244.15m)
      y =  0 (east/west,  corners 1,3)-> z = (TOWER_H_BASE + TOWER_H_TOP) / 2
      y = +R (northernmost, corner 2) -> z = TOWER_H_TOP   (278.90m)

    This makes the roof plane pass through all four roof-level corner points
    forming a flat sloped quadrilateral. The slope is 45 deg in the N-S direction.

    The roof plane equation: z = TOWER_H_BASE + (y + R) / (2R) * dH
    where dH = TOWER_H_TOP - TOWER_H_BASE = 34.75m

    Returns:
        List of ((v0, v1, v2), normal) triangles (watertight solid)
    """
    triangles: List[Triangle3D] = []

    # Unrotated square corners (CCW winding when viewed from +Z)
    corners_unrot = [
        (-HALF_T, -HALF_T),   # 0: SW
        ( HALF_T, -HALF_T),   # 1: SE
        ( HALF_T,  HALF_T),   # 2: NE
        (-HALF_T,  HALF_T),   # 3: NW
    ]

    # Apply 45 deg CCW rotation
    corners = [rotate_point_2d(x, y, TOWER_ROTATION_RAD)
               for x, y in corners_unrot]

    # After rotation:
    #   corners[0] = (0,   -R)   R = HALF_T * sqrt(2)
    #   corners[1] = (+R,   0)
    #   corners[2] = (0,   +R)
    #   corners[3] = (-R,   0)

    # Compute Z elevation for each corner using linear interpolation along Y
    R = HALF_T * math.sqrt(2)  # = 33.835m (circumradius of rotated square)
    dH = TOWER_H_TOP - TOWER_H_BASE  # = 34.75m (height of slanted roof section)

    def roof_z(y_coord: float) -> float:
        """Z height of slanted roof at given Y coordinate."""
        return TOWER_H_BASE + (y_coord + R) / (2.0 * R) * dH

    # Roof corner elevations (Z)
    roof_z_vals = [roof_z(c[1]) for c in corners]
    # corners[0]: y=-R  -> z=TOWER_H_BASE  (244.15m) -- south tip (low)
    # corners[1]: y=0   -> z=midpoint       (261.53m) -- east tip (mid)
    # corners[2]: y=+R  -> z=TOWER_H_TOP    (278.90m) -- north tip (high)
    # corners[3]: y=0   -> z=midpoint       (261.53m) -- west tip (mid)

    # 3D roof corner points
    roof_pts = [(corners[i][0], corners[i][1], roof_z_vals[i]) for i in range(4)]

    # -------------------------------------------------------------------------
    # Main tower body: stilt top to roof corners
    # Build side walls manually so each wall goes from STILT_H to its correct
    # roof Z elevation (which varies per corner because of the slant).
    # -------------------------------------------------------------------------

    # Bottom cap (at STILT_H)
    tris_2d = ear_clip_triangulate(corners)
    for i0, i1, i2 in tris_2d:
        v0 = (corners[i0][0], corners[i0][1], STILT_H)
        v1 = (corners[i1][0], corners[i1][1], STILT_H)
        v2 = (corners[i2][0], corners[i2][1], STILT_H)
        triangles.append(((v0, v2, v1), (0.0, 0.0, -1.0)))  # Downward normal

    # Four trapezoidal side walls (each has 4 corners, 2 at STILT_H and 2 at roof height)
    n = len(corners)
    for i in range(n):
        j = (i + 1) % n
        # Wall base (at stilt level)
        bl = (corners[i][0], corners[i][1], STILT_H)
        br = (corners[j][0], corners[j][1], STILT_H)
        # Wall top (at slant height, varies per corner)
        tl = roof_pts[i]
        tr = roof_pts[j]

        # Outward wall normal (horizontal, perpendicular to wall edge, CCW polygon)
        dx = corners[j][0] - corners[i][0]
        dy = corners[j][1] - corners[i][1]
        wall_n = vec3_normalize((dy, -dx, 0.0))

        # Two triangles per wall quad
        triangles.append(((bl, br, tr), wall_n))
        triangles.append(((bl, tr, tl), wall_n))

    # -------------------------------------------------------------------------
    # Slanted roof surface (planar quadrilateral split into 2 triangles)
    # The four roof corners form a planar quadrilateral because Z = f(Y).
    # Split: (0,1,2) and (0,2,3)  -- maintain CCW winding for upward normal
    # -------------------------------------------------------------------------

    # Verify planarity (all 4 points coplanar -- they are, by construction)
    r0, r1, r2, r3 = roof_pts

    # Normal of the slanted roof plane (should point upward with a northward tilt)
    n_roof = face_normal(r0, r1, r2)
    # For our geometry: slope is in Y direction, normal should be (0, -sin45, cos45)
    # = (0, -0.707, 0.707) approximately

    triangles.append(((r0, r1, r2), n_roof))
    triangles.append(((r0, r2, r3), n_roof))

    return triangles


def generate_citicorp_stilts() -> List[Triangle3D]:
    """
    Generate four concrete stilt columns supporting Citicorp Center.

    Structural specifications (Morgenstern 1995):
      Count:    4 columns
      Section:  7.32m × 7.32m (24 ft square concrete)
      Height:   0 to 34.75m (10 stories above grade)
      Location: Midpoints of the tower's four faces (not at corners)
      Rotation: Same 45° as tower (columns support face-centers)

    This unconventional column placement (midpoint vs corner) was what
    made the original Citicorp design structurally vulnerable to quartering
    winds — the key event in the 1978 crisis.

    Returns:
        List of ((v0, v1, v2), normal) triangles (watertight solid)
    """
    triangles: List[Triangle3D] = []

    for cx_u, cy_u in STILTS_UNROTATED:
        # Rotate stilt center position by 45°
        cx, cy = rotate_point_2d(cx_u, cy_u, TOWER_ROTATION_RAD)

        # Stilt corners (small square, also rotated)
        stilt_fp_unrot = [
            (cx_u - HALF_S, cy_u - HALF_S),
            (cx_u + HALF_S, cy_u - HALF_S),
            (cx_u + HALF_S, cy_u + HALF_S),
            (cx_u - HALF_S, cy_u + HALF_S),
        ]
        stilt_fp = [rotate_point_2d(x, y, TOWER_ROTATION_RAD)
                    for x, y in stilt_fp_unrot]

        triangles.extend(extrude_flat_roof(stilt_fp, 0.0, STILT_H))

    return triangles


# ============================================================================
# STL FILE I/O
# ============================================================================

def write_stl_ascii(
        filepath: str,
        solid_name: str,
        triangles: List[Triangle3D]) -> int:
    """
    Write an ASCII STL file.

    ASCII STL is human-readable and compatible with all STL parsers.
    For large meshes (> 100k triangles), prefer binary STL for size.

    Format:
      solid <name>
        facet normal nx ny nz
          outer loop
            vertex x y z
            vertex x y z
            vertex x y z
          endloop
        endfacet
      endsolid <name>

    Args:
        filepath: Output file path
        solid_name: Name embedded in STL header
        triangles: List of ((v0, v1, v2), normal) tuples

    Returns:
        Number of triangles written
    """
    with open(filepath, 'w', encoding='ascii') as f:
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


def write_stl_binary(
        filepath: str,
        solid_name: str,
        triangles: List[Triangle3D]) -> int:
    """
    Write a binary STL file (50-80% smaller than ASCII).

    Binary STL format:
      80 bytes: Header string (arbitrary text, NOT "solid")
      4 bytes:  uint32 triangle count
      For each triangle (50 bytes):
        12 bytes: float32[3] normal
        12 bytes: float32[3] vertex 0
        12 bytes: float32[3] vertex 1
        12 bytes: float32[3] vertex 2
         2 bytes: uint16 attribute byte count (0)

    Args:
        filepath: Output file path
        solid_name: Text for header (80 bytes max)
        triangles: List of ((v0, v1, v2), normal) tuples

    Returns:
        Number of triangles written
    """
    n_tri = len(triangles)
    header = f"Binary STL: {solid_name}".encode('ascii')[:80].ljust(80, b'\x00')

    with open(filepath, 'wb') as f:
        f.write(header)
        f.write(struct.pack('<I', n_tri))

        for verts, normal in triangles:
            # Normal vector
            f.write(struct.pack('<fff', normal[0], normal[1], normal[2]))
            # Three vertices
            for vx, vy, vz in verts:
                f.write(struct.pack('<fff', vx, vy, vz))
            # Attribute byte count (unused)
            f.write(struct.pack('<H', 0))

    return n_tri


def validate_stl_trimesh(filepath: str) -> Dict[str, Any]:
    """
    Validate STL mesh using trimesh (if available).

    Checks:
      - Is watertight (no open edges)
      - Is winding-consistent (all normals point outward)
      - Volume is positive (solid, not hollow)
      - No degenerate triangles

    Args:
        filepath: Path to STL file

    Returns:
        Dict with validation results, or {'error': ...} if trimesh unavailable
    """
    if not HAS_TRIMESH:
        return {"error": "trimesh not available"}

    try:
        mesh = trimesh.load(filepath, force='mesh')
        return {
            "filepath": filepath,
            "n_triangles": len(mesh.faces),
            "n_vertices": len(mesh.vertices),
            "is_watertight": mesh.is_watertight,
            "is_winding_consistent": mesh.is_winding_consistent,
            "volume_m3": float(mesh.volume) if mesh.is_watertight else None,
            "bounding_box": mesh.bounds.tolist(),
        }
    except Exception as e:
        return {"error": str(e)}


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main() -> None:
    """Main execution function."""
    # -------------------------------------------------------------------------
    # Header
    # -------------------------------------------------------------------------
    print()
    print("=" * 80)
    print("  Citicorp Center LOD2 STL Generator")
    print("  OR 750 Reliability, Safety, and Risk — GMU PhD Program")
    print("=" * 80)
    print()

    # -------------------------------------------------------------------------
    # Parse arguments
    # -------------------------------------------------------------------------
    args = sys.argv[1:]
    offline_mode  = "--offline"  in args
    validate_mode = "--validate" in args
    binary_stl    = "--binary"   in args
    verbose       = "--verbose"  in args

    cityjson_file: Optional[str] = None
    if "--cityjson" in args:
        idx = args.index("--cityjson")
        if idx + 1 < len(args):
            cityjson_file = args[idx + 1]
        else:
            print("Error: --cityjson requires a file path argument")
            sys.exit(1)

    outdir_arg: Optional[str] = None
    if "--outdir" in args:
        idx = args.index("--outdir")
        if idx + 1 < len(args):
            outdir_arg = args[idx + 1]

    # -------------------------------------------------------------------------
    # Output directory
    # -------------------------------------------------------------------------
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if outdir_arg:
        out_dir = os.path.abspath(outdir_arg)
    else:
        out_dir = os.path.join(script_dir, "constant", "triSurface")

    os.makedirs(out_dir, exist_ok=True)

    print(f"Output directory: {out_dir}")
    print()

    # Print capability summary
    print("Available libraries:")
    print(f"  requests:  {'YES' if HAS_REQUESTS else 'NO  (pip install requests)'}")
    print(f"  pyproj:    {'YES' if HAS_PYPROJ else 'NO  (pip install pyproj)'}")
    print(f"  shapely:   {'YES' if HAS_SHAPELY else 'NO  (pip install shapely)'}")
    print(f"  trimesh:   {'YES' if HAS_TRIMESH else 'NO  (pip install trimesh)'}")
    print()

    write_fn = write_stl_binary if binary_stl else write_stl_ascii
    stl_ext_note = "(binary)" if binary_stl else "(ASCII)"

    domain_bbox = (
        DOMAIN_BOUNDS["xmin"], DOMAIN_BOUNDS["xmax"],
        DOMAIN_BOUNDS["ymin"], DOMAIN_BOUNDS["ymax"]
    )

    total_triangles = 0
    all_outputs: Dict[str, int] = {}

    # =========================================================================
    # STEP 1: Citicorp Tower
    # =========================================================================
    print("Step 1: Citicorp Tower (hand-crafted LOD2 geometry)")
    print("-" * 60)

    tower_tris = generate_citicorp_tower()
    tower_path = os.path.join(out_dir, "citicorp_tower.stl")
    n = write_fn(tower_path, "citicorp_tower", tower_tris)
    total_triangles += n
    all_outputs["citicorp_tower.stl"] = n

    print(f"  [OK] citicorp_tower.stl {stl_ext_note}: {n} triangles")
    print(f"       Footprint: {TOWER_W:.2f}m x {TOWER_W:.2f}m (157 ft square)")
    print(f"       Height range: {STILT_H:.2f}m to {TOWER_H_TOP:.2f}m")
    print(f"       Roof eave: {TOWER_H_BASE:.2f}m (801 ft), peak: {TOWER_H_TOP:.2f}m (915 ft)")
    print(f"       Rotation: {TOWER_ROTATION_DEG:.0f}deg from cardinal (Manhattan grid)")
    print()

    # =========================================================================
    # STEP 2: Citicorp Stilts
    # =========================================================================
    print("Step 2: Citicorp Stilts (hand-crafted geometry)")
    print("-" * 60)

    stilt_tris = generate_citicorp_stilts()
    stilt_path = os.path.join(out_dir, "citicorp_stilts.stl")
    n = write_fn(stilt_path, "citicorp_stilts", stilt_tris)
    total_triangles += n
    all_outputs["citicorp_stilts.stl"] = n

    print(f"  [OK] citicorp_stilts.stl {stl_ext_note}: {n} triangles (4 columns)")
    print(f"       Column section: {STILT_W:.2f}m x {STILT_W:.2f}m (24 ft square)")
    print(f"       Height: 0m to {STILT_H:.2f}m (10 stories)")
    print(f"       Location: Face midpoints at 45deg rotation")
    print(f"       Note: Midpoint placement caused 1978 structural vulnerability")
    print()

    # =========================================================================
    # STEP 3: Surrounding Buildings
    # =========================================================================
    print("Step 3: Surrounding Buildings")
    print("-" * 60)

    surr_tris: List[Triangle3D] = []
    n_buildings = 0
    n_skipped = 0
    data_source_used = "None"
    tallest_neighbor = ("", 0.0)

    # ------------------------------------------
    # Attempt 1: CityJSON LOD2 file
    # ------------------------------------------
    if cityjson_file:
        if not os.path.exists(cityjson_file):
            print(f"  [WARN] CityJSON file not found: {cityjson_file}")
        else:
            print(f"  Attempting CityJSON LOD2 parse: {cityjson_file}")
            parser = CityJSONParser(cityjson_file)
            if parser.load():
                data_source_used = f"CityJSON: {os.path.basename(cityjson_file)}"
                lod_counts: Dict[str, int] = {}
                roof_type_counts: Dict[str, int] = {}
                n_lod2_roofs = 0

                for bldg in parser.iter_buildings(domain_bbox=domain_bbox):
                    # Skip Citicorp
                    if _is_citicorp(
                        str(bldg["attributes"].get("bin", "")),
                        bldg["footprint"]
                    ):
                        print(f"  [INFO] Identified Citicorp in CityJSON, skipping")
                        n_skipped += 1
                        continue

                    if bldg["height_m"] < MIN_BLDG_HEIGHT_M:
                        n_skipped += 1
                        continue

                    surr_tris.extend(bldg["triangles"])
                    n_buildings += 1

                    lod = bldg["lod"]
                    lod_counts[lod] = lod_counts.get(lod, 0) + 1

                    rt = bldg["roof_type"]
                    roof_type_counts[rt] = roof_type_counts.get(rt, 0) + 1

                    if bldg["roof_type"] != "flat" and bldg["roof_type"] != "unknown":
                        n_lod2_roofs += 1

                    if bldg["height_m"] > tallest_neighbor[1]:
                        tallest_neighbor = (bldg["id"], bldg["height_m"])

                    if verbose:
                        print(f"    Building {bldg['id']}: "
                              f"h={bldg['height_m']:.1f}m, "
                              f"roof={bldg['roof_type']}, "
                              f"LOD={bldg['lod']}, "
                              f"{len(bldg['triangles'])} tris")

                if n_buildings > 0:
                    print(f"  [OK] Parsed {n_buildings} buildings from CityJSON")
                    print(f"       LOD distribution: {lod_counts}")
                    print(f"       Roof types: {roof_type_counts}")
                    print(f"       Buildings with non-flat LOD2 roofs: {n_lod2_roofs}")
                    print(f"       Skipped: {n_skipped}")
                else:
                    print(f"  [WARN] No buildings extracted from CityJSON")
                    print(f"         Check EPSG code, domain bounds, or file content")
                    surr_tris = []  # Reset for fallback

    # ------------------------------------------
    # Attempt 2: NYC Open Data API (LOD1 fallback)
    # ------------------------------------------
    if not surr_tris and not offline_mode:
        print("  Attempting NYC Open Data Building Footprints API (LOD1)...")
        geojson = fetch_nyc_open_data()

        if geojson:
            features = geojson.get("features", [])
            data_source_used = "NYC Open Data Building Footprints (LOD1 — flat roofs)"
            roof_type_counts: Dict[str, int] = {}

            for feat in features:
                parsed = parse_geojson_building(feat)
                if parsed is None:
                    n_skipped += 1
                    continue

                footprint, height_m, bin_str, props = parsed

                if _is_citicorp(bin_str, footprint):
                    print(f"  [INFO] Identified Citicorp (BIN={bin_str}), using hand-crafted model")
                    n_skipped += 1
                    continue

                # Infer roof type and generate geometry
                roof_type = _infer_roof_type_from_props(props, height_m)
                roof_type_counts[roof_type] = roof_type_counts.get(roof_type, 0) + 1

                # Select roof geometry
                if roof_type == "gabled":
                    # Simplified gabled: eave at 90% height, peak at 100%
                    tris = extrude_gabled_roof(
                        footprint, 0.0,
                        height_m * 0.88,
                        height_m
                    )
                elif roof_type == "hipped":
                    tris = extrude_hipped_roof(
                        footprint, 0.0,
                        height_m * 0.88,
                        height_m
                    )
                elif roof_type == "mansard":
                    tris = extrude_mansard_roof(
                        footprint, 0.0,
                        height_m * 0.75,
                        height_m,
                        setback_dist=max(2.0, min(5.0, height_m * 0.05))
                    )
                else:
                    tris = extrude_flat_roof(footprint, 0.0, height_m)

                if not tris:
                    n_skipped += 1
                    continue

                surr_tris.extend(tris)
                n_buildings += 1

                if height_m > tallest_neighbor[1]:
                    tallest_neighbor = (bin_str, height_m)

                if verbose:
                    print(f"    BIN {bin_str}: h={height_m:.1f}m, roof={roof_type}, "
                          f"{len(tris)} tris")

            if n_buildings > 0:
                print(f"  [OK] Processed {n_buildings} buildings from NYC Open Data")
                print(f"       Skipped: {n_skipped}")
                print(f"       Roof types: {roof_type_counts}")
            else:
                print(f"  [WARN] No buildings from NYC Open Data")
        else:
            print(f"  [FAIL] NYC Open Data API unavailable")

    # ------------------------------------------
    # Attempt 3: Hardcoded fallback geometry
    # ------------------------------------------
    if not surr_tris:
        reason = "--offline" if offline_mode else "all data sources failed"
        print(f"  Using hardcoded surroundings geometry ({reason})")
        surr_tris = make_hardcoded_surroundings()
        n_buildings = len(HARDCODED_SURROUNDINGS)
        data_source_used = f"Hardcoded geometry ({n_buildings} buildings)"

    # Write surrounding buildings
    surr_path = os.path.join(out_dir, "surroundings_lod2.stl")
    n = write_fn(surr_path, "surroundings_lod2", surr_tris)
    total_triangles += n
    all_outputs["surroundings_lod2.stl"] = n

    print()
    print(f"  [OK] surroundings_lod2.stl {stl_ext_note}: {n} triangles, "
          f"{n_buildings} buildings")
    print(f"       Data source: {data_source_used}")
    if tallest_neighbor[1] > 0:
        print(f"       Tallest neighbor: {tallest_neighbor[0]}, "
              f"{tallest_neighbor[1]:.1f}m ({tallest_neighbor[1]/FT_TO_M:.0f} ft)")
    print()

    # =========================================================================
    # Optional: STL validation
    # =========================================================================
    if validate_mode and HAS_TRIMESH:
        print("STL Validation (trimesh)")
        print("-" * 60)
        for fname in ("citicorp_tower.stl", "citicorp_stilts.stl", "surroundings_lod2.stl"):
            fpath = os.path.join(out_dir, fname)
            results = validate_stl_trimesh(fpath)
            if "error" in results:
                print(f"  {fname}: validation failed — {results['error']}")
            else:
                print(f"  {fname}:")
                print(f"    triangles:   {results['n_triangles']}")
                print(f"    watertight:  {results['is_watertight']}")
                print(f"    winding OK:  {results['is_winding_consistent']}")
                if results.get("volume_m3") is not None:
                    print(f"    volume:      {results['volume_m3']:.1f} m3")
        print()

    # =========================================================================
    # Summary
    # =========================================================================
    print("=" * 80)
    print("  Summary")
    print("=" * 80)
    print(f"Total triangles: {total_triangles}")
    print(f"Output directory: {out_dir}")
    print()
    print("Output files:")
    for fname, n_tris in all_outputs.items():
        size_kb = os.path.getsize(os.path.join(out_dir, fname)) / 1024
        print(f"  {fname:<35} {n_tris:>6} triangles  {size_kb:>8.1f} KB")
    print()
    print("Next steps for OpenFOAM:")
    print("  1. Verify geometry in ParaView or MeshLab")
    print("     paraview constant/triSurface/citicorp_tower.stl")
    print()
    print("  2. Run snappyHexMesh:")
    print("     blockMesh && snappyHexMesh -overwrite")
    print()
    print("  3. Check mesh quality:")
    print("     checkMesh -allGeometry -allTopology")
    print()
    print("  4. Run CFD solver:")
    print("     foamRun -solver incompressibleFluid  # Foundation v13")
    print("     simpleFoam -parallel                 # ESI v2512")
    print()
    print("="*80)
    print("  Data Source Documentation")
    print("="*80)
    print()
    print("TUM NYC LOD2 CityGML (Recommended for LOD2 detail):")
    print("  GitHub:  https://github.com/georocket/new-york-city-model-enhanced")
    print("  Format:  CityGML 2.0 with PLUTO integration (~20 GB full city)")
    print("  Tool:    pip install cjio")
    print("  Convert: cjio input.gml upgrade save output.city.json")
    print("  Extract: cjio output.city.json subset --bbox ... save citicorp.city.json")
    print()
    print("TU Delft 3D BAG (CityJSON format, Netherlands primary):")
    print("  URL:    https://3d.bk.tudelft.nl/opendata/3dbag/")
    print("  Tiles:  https://3d.bk.tudelft.nl/opendata/3dbag/v2/tiles/")
    print("  EPSG:   28992 (RD New Dutch Grid)")
    print("  Note:   NYC coverage may be limited")
    print()
    print("NYC Open Data Building Footprints (LOD1 — reliable fallback):")
    print("  API:    https://data.cityofnewyork.us/resource/5zhs-2jue.geojson")
    print("  Docs:   https://dev.socrata.com/foundry/data.cityofnewyork.us/5zhs-2jue")
    print("  EPSG:   4326 (WGS84) via Socrata GeoJSON endpoint")
    print("  Token:  https://dev.socrata.com/ (free, higher rate limits)")
    print()
    print("NYC PLUTO (building attributes, optional join):")
    print("  URL:    https://www.nyc.gov/site/planning/data-maps/open-data/dwn-pluto-mappluto.page")
    print("  Join:   BIN (Building Identification Number)")
    print()
    print("CityJSON Schema Reference:")
    print("  Spec:        https://www.cityjson.org/specs/1.1.3/")
    print("  Tutorials:   https://www.cityjson.org/tutorials/")
    print("  Validator:   https://validator.cityjson.org/")
    print("  Python API:  pip install cjio")
    print()
    print("Coordinate Systems:")
    print("  EPSG:4326  — WGS84 geographic (lon, lat) — NYC Open Data default")
    print("  EPSG:2263  — NAD83 NY State Plane LI (US Survey Feet) — NYC DOITT")
    print("  EPSG:32618 — UTM Zone 18N (meters) — good local projection for NYC")
    print("  EPSG:28992 — RD New Dutch Grid (meters) — TU Delft 3D BAG")
    print("  Local:     — Cartesian, origin Citicorp, X=East, Y=North, Z=Up")
    print()
    print("Conversion commands (pyproj):")
    print("  from pyproj import Transformer")
    print("  t = Transformer.from_crs('EPSG:28992', 'EPSG:4326', always_xy=True)")
    print("  lon, lat = t.transform(x_rdnew, y_rdnew)")
    print()
    print("Usage:")
    print("  python generate_stl_lod2.py                          # Online, NYC API")
    print("  python generate_stl_lod2.py --cityjson file.json     # LOD2 CityJSON")
    print("  python generate_stl_lod2.py --offline                # Hardcoded only")
    print("  python generate_stl_lod2.py --binary                 # Binary STL output")
    print("  python generate_stl_lod2.py --validate               # Validate with trimesh")
    print("  python generate_stl_lod2.py --verbose                # Per-building output")
    print("  python generate_stl_lod2.py --outdir /path/to/dir    # Custom output dir")
    print()


if __name__ == "__main__":
    main()
