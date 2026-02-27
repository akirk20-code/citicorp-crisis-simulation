#!/usr/bin/env python3
"""Generate STL geometry using CityGML (primary) + Socrata API (fallback).

CityGML provides actual building surfaces (LoD2 roof shapes, wall setbacks)
from the NYC 3D Building Model (2014 aerial survey). For buildings built
after 2014, the Socrata Building Footprints API fills the gap with extruded
footprints (LoD1).

Usage:
    python generate_stl_hybrid.py                                    # Citicorp (default)
    python generate_stl_hybrid.py --offline                          # CityGML only, no API
    python generate_stl_hybrid.py --year 1978                        # historical filter
    python generate_stl_hybrid.py --center-lat 40.7489 \
        --center-lon -73.9680 --name un_hq --radius 300 --no-target  # UN HQ area

Dependencies: Python 3.8+ (stdlib only for CityGML; requests for Socrata fallback)
Output: constant/triSurface/{name}_surroundings.stl (+ tower/stilts for Citicorp)
"""

import argparse
import math
import os
import re
import struct
import sys
import time
import zipfile

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

try:
    import urllib.request
    HAS_URLLIB = True
except ImportError:
    HAS_URLLIB = False


# ============================================================
# CITICORP GEOMETRY (meters)
# ============================================================
# Tower rotation: Manhattan grid is 28.7° east of true north.
# In EPSG:2263 (X=east, Y=north), building faces run at 61.3° and 151.3°.
# Measured from CityGML LoD2 wall edges of BIN 1036474.
TOWER_ANGLE_DEG = 28.7   # degrees east of north (= 61.3° from X-axis)

TOWER_W = 47.85       # 157 ft -- square plan
TOWER_H_TOP = 278.9   # 915 ft -- roof height
STILT_H = 34.75       # 114 ft -- stilt height (tower bottom)
# Stilt dimensions refined from photo panel counting:
#   facade = 33 panel-widths across 157 ft → panel = 48/33 = 1.4545 m
#   face width = 5 panels = 7.27 m, depth = 4.5 panels = 6.55 m
STILT_FACE = 7.27     # 5 panels along face (published "24 ft" ≈ 7.32 m)
STILT_DEPTH = 6.55    # 4.5 panels toward center
HALF_T = TOWER_W / 2  # 23.925 m

# Stilt centers + orientation: (cx, cy, face_axis)
#   face_axis='x' → wider in X (south/north faces)
#   face_axis='y' → wider in Y (east/west faces)
#   Stilt outer edge is flush with tower face → center offset inward by depth/2
STILT_INSET = STILT_DEPTH / 2  # 3.275 m inward from face
STILTS = [
    (0,                   -HALF_T + STILT_INSET, 'x'),   # South face
    (HALF_T - STILT_INSET, 0,                    'y'),   # East face
    (0,                    HALF_T - STILT_INSET,  'x'),   # North face
    (-HALF_T + STILT_INSET, 0,                   'y'),   # West face
]

# Center core (elevator/mechanical core visible in street-level photos)
# Square plan with chamfered corners (~8-sided octagon)
# Dimensions from photo panel counting (panel = 48/33 = 1.4545 m):
#   bounding box ≈ 14.5 panels = 21.1 m (each side)
#   flat edges   ≈ 10 panels   = 14.5 m (chamfer cuts corners)
#   chamfer      = (21.1 - 14.5) / 2 = 3.3 m
CORE_W = 21.09        # 14.5 panels × 1.4545 m = 21.09 m bounding box
CORE_D = 21.09        # same -- square plan
CORE_CHAMFER = 3.27   # corner cut so each flat face = 10 panels = 14.545 m

# BINs to exclude from surroundings (all belong to Citicorp complex)
CITICORP_BIN_SET = {"1036474", "1035879", "1087931"}
CITICORP_PROXIMITY_M = 40.0


# ============================================================
# COORDINATE SYSTEM CONSTANTS
# ============================================================
# Default center (Citicorp) — overridden by --center-lat/--center-lon
DEFAULT_CENTER_LAT = 40.7579
DEFAULT_CENTER_LON = -73.9690

# Unit conversions (constant for any NYC location)
US_SURVEY_FT_TO_M = 1200.0 / 3937.0    # 0.30480061 m/ft
M_TO_US_SURVEY_FT = 3937.0 / 1200.0    # 3.28083333 ft/m
FT_TO_M = 0.3048                        # International foot (Socrata heights)

# EPSG:2263 projection constants (constant for Manhattan)
FT_PER_DEG_LAT = 364173.0
FT_PER_DEG_LON = 276540.0

# These are set dynamically in main() from CLI args:
CENTER_LAT = DEFAULT_CENTER_LAT
CENTER_LON = DEFAULT_CENTER_LON
CENTER_X_FT = 992506.0   # recomputed in main()
CENTER_Y_FT = 215620.0   # recomputed in main()
M_PER_DEG_LAT = 111320.0
M_PER_DEG_LON = M_PER_DEG_LAT * math.cos(math.radians(DEFAULT_CENTER_LAT))
DOMAIN_X_MIN = -360
DOMAIN_X_MAX = 360
DOMAIN_Y_MIN = -360
DOMAIN_Y_MAX = 360

# Filtering
MIN_HEIGHT_M = 5.0
MIN_AREA_M2 = 10.0

# Socrata API
NYC_API_URL = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"

# DA district approximate bounding boxes (EPSG:2263 feet)
DA_ENVELOPES = {
    'DA1':  (970000, 990000, 188000, 205000),
    'DA2':  (975000, 995000, 193000, 210000),
    'DA3':  (978000, 998000, 196000, 213000),
    'DA4':  (972000, 992000, 196000, 216000),
    'DA5':  (968000, 988000, 198000, 220000),
    'DA6':  (978000, 1000000, 200000, 216000),
    'DA7':  (982000, 1005000, 210000, 228000),
    'DA8':  (997000, 1015000, 189000, 206000),
    'DA9':  (995000, 1015000, 205000, 225000),
    'DA10': (1000000, 1020000, 220000, 240000),
    'DA11': (985000, 1005000, 225000, 245000),
    'DA12': (978000, 1003000, 194000, 221000),
    'DA13': (978000, 1000000, 240000, 260000),
    'DA14': (1000000, 1025000, 165000, 195000),
    'DA15': (1010000, 1035000, 185000, 210000),
    'DA16': (1020000, 1045000, 205000, 230000),
    'DA17': (1030000, 1060000, 230000, 255000),
    'DA18': (1040000, 1070000, 245000, 270000),
    'DA19': (1050000, 1080000, 195000, 230000),
    'DA20': (1000000, 1035000, 240000, 275000),
}

GITHUB_RELEASE_URL = (
    "https://github.com/georocket/new-york-city-model-enhanced/"
    "releases/download/20v5/{da}_3D_Buildings_Merged.gml.zip"
)

# Hardcoded fallback (18 approximate buildings, used when both sources fail)
FALLBACK_SURROUNDINGS = [
    ("399_Park_Ave",       -200,  -50,  60, 45, 162),
    ("280_Park_Ave",       -200,   50,  55, 40, 135),
    ("780_Third_Ave",       200,  -25,  50, 50, 175),
    ("919_Third_Ave",       200,   60,  55, 45, 182),
    ("885_Third_Ave",       200,  180,  35, 35, 138),
    ("599_Lexington",      -100,   80,  45, 40, 199),
    ("731_Lexington",      -100, -100,  50, 45, 180),
    ("135_E_54th",            0,  100,  40, 35, 150),
    ("153_E_53rd",           80,  -80,  35, 30, 120),
    ("Park_Ave_Tower_S",   -400,  -50,  55, 40, 150),
    ("Park_Ave_Tower_N",   -400,   60,  50, 35, 120),
    ("Second_Ave_S",        400,  -60,  50, 40, 130),
    ("Second_Ave_N",        400,   70,  45, 35, 110),
    ("E_55th_Lex",         -110,  170,  40, 30, 100),
    ("E_52nd_Lex",         -110, -180,  45, 35,  90),
    ("E_55th_Third",        110,  170,  35, 30,  85),
    ("E_52nd_Third",        110, -170,  40, 35,  95),
    ("E_53rd_mid",          -50, -100,  30, 25,  70),
]


# ============================================================
# COORDINATE TRANSFORMS
# ============================================================

def latlon_to_stateplane(lat, lon):
    """Convert WGS84 lat/lon to approximate EPSG:2263 (NY State Plane feet)."""
    x = CENTER_X_FT + (lon - CENTER_LON) * FT_PER_DEG_LON
    y = CENTER_Y_FT + (lat - CENTER_LAT) * FT_PER_DEG_LAT
    return x, y


def wgs84_to_local(lon, lat):
    """Convert WGS84 (lon, lat) to local meters centered on Citicorp.

    Uses the same scale factors as the CityGML (EPSG:2263) path to ensure
    Socrata and CityGML buildings align exactly at the center point.
    Route: WGS84 -> approx EPSG:2263 feet -> local meters.
    """
    x = (lon - CENTER_LON) * FT_PER_DEG_LON * US_SURVEY_FT_TO_M
    y = (lat - CENTER_LAT) * FT_PER_DEG_LAT * US_SURVEY_FT_TO_M
    return x, y


def citygml_to_local(x_ft, y_ft, z_ft, z_ground_ft):
    """Convert EPSG:2263 feet to local meters centered at Citicorp."""
    x_m = (x_ft - CENTER_X_FT) * US_SURVEY_FT_TO_M
    y_m = (y_ft - CENTER_Y_FT) * US_SURVEY_FT_TO_M
    z_m = (z_ft - z_ground_ft) * US_SURVEY_FT_TO_M
    return (x_m, y_m, z_m)


# ============================================================
# POLYGON / TRIANGLE GEOMETRY
# ============================================================

def _normalize(v):
    length = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if length < 1e-15:
        return (0.0, 0.0, 1.0)
    return (v[0] / length, v[1] / length, v[2] / length)


def rotate_triangles(tris, angle_deg):
    """Rotate triangle list around Z axis by angle_deg (CCW from X-axis).

    Each tri is ((v0, v1, v2), normal).
    Returns new list with rotated vertices and recomputed normals.
    """
    theta = math.radians(90.0 - angle_deg)  # convert "degrees E of N" to CCW from X
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)

    def rot(v):
        return (v[0]*cos_t - v[1]*sin_t,
                v[0]*sin_t + v[1]*cos_t,
                v[2])

    rotated = []
    for (v0, v1, v2), _ in tris:
        rv0, rv1, rv2 = rot(v0), rot(v1), rot(v2)
        normal = compute_normal(rv0, rv1, rv2)
        rotated.append(((rv0, rv1, rv2), normal))
    return rotated


def compute_normal(v0, v1, v2):
    """Compute triangle normal via cross product."""
    e1 = (v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2])
    e2 = (v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2])
    nx = e1[1]*e2[2] - e1[2]*e2[1]
    ny = e1[2]*e2[0] - e1[0]*e2[2]
    nz = e1[0]*e2[1] - e1[1]*e2[0]
    length = math.sqrt(nx*nx + ny*ny + nz*nz)
    if length < 1e-12:
        return (0.0, 0.0, 1.0)
    return (nx/length, ny/length, nz/length)


def polygon_area_2d(verts):
    """Signed area of a 2D polygon. Positive = CCW."""
    n = len(verts)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += verts[i][0] * verts[j][1]
        area -= verts[j][0] * verts[i][1]
    return area / 2.0


def point_in_triangle_2d(p, a, b, c):
    """Check if 2D point p is strictly inside triangle (a, b, c)."""
    def sign(p1, p2, p3):
        return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])
    d1 = sign(p, a, b)
    d2 = sign(p, b, c)
    d3 = sign(p, c, a)
    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (has_neg and has_pos)


def triangulate_polygon_earclip(verts):
    """Ear-clipping triangulation of a simple 2D polygon (CCW winding).

    Returns list of (i, j, k) index triples.
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
            nx_idx = idx[(i + 1) % m]

            ax, ay = verts[p]
            bx, by = verts[c]
            cx, cy = verts[nx_idx]

            cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            if cross <= 1e-12:
                continue

            inside = False
            for j in range(m):
                if j in ((i - 1) % m, i, (i + 1) % m):
                    continue
                if point_in_triangle_2d(verts[idx[j]],
                                        verts[p], verts[c], verts[nx_idx]):
                    inside = True
                    break
            if inside:
                continue

            tris.append((p, c, nx_idx))
            idx.pop(i)
            ear_found = True
            break

        if not ear_found:
            for i in range(1, len(idx) - 1):
                tris.append((idx[0], idx[i], idx[i + 1]))
            break

    if len(idx) == 3:
        tris.append((idx[0], idx[1], idx[2]))

    return tris


def fan_triangulate(vertices):
    """Fan triangulation from first vertex (for 3D CityGML polygons)."""
    triangles = []
    if len(vertices) < 3:
        return triangles
    v0 = vertices[0]
    for i in range(1, len(vertices) - 1):
        triangles.append((v0, vertices[i], vertices[i + 1]))
    return triangles


def box_triangles(xmin, ymin, zmin, xmax, ymax, zmax):
    """12 triangles for a watertight axis-aligned box, outward normals."""
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
        ((0, 2, 1), ( 0,  0, -1)),
        ((0, 3, 2), ( 0,  0, -1)),
        ((4, 5, 6), ( 0,  0,  1)),
        ((4, 6, 7), ( 0,  0,  1)),
        ((0, 1, 5), ( 0, -1,  0)),
        ((0, 5, 4), ( 0, -1,  0)),
        ((2, 3, 7), ( 0,  1,  0)),
        ((2, 7, 6), ( 0,  1,  0)),
        ((0, 4, 7), (-1,  0,  0)),
        ((0, 7, 3), (-1,  0,  0)),
        ((1, 2, 6), ( 1,  0,  0)),
        ((1, 6, 5), ( 1,  0,  0)),
    ]
    return [(tuple(v[i] for i in idx), n) for idx, n in faces]


def make_chamfered_rect(cx, cy, w, d, chamfer):
    """Create 8-sided polygon: rectangle with chamfered corners.

    Returns CCW 2D vertices centered at (cx, cy).
    w = width (along x), d = depth (along y), chamfer = corner cut size.
    """
    hw = w / 2
    hd = d / 2
    c = chamfer
    # 8 vertices, starting from top-right going CCW
    return [
        (cx + hw,       cy + hd - c),   # right side, near top-right corner
        (cx + hw,       cy - hd + c),   # right side, near bottom-right corner
        (cx + hw - c,   cy - hd),       # bottom side, near bottom-right corner
        (cx - hw + c,   cy - hd),       # bottom side, near bottom-left corner
        (cx - hw,       cy - hd + c),   # left side, near bottom-left corner
        (cx - hw,       cy + hd - c),   # left side, near top-left corner
        (cx - hw + c,   cy + hd),       # top side, near top-left corner
        (cx + hw - c,   cy + hd),       # top side, near top-right corner
    ]


def extrude_polygon_to_triangles(verts_2d, z_base, z_top):
    """Extrude a CCW 2D polygon into a watertight 3D solid.

    Returns list of ((v0, v1, v2), normal) triangles.
    """
    n = len(verts_2d)
    if n < 3 or z_top <= z_base:
        return []

    tri_idx = triangulate_polygon_earclip(verts_2d)
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

        dx = verts_2d[j][0] - verts_2d[i][0]
        dy = verts_2d[j][1] - verts_2d[i][1]
        wall_n = _normalize((dy, -dx, 0.0))

        result.append(((bl, br, tr), wall_n))
        result.append(((bl, tr, tl), wall_n))

    return result


# ============================================================
# STL WRITERS
# ============================================================

def write_binary_stl(filepath, solid_name, triangles):
    """Write binary STL. triangles = [((v0,v1,v2), normal), ...]"""
    n_tri = len(triangles)
    with open(filepath, 'wb') as f:
        header = solid_name.encode('ascii')[:80].ljust(80, b'\0')
        f.write(header)
        f.write(struct.pack('<I', n_tri))
        for verts, normal in triangles:
            f.write(struct.pack('<fff', *normal))
            for v in verts:
                f.write(struct.pack('<fff', *v))
            f.write(struct.pack('<H', 0))
    return n_tri


def write_ascii_stl(filepath, solid_name, triangles):
    """Write ASCII STL for inspection."""
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


def read_binary_stl(filepath):
    """Read binary STL file and return list of ((v0, v1, v2), normal) tuples."""
    triangles = []
    with open(filepath, 'rb') as f:
        f.read(80)  # header
        n_tri = struct.unpack('<I', f.read(4))[0]
        for _ in range(n_tri):
            data = struct.unpack('<12fH', f.read(50))
            normal = (data[0], data[1], data[2])
            v0 = (data[3], data[4], data[5])
            v1 = (data[6], data[7], data[8])
            v2 = (data[9], data[10], data[11])
            triangles.append(((v0, v1, v2), normal))
    return triangles


# ============================================================
# DA DISTRICT DETECTION + GML FILE ACCESS
# ============================================================

def find_da_districts_for_domain(x_min_ft, x_max_ft, y_min_ft, y_max_ft):
    """Find DA districts overlapping the domain bounding box."""
    matches = []
    for da, (dx_min, dx_max, dy_min, dy_max) in DA_ENVELOPES.items():
        if (x_max_ft >= dx_min and x_min_ft <= dx_max and
                y_max_ft >= dy_min and y_min_ft <= dy_max):
            matches.append(da)
    return sorted(matches)


def get_da_gml_path(da, citygml_dir):
    """Get path to DA GML file, extracting from zip or downloading if needed."""
    da_dir = os.path.join(citygml_dir, da)
    gml_path = os.path.join(da_dir, f"{da}_3D_Buildings_Merged.gml")
    zip_path = os.path.join(citygml_dir, f"{da}.gml.zip")

    if os.path.exists(gml_path):
        return gml_path

    if os.path.exists(zip_path):
        print(f"  Extracting {zip_path}...")
        os.makedirs(da_dir, exist_ok=True)
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(da_dir)
        if os.path.exists(gml_path):
            return gml_path

    url = GITHUB_RELEASE_URL.format(da=da)
    print(f"\n  DA file not found locally. Downloading {da}...")
    print(f"  URL: {url}")

    if HAS_URLLIB:
        try:
            os.makedirs(citygml_dir, exist_ok=True)
            print(f"  Downloading to {zip_path}...")
            urllib.request.urlretrieve(url, zip_path)
            print(f"  Extracting...")
            os.makedirs(da_dir, exist_ok=True)
            with zipfile.ZipFile(zip_path, 'r') as zf:
                zf.extractall(da_dir)
            if os.path.exists(gml_path):
                return gml_path
        except Exception as e:
            print(f"  Download failed: {e}")

    print(f"\n  MANUAL DOWNLOAD REQUIRED:")
    print(f"  1. Download: {url}")
    print(f"  2. Save to:  {zip_path}")
    print(f"  3. Re-run this script")
    return None


# ============================================================
# CITYGML STREAMING EXTRACTION
# ============================================================

def extract_buildings_in_box(gml_path, x_min, x_max, y_min, y_max):
    """Stream through GML file and extract buildings within bounding box.

    Uses a line-by-line state machine to avoid loading the entire file.
    """
    buildings = []
    current_building = []
    in_building = False
    building_in_box = False
    building_count = 0
    match_count = 0

    building_start = re.compile(r'<bldg:Building\s')
    building_end = '</bldg:Building>'
    poslist_pat = re.compile(r'<gml:posList>(.*?)</gml:posList>')

    file_size = os.path.getsize(gml_path)
    print(f"  Scanning {os.path.basename(gml_path)} ({file_size / 1024 / 1024:.0f} MB)")

    t0 = time.time()
    bytes_read = 0
    last_report = 0

    with open(gml_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            bytes_read += len(line.encode('utf-8', errors='replace'))

            if bytes_read - last_report > 50_000_000:
                elapsed = time.time() - t0
                pct = 100 * bytes_read / file_size if file_size > 0 else 0
                rate = bytes_read / 1024 / 1024 / elapsed if elapsed > 0 else 0
                print(f"    {pct:.0f}% ({bytes_read/1024/1024:.0f} MB, "
                      f"{rate:.0f} MB/s, {building_count} bldgs, "
                      f"{match_count} in domain)")
                last_report = bytes_read

            if not in_building:
                if building_start.search(line):
                    in_building = True
                    building_in_box = False
                    current_building = [line]
                    building_count += 1
            else:
                current_building.append(line)

                if not building_in_box:
                    m = poslist_pat.search(line)
                    if m:
                        coords = m.group(1).strip().split()
                        for i in range(0, len(coords) - 2, 3):
                            try:
                                x = float(coords[i])
                                y = float(coords[i + 1])
                                if x_min <= x <= x_max and y_min <= y <= y_max:
                                    building_in_box = True
                                    break
                            except ValueError:
                                continue

                if building_end in line:
                    if building_in_box:
                        buildings.append(''.join(current_building))
                        match_count += 1
                    in_building = False
                    current_building = []

    elapsed = time.time() - t0
    print(f"  Done: {building_count} scanned, "
          f"{match_count} in domain ({elapsed:.1f}s)")

    return buildings


# ============================================================
# CITYGML SURFACE PARSING
# ============================================================

def parse_building_surfaces(gml_text):
    """Parse Wall/Roof/Ground surfaces from a building's GML text."""
    surfaces = []
    pattern = re.compile(
        r'<bldg:(GroundSurface|RoofSurface|WallSurface)\s[^>]*>.*?'
        r'<gml:posList>(.*?)</gml:posList>.*?'
        r'</bldg:\1>',
        re.DOTALL
    )
    for match in pattern.finditer(gml_text):
        stype = match.group(1)
        coords = [float(x) for x in match.group(2).strip().split()]
        n = len(coords) // 3
        if n < 3:
            continue
        vertices = [(coords[i*3], coords[i*3+1], coords[i*3+2]) for i in range(n)]
        if len(vertices) > 1 and vertices[-1] == vertices[0]:
            vertices = vertices[:-1]
        if len(vertices) >= 3:
            surfaces.append({
                'type': stype,
                'vertices_ft': vertices,
                'max_z_ft': max(v[2] for v in vertices),
                'min_z_ft': min(v[2] for v in vertices),
            })
    return surfaces


def get_building_attr(gml_text, attr_name):
    """Extract a string attribute from building GML."""
    m = re.search(
        rf'<gen:stringAttribute name="{attr_name}">\s*<gen:value>(.*?)</gen:value>',
        gml_text
    )
    return m.group(1) if m else None


# ============================================================
# SOCRATA API QUERY (fallback)
# ============================================================

def build_socrata_bbox():
    """Compute WGS84 bounding box from domain bounds + 50m padding."""
    pad = 50
    nw_lat = CENTER_LAT + (DOMAIN_Y_MAX + pad) / M_PER_DEG_LAT
    nw_lon = CENTER_LON + (DOMAIN_X_MIN - pad) / M_PER_DEG_LON
    se_lat = CENTER_LAT + (DOMAIN_Y_MIN - pad) / M_PER_DEG_LAT
    se_lon = CENTER_LON + (DOMAIN_X_MAX + pad) / M_PER_DEG_LON
    return nw_lat, nw_lon, se_lat, se_lon


def fetch_socrata_buildings(year_filter=None):
    """Query Socrata API for buildings in the CFD domain."""
    nw_lat, nw_lon, se_lat, se_lon = build_socrata_bbox()
    where = f"within_box(the_geom, {nw_lat}, {nw_lon}, {se_lat}, {se_lon})"
    if year_filter:
        where += f" AND (cnstrct_yr IS NULL OR cnstrct_yr <= {year_filter})"
    params = {"$where": where, "$limit": 5000}
    print(f"  Querying Socrata API...")
    resp = requests.get(NYC_API_URL, params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()
    if data.get("type") != "FeatureCollection":
        raise ValueError(f"Unexpected response type: {data.get('type')}")
    features = data.get("features", [])
    print(f"  Received {len(features)} building footprints from API")
    return features


def parse_socrata_building(feature):
    """Parse one GeoJSON feature -> (verts_local_m, height_m, bin_str) or None."""
    props = feature.get("properties", {})
    geom = feature.get("geometry", {})

    height_ft = float(props.get("height_roof") or 0)
    height_m = height_ft * FT_TO_M
    if height_m < MIN_HEIGHT_M:
        return None

    bin_str = str(props.get("bin") or "")
    geom_type = geom.get("type", "")
    coords = geom.get("coordinates", [])

    # Extract the largest exterior ring
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

    if abs(polygon_area_2d(verts)) < MIN_AREA_M2:
        return None

    return (verts, height_m, bin_str)


def is_citicorp(bin_str, verts=None):
    """Identify Citicorp Center by BIN or proximity to origin."""
    if bin_str in CITICORP_BIN_SET:
        return True
    if verts:
        cx = sum(v[0] for v in verts) / len(verts)
        cy = sum(v[1] for v in verts) / len(verts)
        if math.sqrt(cx * cx + cy * cy) < CITICORP_PROXIMITY_M:
            return True
    return False


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Hybrid STL generator: CityGML primary + Socrata fallback',
        epilog='Examples:\n'
               '  python generate_stl_hybrid.py                  # Citicorp (default)\n'
               '  python generate_stl_hybrid.py --center-lat 40.7489 '
               '--center-lon -73.9680 --name un_hq --radius 300 --no-target\n',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Location parameters
    parser.add_argument('--center-lat', type=float, default=DEFAULT_CENTER_LAT,
                        help=f'Center latitude (default: {DEFAULT_CENTER_LAT} = Citicorp)')
    parser.add_argument('--center-lon', type=float, default=DEFAULT_CENTER_LON,
                        help=f'Center longitude (default: {DEFAULT_CENTER_LON} = Citicorp)')
    parser.add_argument('--name', type=str, default='citicorp',
                        help='Project name for output files (default: citicorp)')
    parser.add_argument('--radius', type=float, default=None,
                        help='Domain radius in meters (default: 360 for citicorp, symmetric)')
    parser.add_argument('--no-target', action='store_true',
                        help='Skip hand-crafted target building (surroundings only)')
    # Exclusion
    parser.add_argument('--exclude-bins', nargs='*', default=None,
                        help='BINs to exclude from surroundings (default: Citicorp BINs)')
    parser.add_argument('--exclude-radius', type=float, default=None,
                        help='Proximity exclusion radius in meters (default: 40 for Citicorp, 0 otherwise)')
    # Data & output
    parser.add_argument('--data-dir', type=str, default=None,
                        help='CityGML data directory (default: _citygml/ next to script)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='STL output directory (default: constant/triSurface/)')
    parser.add_argument('--offline', action='store_true',
                        help='Skip Socrata API, CityGML only')
    parser.add_argument('--year', type=int, default=None,
                        help='Filter to buildings built by this year (e.g., 1978)')
    parser.add_argument('--ascii', action='store_true',
                        help='Also write ASCII STL copies')
    parser.add_argument('--lod2-tower', action='store_true',
                        help='Use CityGML LoD2 tower (124 tris) instead of simplified box')

    args = parser.parse_args()

    # ---- Dynamic constants from CLI args ----
    global CENTER_LAT, CENTER_LON, CENTER_X_FT, CENTER_Y_FT
    global M_PER_DEG_LAT, M_PER_DEG_LON
    global DOMAIN_X_MIN, DOMAIN_X_MAX, DOMAIN_Y_MIN, DOMAIN_Y_MAX

    CENTER_LAT = args.center_lat
    CENTER_LON = args.center_lon
    CENTER_X_FT = 992506.0 + (CENTER_LON - DEFAULT_CENTER_LON) * FT_PER_DEG_LON
    CENTER_Y_FT = 215620.0 + (CENTER_LAT - DEFAULT_CENTER_LAT) * FT_PER_DEG_LAT
    M_PER_DEG_LAT = 111320.0
    M_PER_DEG_LON = M_PER_DEG_LAT * math.cos(math.radians(CENTER_LAT))

    # Domain sizing
    is_citicorp_mode = (args.name == 'citicorp' and not args.no_target)
    if args.radius is not None:
        r = args.radius
        DOMAIN_X_MIN = -r
        DOMAIN_X_MAX = r
        DOMAIN_Y_MIN = -r
        DOMAIN_Y_MAX = r
    elif is_citicorp_mode:
        # Citicorp default: asymmetric (more downstream)
        DOMAIN_X_MIN = -200
        DOMAIN_X_MAX = 520
        DOMAIN_Y_MIN = -360
        DOMAIN_Y_MAX = 360
    else:
        DOMAIN_X_MIN = -360
        DOMAIN_X_MAX = 360
        DOMAIN_Y_MIN = -360
        DOMAIN_Y_MAX = 360

    # Exclusion logic
    if args.exclude_bins is not None:
        exclude_bins = set(args.exclude_bins)
    elif is_citicorp_mode:
        exclude_bins = CITICORP_BIN_SET
    else:
        exclude_bins = set()

    if args.exclude_radius is not None:
        exclude_proximity = args.exclude_radius
    elif is_citicorp_mode:
        exclude_proximity = CITICORP_PROXIMITY_M
    else:
        exclude_proximity = 0.0

    script_dir = os.path.dirname(os.path.abspath(__file__))
    citygml_dir = args.data_dir or os.path.join(script_dir, '_citygml')
    out_dir = args.output_dir or os.path.join(script_dir, 'constant', 'triSurface')
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 65)
    print(f"Hybrid STL Generator -- {args.name}")
    print("=" * 65)
    print(f"  Center: {CENTER_LAT:.4f}N, {CENTER_LON:.4f}W")
    print(f"  Domain: X=[{DOMAIN_X_MIN}, {DOMAIN_X_MAX}], Y=[{DOMAIN_Y_MIN}, {DOMAIN_Y_MAX}] m")
    if args.no_target:
        print(f"  Mode: surroundings only (no target building)")
    if exclude_bins:
        print(f"  Excluding BINs: {', '.join(sorted(exclude_bins))}")
    if exclude_proximity > 0:
        print(f"  Excluding buildings within {exclude_proximity}m of center")
    if args.year:
        print(f"  Historical filter: buildings built by {args.year}")
    if args.offline:
        print(f"  Offline mode: Socrata API disabled")
    print()

    total_tris = 0
    tower_tris = []
    stilt_tris = []

    # ------------------------------------------------------------------
    # STEP 1: Target building (Citicorp hand-crafted, or skip)
    # ------------------------------------------------------------------
    if args.no_target:
        print("--- Step 1: Target building (skipped -- --no-target) ---")
    else:
        print(f"--- Step 1: Citicorp Tower + Stilts (rotated {TOWER_ANGLE_DEG} deg) ---")

        # Option: use CityGML LoD2 tower (already in correct coordinate frame)
        lod2_path = os.path.join(out_dir, 'citicorp_tower_lod2_bin.stl')
        if args.lod2_tower and os.path.exists(lod2_path):
            tower_tris = read_binary_stl(lod2_path)
            print(f"  Using CityGML LoD2 tower: {lod2_path}")
        else:
            if args.lod2_tower:
                print(f"  WARNING: LoD2 tower not found at {lod2_path}, using simplified box")
            # Generate simplified box in local frame, then rotate to match grid
            tower_tris_local = box_triangles(
                -HALF_T, -HALF_T, STILT_H,
                 HALF_T,  HALF_T, TOWER_H_TOP
            )
            tower_tris = rotate_triangles(tower_tris_local, TOWER_ANGLE_DEG)
        n = write_binary_stl(os.path.join(out_dir, 'citicorp_tower.stl'),
                             'citicorp_tower', tower_tris)
        total_tris += n
        print(f"  citicorp_tower.stl: {n} triangles")

        if args.ascii:
            write_ascii_stl(os.path.join(out_dir, 'citicorp_tower_ascii.stl'),
                            'citicorp_tower', tower_tris)

        stilt_tris_local = []
        # 4 face-midpoint stilts (rectangular: 7.27m face × 6.55m depth)
        hf = STILT_FACE / 2   # half-width along face
        hd = STILT_DEPTH / 2  # half-width toward center
        for cx, cy, axis in STILTS:
            if axis == 'x':
                # South/North: wider in X (face direction)
                stilt_tris_local.extend(box_triangles(
                    cx - hf, cy - hd, 0,
                    cx + hf, cy + hd, STILT_H
                ))
            else:
                # East/West: wider in Y (face direction)
                stilt_tris_local.extend(box_triangles(
                    cx - hd, cy - hf, 0,
                    cx + hd, cy + hf, STILT_H
                ))
        # Center core (elevator/mechanical, octagonal cross-section)
        core_verts = make_chamfered_rect(0, 0, CORE_W, CORE_D, CORE_CHAMFER)
        core_tris = extrude_polygon_to_triangles(core_verts, 0.0, STILT_H)
        stilt_tris_local.extend(core_tris)
        # Rotate entire stilt assembly to match tower
        stilt_tris = rotate_triangles(stilt_tris_local, TOWER_ANGLE_DEG)
        n = write_binary_stl(os.path.join(out_dir, 'citicorp_stilts.stl'),
                             'citicorp_stilts', stilt_tris)
        total_tris += n
        print(f"  citicorp_stilts.stl: {n} triangles "
              f"(4 stilts {STILT_FACE:.1f}m×{STILT_DEPTH:.1f}m + "
              f"core {CORE_W:.0f}m×{CORE_D:.0f}m)")

        if args.ascii:
            write_ascii_stl(os.path.join(out_dir, 'citicorp_stilts_ascii.stl'),
                            'citicorp_stilts', stilt_tris)

    # ------------------------------------------------------------------
    # STEP 2: CityGML surroundings (primary source)
    # ------------------------------------------------------------------
    print("\n--- Step 2: CityGML Surroundings (primary) ---")

    surr_tris = []
    citygml_bins = set()
    citygml_centroids = []   # list of (cx, cy) for spatial dedup
    citygml_count = 0
    citygml_surface_counts = {'WallSurface': 0, 'RoofSurface': 0, 'GroundSurface': 0}
    citygml_ok = False

    # Convert CFD domain corners to EPSG:2263 feet
    x_min_ft = CENTER_X_FT + DOMAIN_X_MIN * M_TO_US_SURVEY_FT
    x_max_ft = CENTER_X_FT + DOMAIN_X_MAX * M_TO_US_SURVEY_FT
    y_min_ft = CENTER_Y_FT + DOMAIN_Y_MIN * M_TO_US_SURVEY_FT
    y_max_ft = CENTER_Y_FT + DOMAIN_Y_MAX * M_TO_US_SURVEY_FT

    print(f"  Domain (EPSG:2263): X=[{x_min_ft:.0f}, {x_max_ft:.0f}], "
          f"Y=[{y_min_ft:.0f}, {y_max_ft:.0f}]")

    districts = find_da_districts_for_domain(x_min_ft, x_max_ft, y_min_ft, y_max_ft)
    if not districts:
        print("  WARNING: No DA districts overlap domain. Skipping CityGML.")
    else:
        print(f"  DA districts to search: {', '.join(districts)}")

        # Extract buildings from each DA district
        all_buildings_gml = []
        for da in districts:
            print(f"\n  --- {da} ---")
            gml_path = get_da_gml_path(da, citygml_dir)
            if gml_path is None:
                print(f"  Skipping {da} (file not available)")
                continue
            buildings = extract_buildings_in_box(gml_path, x_min_ft, x_max_ft,
                                                 y_min_ft, y_max_ft)
            all_buildings_gml.extend(buildings)

        if all_buildings_gml:
            # Find ground elevation from GroundSurface elements
            z_grounds = []
            for bldg_text in all_buildings_gml:
                surfs = parse_building_surfaces(bldg_text)
                for s in surfs:
                    if s['type'] == 'GroundSurface':
                        z_grounds.extend(v[2] for v in s['vertices_ft'])
            z_ground_ft = min(z_grounds) if z_grounds else 25.78
            print(f"\n  Ground elevation: {z_ground_ft:.2f} ft MSL "
                  f"({(z_ground_ft * US_SURVEY_FT_TO_M):.1f} m)")

            # Process each building
            seen_bins = set()
            for bldg_text in all_buildings_gml:
                bin_val = get_building_attr(bldg_text, 'BIN') or ''

                # Deduplicate (building at DA boundary may appear in multiple files)
                if bin_val and bin_val in seen_bins:
                    continue
                if bin_val:
                    seen_bins.add(bin_val)

                # Skip excluded BINs (Citicorp in default mode)
                if bin_val and bin_val in exclude_bins:
                    continue

                # Year filter
                if args.year:
                    year_str = get_building_attr(bldg_text, 'yearbuilt') or '0'
                    try:
                        year_val = int(year_str)
                        if year_val > 0 and year_val > args.year:
                            continue
                    except ValueError:
                        pass

                # Parse surfaces (done once, used for proximity + triangles)
                surfs = parse_building_surfaces(bldg_text)
                if not surfs:
                    continue

                # Skip by proximity to center (catches BINs not in exclude set)
                if exclude_proximity > 0:
                    all_surf_xs = []
                    all_surf_ys = []
                    for s in surfs:
                        for x, y, z in s['vertices_ft']:
                            all_surf_xs.append((x - CENTER_X_FT) * US_SURVEY_FT_TO_M)
                            all_surf_ys.append((y - CENTER_Y_FT) * US_SURVEY_FT_TO_M)
                    if all_surf_xs:
                        bcx = sum(all_surf_xs) / len(all_surf_xs)
                        bcy = sum(all_surf_ys) / len(all_surf_ys)
                        if math.sqrt(bcx*bcx + bcy*bcy) < exclude_proximity:
                            continue

                bldg_tris = []
                for surf in surfs:
                    citygml_surface_counts[surf['type']] = \
                        citygml_surface_counts.get(surf['type'], 0) + 1
                    verts_m = [citygml_to_local(x, y, z, z_ground_ft)
                               for x, y, z in surf['vertices_ft']]
                    tris = fan_triangulate(verts_m)
                    for v0, v1, v2 in tris:
                        normal = compute_normal(v0, v1, v2)
                        bldg_tris.append(((v0, v1, v2), normal))

                if bldg_tris:
                    surr_tris.extend(bldg_tris)
                    if bin_val:
                        citygml_bins.add(bin_val)
                    # Store centroid for spatial dedup
                    all_xs = []
                    all_ys = []
                    for (v0, v1, v2), _ in bldg_tris:
                        for v in (v0, v1, v2):
                            all_xs.append(v[0])
                            all_ys.append(v[1])
                    if all_xs:
                        citygml_centroids.append(
                            (sum(all_xs)/len(all_xs), sum(all_ys)/len(all_ys)))
                    citygml_count += 1

            citygml_ok = True

    print(f"\n  CityGML result: {citygml_count} buildings, "
          f"{len(surr_tris)} triangles")
    for stype, count in sorted(citygml_surface_counts.items()):
        if count > 0:
            print(f"    {stype}: {count}")

    # ------------------------------------------------------------------
    # STEP 3: Socrata fallback (for post-2014 buildings)
    # ------------------------------------------------------------------
    socrata_count = 0
    socrata_tris_before = len(surr_tris)

    if args.offline:
        print("\n--- Step 3: Socrata Fallback (skipped -- offline mode) ---")
    elif not HAS_REQUESTS:
        print("\n--- Step 3: Socrata Fallback (skipped -- requests not installed) ---")
    else:
        print("\n--- Step 3: Socrata Fallback (post-2014 buildings) ---")
        try:
            features = fetch_socrata_buildings(year_filter=args.year)
            skipped_citygml = 0
            skipped_citicorp = 0
            skipped_small = 0

            DEDUP_RADIUS = 15.0  # meters — skip Socrata if centroid within this of CityGML
            skipped_proximity = 0

            for feature in features:
                parsed = parse_socrata_building(feature)
                if parsed is None:
                    skipped_small += 1
                    continue

                verts, height_m, bin_str = parsed

                # Skip excluded BINs or proximity to center
                if bin_str and bin_str in exclude_bins:
                    skipped_citicorp += 1
                    continue
                if exclude_proximity > 0:
                    cx = sum(v[0] for v in verts) / len(verts)
                    cy = sum(v[1] for v in verts) / len(verts)
                    if math.sqrt(cx*cx + cy*cy) < exclude_proximity:
                        skipped_citicorp += 1
                        continue

                if bin_str and bin_str in citygml_bins:
                    skipped_citygml += 1
                    continue

                # Spatial proximity check: skip if centroid near any CityGML building
                if citygml_centroids:
                    cx = sum(v[0] for v in verts) / len(verts)
                    cy = sum(v[1] for v in verts) / len(verts)
                    too_close = False
                    for gcx, gcy in citygml_centroids:
                        if abs(cx - gcx) < DEDUP_RADIUS and abs(cy - gcy) < DEDUP_RADIUS:
                            dist = math.sqrt((cx - gcx)**2 + (cy - gcy)**2)
                            if dist < DEDUP_RADIUS:
                                too_close = True
                                break
                    if too_close:
                        skipped_proximity += 1
                        continue

                tris = extrude_polygon_to_triangles(verts, 0.0, height_m)
                surr_tris.extend(tris)
                socrata_count += 1

            socrata_tri_count = len(surr_tris) - socrata_tris_before
            print(f"  Socrata result: {socrata_count} new buildings, "
                  f"{socrata_tri_count} triangles")
            print(f"    Skipped: {skipped_citygml} BIN match, "
                  f"{skipped_proximity} proximity (<{DEDUP_RADIUS}m), "
                  f"{skipped_citicorp} excluded (target), "
                  f"{skipped_small} below {MIN_HEIGHT_M}m")

        except Exception as e:
            print(f"  Socrata API failed: {e}")
            if not citygml_ok:
                print("  Both sources failed. Using hardcoded fallback.")
                for name, cx, cy, wx, wy, h in FALLBACK_SURROUNDINGS:
                    surr_tris.extend(box_triangles(
                        cx - wx/2, cy - wy/2, 0,
                        cx + wx/2, cy + wy/2, h))
                citygml_count = 0
                socrata_count = len(FALLBACK_SURROUNDINGS)

    # ------------------------------------------------------------------
    # STEP 4: Write surroundings STL
    # ------------------------------------------------------------------
    if is_citicorp_mode:
        surr_filename = 'surroundings.stl'
    else:
        surr_filename = f'{args.name}_surroundings.stl'

    print(f"\n--- Step 4: Write {surr_filename} ---")

    n = write_binary_stl(os.path.join(out_dir, surr_filename),
                         f'{args.name}_surroundings', surr_tris)
    total_tris += n
    print(f"  {surr_filename}: {n} triangles")

    if args.ascii:
        surr_ascii = surr_filename.replace('.stl', '_ascii.stl')
        write_ascii_stl(os.path.join(out_dir, surr_ascii),
                        f'{args.name}_surroundings', surr_tris)

    # Bounding box
    if surr_tris:
        all_verts = []
        for verts_tuple, _ in surr_tris:
            all_verts.extend(verts_tuple)
        xs = [v[0] for v in all_verts]
        ys = [v[1] for v in all_verts]
        zs = [v[2] for v in all_verts]
        print(f"  Bounding box:")
        print(f"    X: {min(xs):.1f} to {max(xs):.1f} m")
        print(f"    Y: {min(ys):.1f} to {max(ys):.1f} m")
        print(f"    Z: {min(zs):.1f} to {max(zs):.1f} m")

    # ------------------------------------------------------------------
    # PROVENANCE REPORT
    # ------------------------------------------------------------------
    total_buildings = citygml_count + socrata_count
    print("\n" + "=" * 65)
    print(f"PROVENANCE REPORT — {args.name}")
    print("=" * 65)
    print(f"  Center:                {CENTER_LAT:.4f}N, {CENTER_LON:.4f}W")
    print(f"  CityGML (primary):     {citygml_count} buildings "
          f"({citygml_count*100//max(total_buildings,1)}%)")
    print(f"  Socrata (fallback):    {socrata_count} buildings "
          f"({socrata_count*100//max(total_buildings,1)}%)")
    print(f"  Total surroundings:    {total_buildings} buildings, "
          f"{len(surr_tris)} triangles")
    if exclude_bins:
        print(f"  Excluded BINs:         {', '.join(sorted(exclude_bins))}")
    if exclude_proximity > 0:
        print(f"  Excluded proximity:    <{exclude_proximity}m from center")
    if args.year:
        print(f"  Historical filter:     buildings built by {args.year}")
    else:
        print(f"  Historical filter:     None (present-day)")
    print()
    print(f"  Output files in {out_dir}/:")
    if tower_tris:
        print(f"    citicorp_tower.stl:  {len(tower_tris)} triangles")
    if stilt_tris:
        print(f"    citicorp_stilts.stl: {len(stilt_tris)} triangles "
              f"(4 stilts + center core)")
    print(f"    {surr_filename}:{' ' * max(1, 20-len(surr_filename))}{len(surr_tris)} triangles")
    print(f"    TOTAL:               {total_tris} triangles")
    print("=" * 65)


if __name__ == '__main__':
    main()
