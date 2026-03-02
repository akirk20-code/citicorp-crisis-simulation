#!/usr/bin/env python3
"""Generate STL geometry for Citicorp Center CFD simulation.

Fetches real building footprints from the NYC Open Data API (Socrata)
and extrudes them to surveyed roof heights. Citicorp tower and stilts
use hand-crafted detailed geometry.

Dependencies: requests (pip install requests)
Fallback: If the API is unreachable or requests is not installed,
          uses hardcoded approximate geometry for 18 surrounding buildings.

All dimensions in meters. Citicorp Center at origin (0, 0, 0).
"""

import os
import sys
import math
import struct
import json

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# ============================================================
# CITICORP GEOMETRY (meters)
# ============================================================
TOWER_W = 47.85       # 157 ft — square plan
TOWER_H_TOP = 278.9   # 915 ft — roof height
STILT_H = 34.75       # 114 ft — stilt height (tower bottom)
STILT_W = 7.32        # 24 ft — stilt cross-section
HALF_T = TOWER_W / 2  # 23.925 m
HALF_S = STILT_W / 2  # 3.66 m

# Stilt centers (at face midpoints, NOT corners)
STILTS = [
    (0,        -HALF_T),   # South
    (HALF_T,    0),        # East
    (0,         HALF_T),   # North
    (-HALF_T,   0),        # West
]

# ============================================================
# NYC OPEN DATA API PARAMETERS
# ============================================================
CITICORP_LAT = 40.7579
CITICORP_LON = -73.9690

API_URL = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"

# Domain bounds (meters relative to Citicorp center)
DOMAIN_X_MIN = -200
DOMAIN_X_MAX = 520
DOMAIN_Y_MIN = -360
DOMAIN_Y_MAX = 360

# Projection constants
M_PER_DEG_LAT = 111320.0
M_PER_DEG_LON = M_PER_DEG_LAT * math.cos(math.radians(CITICORP_LAT))

FT_TO_M = 0.3048
MIN_HEIGHT_M = 5.0        # skip buildings shorter than this
CITICORP_PROXIMITY_M = 35  # radius to identify Citicorp by location

# ============================================================
# COORDINATE TRANSFORM
# ============================================================

def wgs84_to_local(lon, lat):
    """Convert WGS84 (lon, lat) to local meters centered on Citicorp."""
    x = (lon - CITICORP_LON) * M_PER_DEG_LON
    y = (lat - CITICORP_LAT) * M_PER_DEG_LAT
    return x, y


# ============================================================
# POLYGON TRIANGULATION (ear-clipping)
# ============================================================

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


def triangulate_polygon(verts):
    """Ear-clipping triangulation of a simple 2D polygon.

    Returns list of (i, j, k) index triples into the original vertex list.
    Assumes CCW winding.
    """
    n = len(verts)
    if n < 3:
        return []
    if n == 3:
        return [(0, 1, 2)]

    # Work with a copy of indices
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

            # Is triangle (p, c, nx) a valid ear?
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
                if point_in_triangle_2d(verts[idx[j]],
                                        verts[p], verts[c], verts[nx]):
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


# ============================================================
# STL PRIMITIVES
# ============================================================

def _normalize(v):
    length = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if length < 1e-15:
        return (0.0, 0.0, 1.0)
    return (v[0] / length, v[1] / length, v[2] / length)


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


def extrude_polygon_to_triangles(verts_2d, z_base, z_top):
    """Extrude a CCW 2D polygon into a watertight 3D solid.

    Returns list of ((v0, v1, v2), normal) triangles.
    """
    n = len(verts_2d)
    if n < 3 or z_top <= z_base:
        return []

    tri_idx = triangulate_polygon(verts_2d)
    result = []

    # --- Top face (normal +z) ---
    for i0, i1, i2 in tri_idx:
        v0 = (verts_2d[i0][0], verts_2d[i0][1], z_top)
        v1 = (verts_2d[i1][0], verts_2d[i1][1], z_top)
        v2 = (verts_2d[i2][0], verts_2d[i2][1], z_top)
        result.append(((v0, v1, v2), (0.0, 0.0, 1.0)))

    # --- Bottom face (normal -z, reverse winding) ---
    for i0, i1, i2 in tri_idx:
        v0 = (verts_2d[i0][0], verts_2d[i0][1], z_base)
        v1 = (verts_2d[i1][0], verts_2d[i1][1], z_base)
        v2 = (verts_2d[i2][0], verts_2d[i2][1], z_base)
        result.append(((v0, v2, v1), (0.0, 0.0, -1.0)))

    # --- Side walls ---
    for i in range(n):
        j = (i + 1) % n
        bl = (verts_2d[i][0], verts_2d[i][1], z_base)
        br = (verts_2d[j][0], verts_2d[j][1], z_base)
        tl = (verts_2d[i][0], verts_2d[i][1], z_top)
        tr = (verts_2d[j][0], verts_2d[j][1], z_top)

        # Outward normal for CCW polygon: (dy, -dx, 0)
        dx = verts_2d[j][0] - verts_2d[i][0]
        dy = verts_2d[j][1] - verts_2d[i][1]
        wall_n = _normalize((dy, -dx, 0.0))

        result.append(((bl, br, tr), wall_n))
        result.append(((bl, tr, tl), wall_n))

    return result


# ============================================================
# STL FILE WRITERS
# ============================================================

def write_stl_ascii(filepath, solid_name, triangles):
    """Write ASCII STL. Returns triangle count."""
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


# ============================================================
# NYC OPEN DATA QUERY
# ============================================================

def _build_bbox():
    """Compute WGS84 bounding box from domain bounds + 50 m padding."""
    pad = 50  # meters
    nw_lat = CITICORP_LAT + (DOMAIN_Y_MAX + pad) / M_PER_DEG_LAT
    nw_lon = CITICORP_LON + (DOMAIN_X_MIN - pad) / M_PER_DEG_LON
    se_lat = CITICORP_LAT + (DOMAIN_Y_MIN - pad) / M_PER_DEG_LAT
    se_lon = CITICORP_LON + (DOMAIN_X_MAX + pad) / M_PER_DEG_LON
    return nw_lat, nw_lon, se_lat, se_lon


def fetch_nyc_buildings():
    """Query NYC Open Data for building footprints inside the CFD domain."""
    nw_lat, nw_lon, se_lat, se_lon = _build_bbox()
    params = {
        "$where": f"within_box(the_geom, {nw_lat}, {nw_lon}, {se_lat}, {se_lon})",
        "$limit": 5000,
    }
    print(f"  Querying NYC Open Data API...")
    print(f"  Bounding box: NW ({nw_lat:.5f}, {nw_lon:.5f})  "
          f"SE ({se_lat:.5f}, {se_lon:.5f})")

    resp = requests.get(API_URL, params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    if data.get("type") != "FeatureCollection":
        raise ValueError(f"Unexpected response type: {data.get('type')}")

    features = data.get("features", [])
    print(f"  Received {len(features)} building footprints from API")
    return features


def parse_building(feature):
    """Parse one GeoJSON feature → (verts_local_m, height_m, bin_str) or None."""
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

    if ring is None or len(ring) < 4:  # need at least 3 unique + closing
        return None

    # Convert to local meters
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

    return (verts, height_m, bin_str)


def is_citicorp(bin_str, verts):
    """Identify Citicorp Center by BIN or proximity to origin."""
    if bin_str in ("1035879", "1087931"):
        return True
    if verts:
        cx = sum(v[0] for v in verts) / len(verts)
        cy = sum(v[1] for v in verts) / len(verts)
        if math.sqrt(cx * cx + cy * cy) < CITICORP_PROXIMITY_M:
            return True
    return False


# ============================================================
# FALLBACK HARDCODED GEOMETRY
# ============================================================

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


def make_fallback_surroundings():
    tris = []
    for name, cx, cy, wx, wy, h in FALLBACK_SURROUNDINGS:
        tris.extend(box_triangles(cx - wx/2, cy - wy/2, 0,
                                  cx + wx/2, cy + wy/2, h))
    return tris, len(FALLBACK_SURROUNDINGS)


# ============================================================
# MAIN
# ============================================================

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(script_dir, "constant", "triSurface")
    os.makedirs(out_dir, exist_ok=True)

    offline = "--offline" in sys.argv
    total_tris = 0

    # --- Citicorp Tower (always use detailed hand-crafted geometry) ---
    print("\n--- Citicorp Tower ---")
    tower_tris = box_triangles(
        -HALF_T, -HALF_T, STILT_H,
         HALF_T,  HALF_T, TOWER_H_TOP
    )
    n = write_stl_ascii(os.path.join(out_dir, "citicorp_tower.stl"),
                        "citicorp_tower", tower_tris)
    total_tris += n
    print(f"  citicorp_tower.stl: {n} triangles")
    print(f"    Bounds: ({-HALF_T:.1f}, {-HALF_T:.1f}, {STILT_H:.1f}) "
          f"to ({HALF_T:.1f}, {HALF_T:.1f}, {TOWER_H_TOP:.1f})")

    # --- Stilt Columns ---
    stilt_tris = []
    for cx, cy in STILTS:
        stilt_tris.extend(box_triangles(
            cx - HALF_S, cy - HALF_S, 0,
            cx + HALF_S, cy + HALF_S, STILT_H
        ))
    n = write_stl_ascii(os.path.join(out_dir, "citicorp_stilts.stl"),
                        "citicorp_stilts", stilt_tris)
    total_tris += n
    print(f"  citicorp_stilts.stl: {n} triangles (4 stilts)")

    # --- Surrounding Buildings ---
    print("\n--- Surrounding Buildings ---")
    surr_tris = []
    building_count = 0
    used_api = False

    if HAS_REQUESTS and not offline:
        try:
            features = fetch_nyc_buildings()
            citicorp_found = False
            skipped = 0
            tallest = ("", 0)

            for feature in features:
                parsed = parse_building(feature)
                if parsed is None:
                    skipped += 1
                    continue

                verts, height_m, bin_str = parsed

                if is_citicorp(bin_str, verts):
                    citicorp_found = True
                    print(f"  Identified Citicorp (BIN={bin_str}), "
                          f"using hand-crafted model instead")
                    continue

                tris = extrude_polygon_to_triangles(verts, 0.0, height_m)
                surr_tris.extend(tris)
                building_count += 1

                if height_m > tallest[1]:
                    tallest = (bin_str, height_m)

            used_api = True
            if not citicorp_found:
                print("  Warning: Citicorp not auto-identified in API data")
            print(f"  Processed {building_count} buildings "
                  f"(skipped {skipped} below {MIN_HEIGHT_M}m)")
            if tallest[1] > 0:
                print(f"  Tallest neighbor: BIN {tallest[0]}, "
                      f"{tallest[1]:.1f} m ({tallest[1]/FT_TO_M:.0f} ft)")

        except Exception as e:
            print(f"  API query failed: {e}")
            print("  Falling back to hardcoded geometry...")
            surr_tris, building_count = make_fallback_surroundings()
    else:
        reason = "--offline flag" if offline else "requests not installed"
        print(f"  Using hardcoded geometry ({reason})")
        surr_tris, building_count = make_fallback_surroundings()

    n = write_stl_ascii(os.path.join(out_dir, "surroundings.stl"),
                        "surroundings", surr_tris)
    total_tris += n
    source = "NYC Open Data API" if used_api else "hardcoded fallback"
    print(f"  surroundings.stl: {n} triangles "
          f"({building_count} buildings, source: {source})")

    # --- Summary ---
    print(f"\nTotal: {total_tris} triangles in {out_dir}/")
    print("STL generation complete.")


if __name__ == "__main__":
    main()
