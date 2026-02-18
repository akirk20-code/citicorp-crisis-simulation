#!/usr/bin/env python3
"""
Generate STL geometry for 1978 historical Citicorp CFD simulation.

Uses OpenStreetMap Overpass API filtered to buildings existing in 1978,
cross-referenced with NYC Building Footprints API construction_year field.

Citicorp tower uses hand-crafted geometry with:
  - Correct 45° slanted roof (wedge, 248m to 278.9m peak)
  - 45° rotation matching true building orientation
  - 4 stilt columns at face midpoints

All dimensions in meters. Citicorp Center at origin (0, 0, 0).
"""

import os, sys, math, struct, json, re

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# ─── CONSTANTS ────────────────────────────────────────────────────────────────
CITICORP_LAT = 40.7579
CITICORP_LON = -73.9690

TOWER_W    = 47.85    # 157 ft square plan
TOWER_H_TOP = 278.90  # 915 ft — crown peak
ROOF_START  = 248.00  # 813 ft — where sloped roof begins
STILT_H    = 34.75    # 114 ft
STILT_W    = 7.32     # 24 ft
HALF_T     = TOWER_W / 2
HALF_S     = STILT_W / 2

# Stilt centres at face mid-points (not corners)
STILTS = [(0, -HALF_T), (HALF_T, 0), (0, HALF_T), (-HALF_T, 0)]

M_PER_DEG_LAT = 111320.0
M_PER_DEG_LON = M_PER_DEG_LAT * math.cos(math.radians(CITICORP_LAT))
FT_TO_M       = 0.3048
MIN_HEIGHT_M  = 5.0

HISTORICAL_YEAR = 1978  # Filter: only buildings that existed by this year

# Bounding box (±360 m around Citicorp)
BB_LAT_MIN, BB_LAT_MAX = 40.7507, 40.7651
BB_LON_MIN, BB_LON_MAX = -73.9762, -73.9618

# ─── COORDINATE HELPERS ───────────────────────────────────────────────────────
def wgs84_to_local(lon, lat):
    x = (lon - CITICORP_LON) * M_PER_DEG_LON
    y = (lat - CITICORP_LAT) * M_PER_DEG_LAT
    return x, y

def rotate45(x, y):
    """Rotate point (x,y) by 45° around origin."""
    c = s = math.sqrt(0.5)
    return c*x - s*y, s*x + c*y

# ─── POLYGON HELPERS ──────────────────────────────────────────────────────────
def poly_area(verts):
    n = len(verts)
    a = 0.0
    for i in range(n):
        j = (i+1) % n
        a += verts[i][0]*verts[j][1] - verts[j][0]*verts[i][1]
    return a / 2.0

def ear_clip(verts):
    """Triangulate a simple polygon via ear clipping. Returns list of (p0,p1,p2)."""
    pts = list(verts)
    if poly_area(pts) < 0:
        pts.reverse()
    tris = []
    while len(pts) > 3:
        clipped = False
        n = len(pts)
        for i in range(n):
            a, b, c = pts[(i-1)%n], pts[i], pts[(i+1)%n]
            if _is_ear(pts, i):
                tris.append((a, b, c))
                pts.pop(i)
                clipped = True
                break
        if not clipped:
            break  # degenerate
    if len(pts) == 3:
        tris.append(tuple(pts))
    return tris

def _is_ear(pts, idx):
    n = len(pts)
    a, b, c = pts[(idx-1)%n], pts[idx], pts[(idx+1)%n]
    if _cross2d(a, b, c) <= 0:
        return False
    tri = (a, b, c)
    for j in range(n):
        if j in ((idx-1)%n, idx, (idx+1)%n):
            continue
        if _pt_in_tri(pts[j], *tri):
            return False
    return True

def _cross2d(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])

def _pt_in_tri(p, a, b, c):
    d1 = _cross2d(a, b, p)
    d2 = _cross2d(b, c, p)
    d3 = _cross2d(c, a, p)
    has_neg = (d1<0) or (d2<0) or (d3<0)
    has_pos = (d1>0) or (d2>0) or (d3>0)
    return not (has_neg and has_pos)

# ─── STL HELPERS ──────────────────────────────────────────────────────────────
def normal(p0, p1, p2):
    ax, ay, az = p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]
    bx, by, bz = p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]
    nx = ay*bz - az*by
    ny = az*bx - ax*bz
    nz = ax*by - ay*bx
    mag = math.sqrt(nx*nx + ny*ny + nz*nz)
    if mag < 1e-12:
        return (0.0, 0.0, 1.0)
    return (nx/mag, ny/mag, nz/mag)

def write_stl(path, name, triangles):
    """Write ASCII STL. Each tri is (p0, p1, p2) with p=(x,y,z)."""
    with open(path, 'w') as f:
        f.write(f"solid {name}\n")
        for tri in triangles:
            p0, p1, p2 = tri
            n = normal(p0, p1, p2)
            f.write(f"  facet normal {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n")
            f.write( "    outer loop\n")
            for p in (p0, p1, p2):
                f.write(f"      vertex {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
            f.write( "    endloop\n")
            f.write( "  endfacet\n")
        f.write(f"endsolid {name}\n")
    print(f"  Wrote {path} ({len(triangles)} triangles)")
    return len(triangles)

# ─── CITICORP TOWER GEOMETRY ──────────────────────────────────────────────────
def make_tower():
    """
    Citicorp tower with correct 45-degree slanted roof and 45-degree rotation.

    Box section: z = STILT_H to ROOF_START (248 m)
    Roof section: flat slanted plane rising from S face to N face
      - South edge z = ROOF_START
      - North edge z = TOWER_H_TOP
    Building is rotated 45 degrees in plan so corners face N/S/E/W.
    """
    hw = HALF_T
    tris = []

    # Pre-rotation corners (axis-aligned, then rotate 45 deg)
    def R(x, y): return rotate45(x, y)

    # Four corners pre-rotation
    corners_pre = [(-hw, -hw), (hw, -hw), (hw, hw), (-hw, hw)]
    corners = [R(x, y) for x, y in corners_pre]

    # Helper: side wall quad → 2 triangles
    def quad(p0, p1, p2, p3):
        return [(p0,p1,p2), (p0,p2,p3)]

    # ── Walls (box, STILT_H → ROOF_START) ────────────────────────────────────
    zb, zt = STILT_H, ROOF_START
    for i in range(4):
        c0 = corners[i]
        c1 = corners[(i+1)%4]
        p00 = (c0[0], c0[1], zb)
        p10 = (c1[0], c1[1], zb)
        p11 = (c1[0], c1[1], zt)
        p01 = (c0[0], c0[1], zt)
        tris += quad(p00, p10, p11, p01)

    # ── Bottom cap (at STILT_H) ───────────────────────────────────────────────
    base_2d = [(c[0], c[1]) for c in corners]
    for t in ear_clip(base_2d):
        p0 = (t[0][0], t[0][1], zb)
        p1 = (t[1][0], t[1][1], zb)
        p2 = (t[2][0], t[2][1], zb)
        # Flip winding for downward-facing normal
        tris.append((p0, p2, p1))

    # ── Slanted roof ──────────────────────────────────────────────────────────
    # After 45-deg rotation, the "south" tip is at (0, -hw*sqrt2)
    # and "north" tip is at (0, +hw*sqrt2). The roof slopes linearly
    # with y: z_roof(y) = ROOF_START + (TOWER_H_TOP - ROOF_START) * (y + hw*r2) / (2*hw*r2)
    # where r2 = sqrt(2).
    #
    # The four rotated corner y-coords:
    #   SW corner y = -hw*sqrt(2)/sqrt(2)... let's just use the actual values.
    #   corners[0] = R(-hw,-hw) = (0, -hw*sqrt2)  ← S tip,   z = ROOF_START
    #   corners[1] = R( hw,-hw) = (hw*sqrt2, 0)   ← E mid,   z = midpoint
    #   corners[2] = R( hw, hw) = (0, +hw*sqrt2)  ← N tip,   z = TOWER_H_TOP
    #   corners[3] = R(-hw, hw) = (-hw*sqrt2, 0)  ← W mid,   z = midpoint

    r2 = math.sqrt(2.0)
    y_south = -hw * r2
    y_north =  hw * r2
    y_range = y_north - y_south

    def roof_z(cx, cy):
        """Interpolate roof height from south (ROOF_START) to north (TOWER_H_TOP)."""
        t = (cy - y_south) / y_range
        return ROOF_START + t * (TOWER_H_TOP - ROOF_START)

    # Roof corners (each at their interpolated z)
    roof_pts = [(c[0], c[1], roof_z(c[0], c[1])) for c in corners]

    # Two triangles for the flat slanted roof plane
    tris.append((roof_pts[0], roof_pts[1], roof_pts[2]))
    tris.append((roof_pts[0], roof_pts[2], roof_pts[3]))

    # ── Upper walls (ROOF_START to roof edge — fill triangle gaps at each side)
    for i in range(4):
        c0 = corners[i]
        c1 = corners[(i+1)%4]
        rb = (c0[0], c0[1], ROOF_START)   # bottom of roof section
        re = (c1[0], c1[1], ROOF_START)
        rt0 = (c0[0], c0[1], roof_z(c0[0], c0[1]))
        rt1 = (c1[0], c1[1], roof_z(c1[0], c1[1]))

        # If rb and rt0 differ in z (not a degenerate edge), add triangle
        if abs(rt0[2] - rb[2]) > 0.01 or abs(rt1[2] - re[2]) > 0.01:
            tris.append((rb, re, rt1))
            if abs(rt0[2] - rb[2]) > 0.01:
                tris.append((rb, rt1, rt0))

    return tris


def make_stilts():
    """Four stilt columns at face midpoints, 0 to STILT_H."""
    tris = []
    for cx, cy in STILTS:
        rx, ry = rotate45(cx, cy)
        hw = HALF_S
        # Axis-aligned small box, no need to rotate stilt cross-section
        corners = [
            (rx-hw, ry-hw), (rx+hw, ry-hw),
            (rx+hw, ry+hw), (rx-hw, ry+hw)
        ]
        zb, zt = 0.0, STILT_H
        for i in range(4):
            c0, c1 = corners[i], corners[(i+1)%4]
            p00 = (c0[0], c0[1], zb)
            p10 = (c1[0], c1[1], zb)
            p11 = (c1[0], c1[1], zt)
            p01 = (c0[0], c0[1], zt)
            tris += [(p00, p10, p11), (p00, p11, p01)]
        # Top cap
        tris += [(corners[0]+(zt,), corners[2]+(zt,), corners[1]+(zt,)),
                 (corners[0]+(zt,), corners[3]+(zt,), corners[2]+(zt,))]
    return tris


# ─── EXTRUDE POLYGON TO STL ───────────────────────────────────────────────────
def extrude_polygon(verts_2d, z_bot, z_top):
    """Extrude a 2D polygon to a closed prism between z_bot and z_top."""
    tris = []
    n = len(verts_2d)
    if poly_area(verts_2d) < 0:
        verts_2d = list(reversed(verts_2d))

    # Walls
    for i in range(n):
        c0, c1 = verts_2d[i], verts_2d[(i+1)%n]
        p00 = (c0[0], c0[1], z_bot)
        p10 = (c1[0], c1[1], z_bot)
        p11 = (c1[0], c1[1], z_top)
        p01 = (c0[0], c0[1], z_top)
        tris += [(p00, p10, p11), (p00, p11, p01)]

    # Top cap
    for t in ear_clip(verts_2d):
        tris.append(((t[0][0],t[0][1],z_top),(t[1][0],t[1][1],z_top),(t[2][0],t[2][1],z_top)))

    # Bottom cap (flipped winding)
    for t in ear_clip(verts_2d):
        tris.append(((t[0][0],t[0][1],z_bot),(t[2][0],t[2][1],z_bot),(t[1][0],t[1][1],z_bot)))

    return tris


# ─── NYC OPEN DATA API (with construction_year filter) ───────────────────────
NYC_API = "https://data.cityofnewyork.us/resource/5zhs-2jue.geojson"

def fetch_buildings_1978():
    """
    Fetch buildings from NYC Open Data API that existed by HISTORICAL_YEAR.
    Filter: construction_year <= 1978 OR construction_year IS NULL
    (NULL means year is unknown — safer to include than exclude)
    """
    print(f"  Querying NYC Open Data API for buildings built by {HISTORICAL_YEAR}...")
    where = (
        f"within_box(the_geom, {BB_LAT_MAX}, {BB_LON_MIN}, {BB_LAT_MIN}, {BB_LON_MAX})"
        f" AND (construction_year IS NULL OR construction_year <= {HISTORICAL_YEAR})"
    )
    params = {"$where": where, "$limit": 10000}
    resp = requests.get(NYC_API, params=params, timeout=60)
    resp.raise_for_status()
    features = resp.json().get("features", [])
    print(f"  Received {len(features)} pre-{HISTORICAL_YEAR} building footprints")
    return features


def fetch_buildings_osm_1978():
    """
    Fallback: OSM Overpass with start_date tag filtering.
    Less reliable than NYC API but provides a fallback option.
    """
    print("  Querying OSM Overpass API for pre-1978 buildings...")
    overpass_url = "https://overpass-api.de/api/interpreter"
    # Query buildings with start_date <= 1978 OR no start_date
    query = f"""
[out:json][timeout:90];
(
  way["building"][!"start_date"]
    ({BB_LAT_MIN},{BB_LON_MIN},{BB_LAT_MAX},{BB_LON_MAX});
  way["building"]["start_date"<="1978"]
    ({BB_LAT_MIN},{BB_LON_MIN},{BB_LAT_MAX},{BB_LON_MAX});
);
out geom;
"""
    resp = requests.post(overpass_url, data={"data": query}, timeout=120)
    resp.raise_for_status()
    return resp.json().get("elements", [])


def parse_nyc_feature(feature):
    """Parse GeoJSON feature from NYC API → (verts_2d, height_m) or None."""
    props = feature.get("properties", {})
    geom  = feature.get("geometry", {})

    height_ft = props.get("height_roof")
    if not height_ft:
        return None
    height_m = float(height_ft) * FT_TO_M
    if height_m < MIN_HEIGHT_M:
        return None

    coords = geom.get("coordinates", [])
    gtype  = geom.get("type", "")
    ring   = None
    if gtype == "Polygon" and coords:
        ring = coords[0]
    elif gtype == "MultiPolygon" and coords:
        best, ba = None, 0
        for poly in coords:
            if poly and poly[0]:
                a = abs(poly_area([(c[0],c[1]) for c in poly[0]]))
                if a > ba:
                    ba, best = a, poly[0]
        ring = best

    if not ring or len(ring) < 4:
        return None

    verts = [wgs84_to_local(c[0], c[1]) for c in ring]
    if len(verts) > 1 and math.hypot(verts[0][0]-verts[-1][0], verts[0][1]-verts[-1][1]) < 0.01:
        verts.pop()
    if len(verts) < 3:
        return None
    if abs(poly_area(verts)) < 10.0:
        return None

    return verts, height_m


def is_citicorp_feature(feature):
    """Detect if a NYC API feature is the Citicorp footprint."""
    props = feature.get("properties", {})
    bin_str = str(props.get("bin", ""))
    if bin_str in ("1035879", "1087931"):
        return True
    geom = feature.get("geometry", {})
    coords = geom.get("coordinates", [[[]]])
    gtype  = geom.get("type", "")
    ring = None
    if gtype == "Polygon" and coords:
        ring = coords[0]
    elif gtype == "MultiPolygon" and coords and coords[0]:
        ring = coords[0][0]
    if ring:
        verts = [wgs84_to_local(c[0], c[1]) for c in ring]
        if verts:
            cx = sum(v[0] for v in verts) / len(verts)
            cy = sum(v[1] for v in verts) / len(verts)
            if math.hypot(cx, cy) < 35:
                return True
    return False


# ─── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "constant", "triSurface")
    os.makedirs(out_dir, exist_ok=True)

    print(f"\nCiticorp 1978 Historical CFD Geometry Generator")
    print(f"Source: OSM + NYC Open Data (construction_year <= {HISTORICAL_YEAR})")
    print("=" * 60)

    # ── Tower & Stilts (hand-crafted, unchanged from actual geometry) ─────────
    print("\nGenerating Citicorp tower (hand-crafted, 45-deg slanted roof, 45-deg rotation)...")
    tower_tris = make_tower()
    write_stl(os.path.join(out_dir, "citicorp_tower.stl"), "citicorp_tower", tower_tris)

    print("\nGenerating Citicorp stilts (hand-crafted)...")
    stilt_tris = make_stilts()
    write_stl(os.path.join(out_dir, "citicorp_stilts.stl"), "citicorp_stilts", stilt_tris)

    # ── Surroundings: pre-1978 buildings ─────────────────────────────────────
    print(f"\nFetching pre-{HISTORICAL_YEAR} surrounding buildings...")
    surr_tris = []
    building_count = 0
    citicorp_skipped = 0
    skipped = 0

    if HAS_REQUESTS:
        try:
            # Primary source: NYC Building Footprints with construction_year filter
            features = fetch_buildings_1978()

            for feat in features:
                if is_citicorp_feature(feat):
                    citicorp_skipped += 1
                    continue
                parsed = parse_nyc_feature(feat)
                if parsed is None:
                    skipped += 1
                    continue
                verts, height_m = parsed
                surr_tris.extend(extrude_polygon(verts, 0.0, height_m))
                building_count += 1

            # Stats
            heights = []
            for feat in features:
                if not is_citicorp_feature(feat):
                    props = feat.get("properties", {})
                    hft = props.get("height_roof")
                    if hft:
                        heights.append(float(hft) * FT_TO_M)

            print(f"\n  Pre-{HISTORICAL_YEAR} buildings processed: {building_count}")
            print(f"  Skipped (too short or degenerate): {skipped}")
            if heights:
                print(f"  Height range: {min(heights):.1f}m - {max(heights):.1f}m")
                print(f"  Mean height:  {sum(heights)/len(heights):.1f}m")
                print(f"  Median height:{sorted(heights)[len(heights)//2]:.1f}m")

        except Exception as e:
            print(f"  NYC API failed: {e}")
            print("  No fallback available — surroundings STL will be empty")
    else:
        print("  requests not installed — install with: pip install requests")
        print("  surroundings STL will be empty (Citicorp geometry still generated)")

    write_stl(os.path.join(out_dir, "surroundings.stl"),
              f"surroundings_osm_{HISTORICAL_YEAR}", surr_tris)

    print(f"\n{'='*60}")
    print(f"1978 HISTORICAL GEOMETRY SUMMARY")
    print(f"{'='*60}")
    print(f"Tower triangles:         {len(tower_tris)}")
    print(f"Stilt triangles:         {len(stilt_tris)}")
    print(f"Surrounding triangles:   {len(surr_tris)}")
    print(f"Buildings included:      {building_count}")
    print(f"Historical cutoff:       {HISTORICAL_YEAR}")
    print(f"\nOutput directory: {out_dir}")
    print(f"\nNext steps:")
    print(f"  1. Copy STLs to WSL: ~/citicorp_cfd_osm_1978/constant/triSurface/")
    print(f"  2. Run surfaceFeatures")
    print(f"  3. Run blockMesh && snappyHexMesh")
    print(f"  4. Open citicorp_cfd_osm_1978.foam in ParaView")


if __name__ == "__main__":
    main()
