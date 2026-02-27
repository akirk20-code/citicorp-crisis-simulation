#!/usr/bin/env python3
"""
Extract LOD2 buildings from NYC 3D Building Model for any location.

Usage:
    python extract_buildings.py --lat 40.7580 --lon -73.9690 --radius 300
    python extract_buildings.py --lat 40.7580 --lon -73.9690 --radius 300 --format stl
    python extract_buildings.py --lat 40.7484 --lon -73.9857 --radius 200 --name "empire_state"

Data source: NYC 3D Building Model (georocket enhanced, v20v5)
    CityGML LOD2 with RoofSurface, WallSurface, GroundSurface
    Coordinate system: EPSG:2263 (NY State Plane Long Island, US Survey Feet)
    20 DA districts covering all of NYC

How it works:
    1. Convert lat/lon to EPSG:2263 (NY State Plane feet)
    2. Auto-detect which DA district file(s) to search
    3. Stream through GML, extract buildings within radius
    4. Output: CityGML fragments + optional STL conversion

Requirements: Python 3.8+ (stdlib only — no pip installs needed)
"""

import argparse
import math
import os
import re
import struct
import sys
import zipfile
import time

# If requests is available, use it for auto-download; otherwise manual instructions
try:
    import urllib.request
    HAS_URLLIB = True
except ImportError:
    HAS_URLLIB = False


# ============================================================
# Constants
# ============================================================

US_SURVEY_FT_TO_M = 1200.0 / 3937.0   # 0.30480061 m/ft
M_TO_US_SURVEY_FT = 3937.0 / 1200.0   # 3.28083333 ft/m

# EPSG:2263 projection parameters (NY State Plane Long Island, NAD83)
# Using simplified equirectangular approximation from WGS84 lat/lon
# This is accurate to ~0.1% within NYC (good enough for bounding box queries)
# For exact conversion, use pyproj — but we avoid external dependencies.
#
# Reference point for the approximation (center of Manhattan):
REF_LAT = 40.7580   # degrees
REF_LON = -73.9700  # degrees
REF_X = 992506.0    # EPSG:2263 X (feet) — Citicorp Center
REF_Y = 215620.0    # EPSG:2263 Y (feet) — Citicorp Center

# Approximate scale factors at NYC latitude
# 1 degree latitude  ~ 111,000 m ~ 364,173 ft
# 1 degree longitude ~ 111,000 * cos(40.758) ~ 84,290 m ~ 276,540 ft
FT_PER_DEG_LAT = 364173.0
FT_PER_DEG_LON = 276540.0

# DA district approximate bounding boxes (EPSG:2263 X range, Y range)
# These are rough envelopes to quickly identify which DA file to search.
# Built from the envelope elements in each DA GML file.
DA_ENVELOPES = {
    # DA: (x_min, x_max, y_min, y_max)  -- all in EPSG:2263 feet
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


# ============================================================
# Coordinate conversion (WGS84 lat/lon <-> EPSG:2263)
# ============================================================

def latlon_to_stateplane(lat, lon):
    """Convert WGS84 lat/lon to approximate EPSG:2263 (NY State Plane feet).

    Uses linear approximation around the reference point. Accurate to ~1m
    within 10km of the reference point, which is sufficient for building queries.
    """
    x = REF_X + (lon - REF_LON) * FT_PER_DEG_LON
    y = REF_Y + (lat - REF_LAT) * FT_PER_DEG_LAT
    return x, y


def stateplane_to_latlon(x, y):
    """Convert EPSG:2263 feet to approximate WGS84 lat/lon."""
    lon = REF_LON + (x - REF_X) / FT_PER_DEG_LON
    lat = REF_LAT + (y - REF_Y) / FT_PER_DEG_LAT
    return lat, lon


# ============================================================
# DA district detection
# ============================================================

def find_da_districts(cx, cy, radius_ft):
    """Find which DA district(s) contain the query bounding box."""
    x_min = cx - radius_ft
    x_max = cx + radius_ft
    y_min = cy - radius_ft
    y_max = cy + radius_ft

    matches = []
    for da, (dx_min, dx_max, dy_min, dy_max) in DA_ENVELOPES.items():
        # Check overlap
        if x_max >= dx_min and x_min <= dx_max and y_max >= dy_min and y_min <= dy_max:
            matches.append(da)

    return matches


def get_da_gml_path(da, citygml_dir):
    """Get path to DA GML file, downloading if needed."""
    da_dir = os.path.join(citygml_dir, da)
    gml_path = os.path.join(da_dir, f"{da}_3D_Buildings_Merged.gml")
    zip_path = os.path.join(citygml_dir, f"{da}.gml.zip")

    if os.path.exists(gml_path):
        return gml_path

    # Check if zip exists but not extracted
    if os.path.exists(zip_path):
        print(f"  Extracting {zip_path}...")
        os.makedirs(da_dir, exist_ok=True)
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(da_dir)
        if os.path.exists(gml_path):
            return gml_path

    # Need to download
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
# Building extraction (streaming parser)
# ============================================================

def extract_buildings_in_box(gml_path, x_min, x_max, y_min, y_max):
    """Stream through GML file and extract buildings within bounding box.

    Uses a line-by-line state machine to avoid loading the entire file into memory.
    Each DA GML file can be 200-800MB uncompressed.
    """
    buildings = []
    current_building = []
    in_building = False
    building_in_box = False
    building_count = 0
    match_count = 0

    # Pre-compile patterns for speed
    building_start = re.compile(r'<bldg:Building\s')
    building_end = '</bldg:Building>'
    poslist_pat = re.compile(r'<gml:posList>(.*?)</gml:posList>')

    file_size = os.path.getsize(gml_path)
    print(f"  Scanning {gml_path}")
    print(f"  File size: {file_size / 1024 / 1024:.0f} MB")
    print(f"  Query box: X=[{x_min:.0f}, {x_max:.0f}], Y=[{y_min:.0f}, {y_max:.0f}]")

    t0 = time.time()
    bytes_read = 0
    last_report = 0

    with open(gml_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            bytes_read += len(line.encode('utf-8', errors='replace'))

            # Progress reporting every ~50MB
            if bytes_read - last_report > 50_000_000:
                elapsed = time.time() - t0
                pct = 100 * bytes_read / file_size if file_size > 0 else 0
                rate = bytes_read / 1024 / 1024 / elapsed if elapsed > 0 else 0
                print(f"    {pct:.0f}% ({bytes_read/1024/1024:.0f} MB, "
                      f"{rate:.0f} MB/s, {building_count} buildings, "
                      f"{match_count} matches)")
                last_report = bytes_read

            if not in_building:
                if building_start.search(line):
                    in_building = True
                    building_in_box = False
                    current_building = [line]
                    building_count += 1
            else:
                current_building.append(line)

                # Check coordinates in posList elements
                if not building_in_box:
                    m = poslist_pat.search(line)
                    if m:
                        coords = m.group(1).strip().split()
                        # Coordinates are X Y Z triplets
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
    print(f"  Done: {building_count} buildings scanned, "
          f"{match_count} in box ({elapsed:.1f}s)")

    return buildings


# ============================================================
# CityGML to STL conversion (reused from _convert_citygml_to_stl.py)
# ============================================================

def parse_building_surfaces(gml_text):
    """Parse surfaces from a single building's GML text."""
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


def triangulate_polygon(vertices):
    """Fan triangulation from first vertex."""
    triangles = []
    if len(vertices) < 3:
        return triangles
    v0 = vertices[0]
    for i in range(1, len(vertices) - 1):
        triangles.append((v0, vertices[i], vertices[i + 1]))
    return triangles


def compute_normal(v0, v1, v2):
    """Compute triangle normal."""
    e1 = (v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2])
    e2 = (v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2])
    nx = e1[1]*e2[2] - e1[2]*e2[1]
    ny = e1[2]*e2[0] - e1[0]*e2[2]
    nz = e1[0]*e2[1] - e1[1]*e2[0]
    length = math.sqrt(nx*nx + ny*ny + nz*nz)
    if length < 1e-12:
        return (0.0, 0.0, 1.0)
    return (nx/length, ny/length, nz/length)


def convert_coords(vertices_ft, x0, y0, z_ground):
    """Convert EPSG:2263 feet to local meters."""
    return [((x-x0)*US_SURVEY_FT_TO_M, (y-y0)*US_SURVEY_FT_TO_M,
             (z-z_ground)*US_SURVEY_FT_TO_M)
            for x, y, z in vertices_ft]


def buildings_to_stl(buildings_gml, center_x, center_y, output_path, solid_name):
    """Convert extracted buildings to binary STL centered at given point."""
    # Find ground elevation from GroundSurface elements
    z_grounds = []
    for bldg_text in buildings_gml:
        surfs = parse_building_surfaces(bldg_text)
        for s in surfs:
            if s['type'] == 'GroundSurface':
                z_grounds.extend(v[2] for v in s['vertices_ft'])
    z_ground = min(z_grounds) if z_grounds else 0.0

    # Convert all buildings to triangles
    all_triangles = []
    for bldg_text in buildings_gml:
        surfs = parse_building_surfaces(bldg_text)
        for surf in surfs:
            verts_m = convert_coords(surf['vertices_ft'], center_x, center_y, z_ground)
            tris = triangulate_polygon(verts_m)
            for v0, v1, v2 in tris:
                normal = compute_normal(v0, v1, v2)
                all_triangles.append((normal, v0, v1, v2))

    if not all_triangles:
        print(f"  No triangles to write!")
        return 0

    # Write binary STL
    n_tri = len(all_triangles)
    with open(output_path, 'wb') as f:
        header = solid_name.encode('ascii')[:80].ljust(80, b'\0')
        f.write(header)
        f.write(struct.pack('<I', n_tri))
        for normal, v0, v1, v2 in all_triangles:
            f.write(struct.pack('<fff', *normal))
            f.write(struct.pack('<fff', *v0))
            f.write(struct.pack('<fff', *v1))
            f.write(struct.pack('<fff', *v2))
            f.write(struct.pack('<H', 0))

    print(f"  Wrote {output_path}: {n_tri} triangles")
    return n_tri


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Extract LOD2 buildings from NYC 3D Building Model',
        epilog='Example: python extract_buildings.py --lat 40.7580 --lon -73.9690 --radius 300'
    )
    parser.add_argument('--lat', type=float, required=True,
                        help='Center latitude (WGS84)')
    parser.add_argument('--lon', type=float, required=True,
                        help='Center longitude (WGS84)')
    parser.add_argument('--radius', type=float, default=300,
                        help='Extraction radius in meters (default: 300)')
    parser.add_argument('--name', type=str, default='buildings',
                        help='Output name prefix (default: buildings)')
    parser.add_argument('--format', choices=['gml', 'stl', 'both'], default='both',
                        help='Output format (default: both)')
    parser.add_argument('--data-dir', type=str, default=None,
                        help='CityGML data directory (default: _citygml/ next to script)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory (default: current directory)')

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    citygml_dir = args.data_dir or os.path.join(script_dir, '_citygml')
    output_dir = args.output_dir or os.getcwd()

    print("=" * 65)
    print("NYC 3D Building Model -- LOD2 Building Extractor")
    print("=" * 65)
    print(f"  Center:  {args.lat:.6f}, {args.lon:.6f}")
    print(f"  Radius:  {args.radius:.0f} m")
    print(f"  Name:    {args.name}")
    print(f"  Format:  {args.format}")

    # Convert to State Plane
    cx, cy = latlon_to_stateplane(args.lat, args.lon)
    radius_ft = args.radius * M_TO_US_SURVEY_FT
    print(f"\n  EPSG:2263: X={cx:.0f}, Y={cy:.0f} ft")
    print(f"  Radius:    {radius_ft:.0f} ft")

    # Find DA districts
    x_min = cx - radius_ft
    x_max = cx + radius_ft
    y_min = cy - radius_ft
    y_max = cy + radius_ft

    districts = find_da_districts(cx, cy, radius_ft)
    if not districts:
        print("\n  ERROR: Location is outside NYC 3D Building Model coverage.")
        print("  The model covers NYC's 5 boroughs only.")
        sys.exit(1)

    print(f"\n  DA districts to search: {', '.join(districts)}")

    # Extract buildings from each district
    all_buildings = []
    for da in districts:
        print(f"\n--- {da} ---")
        gml_path = get_da_gml_path(da, citygml_dir)
        if gml_path is None:
            print(f"  Skipping {da} (file not available)")
            continue
        buildings = extract_buildings_in_box(gml_path, x_min, x_max, y_min, y_max)
        all_buildings.extend(buildings)

    if not all_buildings:
        print("\n  No buildings found in the search area.")
        sys.exit(0)

    # Print building summary
    print(f"\n{'=' * 65}")
    print(f"Found {len(all_buildings)} buildings")
    print(f"{'=' * 65}")

    heights = []
    for bldg_text in all_buildings:
        bin_val = get_building_attr(bldg_text, 'BIN')
        name_val = get_building_attr(bldg_text, 'name') or ''
        year_val = get_building_attr(bldg_text, 'yearbuilt') or '?'
        floors_val = get_building_attr(bldg_text, 'numfloors') or '?'
        landmark = get_building_attr(bldg_text, 'landmark') or ''

        surfs = parse_building_surfaces(bldg_text)
        if surfs:
            z_vals = [v[2] for s in surfs for v in s['vertices_ft']]
            h_ft = max(z_vals) - min(z_vals)
            h_m = h_ft * US_SURVEY_FT_TO_M
            heights.append(h_m)
        else:
            h_m = 0
            heights.append(0)

        # Print notable buildings (tall or named)
        if h_m > 50 or landmark or (name_val and name_val != '0'):
            tag = f" ** {landmark}" if landmark else ""
            print(f"  BIN {bin_val}: {h_m:.0f}m, {floors_val} floors, "
                  f"built {year_val}{tag}")

    if heights:
        print(f"\n  Height stats: min={min(heights):.0f}m, "
              f"max={max(heights):.0f}m, "
              f"mean={sum(heights)/len(heights):.0f}m, "
              f"median={sorted(heights)[len(heights)//2]:.0f}m")

    # Write GML output
    if args.format in ('gml', 'both'):
        gml_out = os.path.join(output_dir, f"{args.name}.gml")
        with open(gml_out, 'w', encoding='utf-8') as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write(f'<!-- {len(all_buildings)} buildings extracted from '
                    f'NYC 3D Building Model -->\n')
            f.write(f'<!-- Center: {args.lat:.6f}, {args.lon:.6f}, '
                    f'radius: {args.radius:.0f}m -->\n')
            f.write('<buildings>\n')
            for bldg_text in all_buildings:
                f.write(bldg_text)
                f.write('\n')
            f.write('</buildings>\n')
        print(f"\n  GML: {gml_out}")

    # Write STL output
    if args.format in ('stl', 'both'):
        stl_out = os.path.join(output_dir, f"{args.name}.stl")
        n_tri = buildings_to_stl(all_buildings, cx, cy, stl_out, args.name)

        # Also compute bounding box
        if n_tri > 0:
            # Read back for stats
            with open(stl_out, 'rb') as f:
                f.read(80)  # header
                nt = struct.unpack('<I', f.read(4))[0]
                xs, ys, zs = [], [], []
                for _ in range(nt):
                    f.read(12)  # normal
                    for _ in range(3):
                        vx, vy, vz = struct.unpack('<fff', f.read(12))
                        xs.append(vx)
                        ys.append(vy)
                        zs.append(vz)
                    f.read(2)  # attribute
            print(f"  STL bounding box (local meters from center):")
            print(f"    X: {min(xs):.1f} to {max(xs):.1f} m")
            print(f"    Y: {min(ys):.1f} to {max(ys):.1f} m")
            print(f"    Z: {min(zs):.1f} to {max(zs):.1f} m")

    print(f"\nDone! {len(all_buildings)} buildings extracted.")


if __name__ == '__main__':
    main()
