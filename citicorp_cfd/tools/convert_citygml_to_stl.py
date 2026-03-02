#!/usr/bin/env python3
"""
Convert Citicorp CityGML LOD2 building to STL for OpenFOAM.

Source: NYC 3D Building Model (georocket enhanced), DA12
        BIN 1036474 — Citigroup Center (1978, 60 floors, Individual Landmark)

Coordinate system:
  Input:  EPSG:2263 (NY State Plane Long Island, US Survey Feet)
  Output: Local meters, centered at tower centroid, Z = 0 at ground

The CityGML covers the entire Citicorp complex (tower + church + atrium).
We separate it into:
  - citicorp_tower.stl   = main tower body (walls/roof above stilt level)
  - citicorp_complex.stl  = full complex including lower structures
  - citicorp_stilts.stl   = kept from existing geometry (not in aerial CityGML)
"""

import re
import struct
import math
import os

# --- Constants ---
US_SURVEY_FT_TO_M = 1200.0 / 3937.0  # = 0.3048006096... m per ft

# Ground elevation from CityGML GroundSurface
Z_GROUND = 25.78  # feet above sea level

# Tower geometry thresholds
STILT_TOP_FT = Z_GROUND + 114.0  # 114 ft stilt height = ~34.75 m
TOWER_MIN_Z_FT = 800.0  # Walls reaching above this are "tower"

# --- Parse CityGML ---

def parse_citygml(filepath):
    """Parse CityGML file and extract all polygons with surface type."""
    with open(filepath, 'r') as f:
        content = f.read()

    surfaces = []

    # Find all boundedBy blocks with surface type
    # Pattern: <bldg:XXXSurface ...> ... <gml:posList>coords</gml:posList> ...
    pattern = re.compile(
        r'<bldg:(GroundSurface|RoofSurface|WallSurface)\s[^>]*>.*?'
        r'<gml:posList>(.*?)</gml:posList>.*?'
        r'</bldg:\1>',
        re.DOTALL
    )

    for match in pattern.finditer(content):
        stype = match.group(1)
        coords_str = match.group(2).strip()
        coords = [float(x) for x in coords_str.split()]

        # Coordinates are X Y Z triplets
        n = len(coords) // 3
        if n < 3:
            continue

        vertices = []
        for i in range(n):
            x_ft = coords[i * 3]
            y_ft = coords[i * 3 + 1]
            z_ft = coords[i * 3 + 2]
            vertices.append((x_ft, y_ft, z_ft))

        # Remove closing vertex if it duplicates the first
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


def convert_coords(vertices_ft, x0, y0, z_ground):
    """Convert EPSG:2263 feet to local meters centered at (x0, y0), z=0 at ground."""
    result = []
    for x_ft, y_ft, z_ft in vertices_ft:
        x_m = (x_ft - x0) * US_SURVEY_FT_TO_M
        y_m = (y_ft - y0) * US_SURVEY_FT_TO_M
        z_m = (z_ft - z_ground) * US_SURVEY_FT_TO_M
        result.append((x_m, y_m, z_m))
    return result


def triangulate_polygon(vertices):
    """Fan triangulation from first vertex. Works for convex and mildly concave polygons."""
    triangles = []
    if len(vertices) < 3:
        return triangles
    v0 = vertices[0]
    for i in range(1, len(vertices) - 1):
        triangles.append((v0, vertices[i], vertices[i + 1]))
    return triangles


def compute_normal(v0, v1, v2):
    """Compute outward-facing normal for a triangle."""
    # Edge vectors
    e1 = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
    e2 = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])
    # Cross product
    nx = e1[1] * e2[2] - e1[2] * e2[1]
    ny = e1[2] * e2[0] - e1[0] * e2[2]
    nz = e1[0] * e2[1] - e1[1] * e2[0]
    # Normalize
    length = math.sqrt(nx * nx + ny * ny + nz * nz)
    if length < 1e-12:
        return (0.0, 0.0, 1.0)
    return (nx / length, ny / length, nz / length)


def write_binary_stl(filepath, solid_name, triangles_with_normals):
    """Write binary STL file."""
    n_tri = len(triangles_with_normals)
    with open(filepath, 'wb') as f:
        # 80-byte header
        header = solid_name.encode('ascii')[:80].ljust(80, b'\0')
        f.write(header)
        # Number of triangles
        f.write(struct.pack('<I', n_tri))
        # Each triangle: normal (3 floats) + 3 vertices (9 floats) + attribute (2 bytes)
        for normal, v0, v1, v2 in triangles_with_normals:
            f.write(struct.pack('<fff', *normal))
            f.write(struct.pack('<fff', *v0))
            f.write(struct.pack('<fff', *v1))
            f.write(struct.pack('<fff', *v2))
            f.write(struct.pack('<H', 0))  # attribute byte count

    print(f"  {filepath}: {n_tri} triangles")


def write_ascii_stl(filepath, solid_name, triangles_with_normals):
    """Write ASCII STL file (for debugging / readability)."""
    n_tri = len(triangles_with_normals)
    with open(filepath, 'w') as f:
        f.write(f"solid {solid_name}\n")
        for normal, v0, v1, v2 in triangles_with_normals:
            f.write(f"  facet normal {normal[0]:.6f} {normal[1]:.6f} {normal[2]:.6f}\n")
            f.write("    outer loop\n")
            for v in [v0, v1, v2]:
                f.write(f"      vertex {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {solid_name}\n")

    print(f"  {filepath}: {n_tri} triangles (ASCII)")


def surfaces_to_triangles(surfaces, x0, y0, z_ground):
    """Convert parsed surfaces to triangles in local meters."""
    all_triangles = []
    for surf in surfaces:
        verts_m = convert_coords(surf['vertices_ft'], x0, y0, z_ground)
        tris = triangulate_polygon(verts_m)
        for v0, v1, v2 in tris:
            normal = compute_normal(v0, v1, v2)
            all_triangles.append((normal, v0, v1, v2))
    return all_triangles


def classify_surfaces(surfaces):
    """Classify surfaces into tower vs lower structures."""
    tower_surfaces = []
    lower_surfaces = []

    for surf in surfaces:
        # Tower = walls or roof surfaces where the max Z > 800 ft
        # This captures the main tower shaft and crown
        if surf['max_z_ft'] > TOWER_MIN_Z_FT:
            tower_surfaces.append(surf)
        else:
            lower_surfaces.append(surf)

    return tower_surfaces, lower_surfaces


def find_tower_centroid(tower_surfaces):
    """Find the tower center as the bounding-box midpoint of tower walls."""
    # Use wall surfaces with high max_z (the tall tower walls)
    tall_walls = [s for s in tower_surfaces
                  if s['type'] == 'WallSurface' and s['max_z_ft'] > TOWER_MIN_Z_FT]

    if not tall_walls:
        tall_walls = tower_surfaces

    # Collect all unique XY positions from tall walls
    all_x = []
    all_y = []
    for s in tall_walls:
        for v in s['vertices_ft']:
            all_x.append(v[0])
            all_y.append(v[1])

    # Bounding box midpoint (more robust than centroid for irregular shapes)
    cx = (min(all_x) + max(all_x)) / 2.0
    cy = (min(all_y) + max(all_y)) / 2.0
    return cx, cy


def print_surface_stats(surfaces, label):
    """Print statistics about a set of surfaces."""
    if not surfaces:
        print(f"  {label}: 0 surfaces")
        return

    types = {}
    z_min = float('inf')
    z_max = float('-inf')

    for s in surfaces:
        types[s['type']] = types.get(s['type'], 0) + 1
        z_min = min(z_min, s['min_z_ft'])
        z_max = max(z_max, s['max_z_ft'])

    z_min_m = (z_min - Z_GROUND) * US_SURVEY_FT_TO_M
    z_max_m = (z_max - Z_GROUND) * US_SURVEY_FT_TO_M

    print(f"  {label}: {len(surfaces)} surfaces")
    for t, c in sorted(types.items()):
        print(f"    {t}: {c}")
    print(f"    Z range: {z_min:.1f} - {z_max:.1f} ft ({z_min_m:.1f} - {z_max_m:.1f} m above ground)")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    gml_path = os.path.join(script_dir, '_citygml', 'citicorp_building.gml')
    stl_dir = os.path.join(script_dir, 'constant', 'triSurface')

    print("=" * 60)
    print("CityGML LOD2 -> STL Converter for Citicorp Center")
    print("=" * 60)

    # Parse
    print(f"\nParsing: {gml_path}")
    surfaces = parse_citygml(gml_path)
    print(f"  Found {len(surfaces)} polygons total")

    # Classify
    print("\nClassifying surfaces:")
    tower_surfs, lower_surfs = classify_surfaces(surfaces)
    print_surface_stats(tower_surfs, "Tower (max Z > 800 ft)")
    print_surface_stats(lower_surfs, "Lower structures (max Z <= 800 ft)")

    # Find tower centroid for centering
    x0, y0 = find_tower_centroid(tower_surfs)
    print(f"\n  Tower centroid (EPSG:2263): X={x0:.1f}, Y={y0:.1f} ft")
    print(f"  Ground elevation: {Z_GROUND:.2f} ft")

    # Convert to triangles
    print("\nTriangulating...")
    tower_tris = surfaces_to_triangles(tower_surfs, x0, y0, Z_GROUND)
    lower_tris = surfaces_to_triangles(lower_surfs, x0, y0, Z_GROUND)
    all_tris = surfaces_to_triangles(surfaces, x0, y0, Z_GROUND)

    print(f"  Tower: {len(tower_tris)} triangles")
    print(f"  Lower: {len(lower_tris)} triangles")
    print(f"  Total: {len(all_tris)} triangles")

    # Compute bounding boxes
    print("\nBounding boxes (local meters, centered at tower):")
    for label, tris in [("Tower", tower_tris), ("Full complex", all_tris)]:
        if not tris:
            continue
        all_verts = []
        for _, v0, v1, v2 in tris:
            all_verts.extend([v0, v1, v2])
        xs = [v[0] for v in all_verts]
        ys = [v[1] for v in all_verts]
        zs = [v[2] for v in all_verts]
        print(f"  {label}:")
        print(f"    X: {min(xs):.1f} to {max(xs):.1f} m (span {max(xs)-min(xs):.1f} m)")
        print(f"    Y: {min(ys):.1f} to {max(ys):.1f} m (span {max(ys)-min(ys):.1f} m)")
        print(f"    Z: {min(zs):.1f} to {max(zs):.1f} m (height {max(zs):.1f} m)")

    # Write STL files
    print(f"\nWriting STL files to: {stl_dir}")
    os.makedirs(stl_dir, exist_ok=True)

    # Binary STLs (for OpenFOAM)
    write_binary_stl(
        os.path.join(stl_dir, 'citicorp_tower_lod2_bin.stl'),
        'citicorp_tower',
        tower_tris
    )
    write_binary_stl(
        os.path.join(stl_dir, 'citicorp_complex_lod2_bin.stl'),
        'citicorp_complex',
        all_tris
    )

    # ASCII STLs (for inspection)
    write_ascii_stl(
        os.path.join(stl_dir, 'citicorp_tower_lod2.stl'),
        'citicorp_tower',
        tower_tris
    )
    write_ascii_stl(
        os.path.join(stl_dir, 'citicorp_complex_lod2.stl'),
        'citicorp_complex',
        all_tris
    )

    # Summary comparison with current simplified geometry
    print("\n" + "=" * 60)
    print("Comparison with current simplified STL:")
    print("=" * 60)
    print(f"  Current citicorp_tower.stl:  19 triangles (4 walls + diamond roof + crown)")
    print(f"  Current citicorp_stilts.stl: 40 triangles (4 axis-aligned boxes)")
    print(f"  LOD2 tower:                  {len(tower_tris)} triangles (actual building outline)")
    print(f"  LOD2 full complex:           {len(all_tris)} triangles (tower + church + atrium)")
    print()
    print("  NOTE: citicorp_stilts.stl is NOT in CityGML (hidden from aerial survey).")
    print("        Keep existing stilts or improve them separately.")

    print("\nDone!")


if __name__ == '__main__':
    main()
