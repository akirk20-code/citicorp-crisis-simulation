#!/usr/bin/env python3
"""
Compare current simplified STL with NYC LOD2 CityGML-derived STL.
Produces a text report showing differences in geometry.
"""

import os
import re
import math

def parse_ascii_stl(filepath):
    """Parse ASCII STL and return list of triangles with normals."""
    triangles = []
    with open(filepath, 'r') as f:
        content = f.read()

    # Find all facet blocks
    pattern = re.compile(
        r'facet normal\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+'
        r'outer loop\s+'
        r'vertex\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+'
        r'vertex\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+'
        r'vertex\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+'
        r'endloop',
        re.IGNORECASE
    )

    for m in pattern.finditer(content):
        vals = [float(m.group(i)) for i in range(1, 13)]
        normal = (vals[0], vals[1], vals[2])
        v0 = (vals[3], vals[4], vals[5])
        v1 = (vals[6], vals[7], vals[8])
        v2 = (vals[9], vals[10], vals[11])
        triangles.append((normal, v0, v1, v2))

    return triangles


def bbox(triangles):
    """Compute bounding box from triangles."""
    all_v = []
    for _, v0, v1, v2 in triangles:
        all_v.extend([v0, v1, v2])
    xs = [v[0] for v in all_v]
    ys = [v[1] for v in all_v]
    zs = [v[2] for v in all_v]
    return {
        'xmin': min(xs), 'xmax': max(xs),
        'ymin': min(ys), 'ymax': max(ys),
        'zmin': min(zs), 'zmax': max(zs),
    }


def unique_z_levels(triangles, tol=0.5):
    """Find unique Z levels in the triangles."""
    z_vals = set()
    for _, v0, v1, v2 in triangles:
        for v in [v0, v1, v2]:
            z_vals.add(round(v[2], 1))
    # Cluster nearby values
    z_sorted = sorted(z_vals)
    clusters = []
    for z in z_sorted:
        if not clusters or z - clusters[-1] > tol:
            clusters.append(z)
    return clusters


def footprint_at_z(triangles, z_target, tol=1.0):
    """Extract XY vertices of triangles near a given Z level."""
    verts = set()
    for _, v0, v1, v2 in triangles:
        for v in [v0, v1, v2]:
            if abs(v[2] - z_target) < tol:
                verts.add((round(v[0], 2), round(v[1], 2)))
    return sorted(verts)


def main():
    stl_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'constant', 'triSurface')

    files = {
        'Current tower': os.path.join(stl_dir, 'citicorp_tower.stl'),
        'Current stilts': os.path.join(stl_dir, 'citicorp_stilts.stl'),
        'LOD2 tower': os.path.join(stl_dir, 'citicorp_tower_lod2.stl'),
        'LOD2 complex': os.path.join(stl_dir, 'citicorp_complex_lod2.stl'),
    }

    print("=" * 70)
    print("STL Geometry Comparison: Current Simplified vs NYC LOD2")
    print("=" * 70)

    data = {}
    for name, path in files.items():
        if not os.path.exists(path):
            print(f"\n  {name}: FILE NOT FOUND ({path})")
            continue
        tris = parse_ascii_stl(path)
        bb = bbox(tris)
        z_levels = unique_z_levels(tris)
        data[name] = {'tris': tris, 'bbox': bb, 'z_levels': z_levels}

        print(f"\n--- {name} ---")
        print(f"  Triangles: {len(tris)}")
        print(f"  X: {bb['xmin']:.2f} to {bb['xmax']:.2f} m  (span {bb['xmax']-bb['xmin']:.2f} m)")
        print(f"  Y: {bb['ymin']:.2f} to {bb['ymax']:.2f} m  (span {bb['ymax']-bb['ymin']:.2f} m)")
        print(f"  Z: {bb['zmin']:.2f} to {bb['zmax']:.2f} m  (height {bb['zmax']-bb['zmin']:.2f} m)")
        print(f"  Z levels: {len(z_levels)}")
        for z in z_levels:
            print(f"    {z:.1f} m", end="")
            if z < 1:
                print(" (ground)")
            elif abs(z - 34.75) < 1:
                print(" (stilt top)")
            elif abs(z - 248) < 5:
                print(" (main roof)")
            elif z > 270:
                print(" (crown)")
            else:
                print()

    # Detailed comparison
    if 'Current tower' in data and 'LOD2 tower' in data:
        print("\n" + "=" * 70)
        print("Detailed Tower Comparison")
        print("=" * 70)

        curr = data['Current tower']
        lod2 = data['LOD2 tower']

        print("\n  Property              Current         LOD2          Difference")
        print("  " + "-" * 65)

        # Width comparisons
        curr_w = curr['bbox']['xmax'] - curr['bbox']['xmin']
        lod2_w = lod2['bbox']['xmax'] - lod2['bbox']['xmin']
        print(f"  X span (m)            {curr_w:>8.2f}       {lod2_w:>8.2f}       {lod2_w-curr_w:>+8.2f}")

        curr_d = curr['bbox']['ymax'] - curr['bbox']['ymin']
        lod2_d = lod2['bbox']['ymax'] - lod2['bbox']['ymin']
        print(f"  Y span (m)            {curr_d:>8.2f}       {lod2_d:>8.2f}       {lod2_d-curr_d:>+8.2f}")

        curr_h = curr['bbox']['zmax'] - curr['bbox']['zmin']
        lod2_h = lod2['bbox']['zmax'] - lod2['bbox']['zmin']
        print(f"  Height (m)            {curr_h:>8.2f}       {lod2_h:>8.2f}       {lod2_h-curr_h:>+8.2f}")

        print(f"  Triangle count        {len(curr['tris']):>8d}       {len(lod2['tris']):>8d}       {len(lod2['tris'])-len(curr['tris']):>+8d}")
        print(f"  Z levels              {len(curr['z_levels']):>8d}       {len(lod2['z_levels']):>8d}       {len(lod2['z_levels'])-len(curr['z_levels']):>+8d}")

        # Rotation check: current tower is a 45-degree rotated square
        # Vertices should be at (+-d, 0) and (0, +-d)
        print("\n  Shape analysis:")
        curr_fp = footprint_at_z(curr['tris'], 100, 200)  # wide Z range to get wall vertices
        lod2_fp = footprint_at_z(lod2['tris'], 150, 100)  # mid-height for consistent cross-section

        # Check for diamond pattern (45 deg rotation)
        if curr_fp:
            # Find extreme X and Y vertices
            max_x_v = max(curr_fp, key=lambda v: v[0])
            min_x_v = min(curr_fp, key=lambda v: v[0])
            max_y_v = max(curr_fp, key=lambda v: v[1])
            min_y_v = min(curr_fp, key=lambda v: v[1])
            print(f"    Current: diamond with corners at")
            print(f"      (+X): ({max_x_v[0]:>7.2f}, {max_x_v[1]:>7.2f})")
            print(f"      (-X): ({min_x_v[0]:>7.2f}, {min_x_v[1]:>7.2f})")
            print(f"      (+Y): ({max_y_v[0]:>7.2f}, {max_y_v[1]:>7.2f})")
            print(f"      (-Y): ({min_y_v[0]:>7.2f}, {min_y_v[1]:>7.2f})")

        if lod2_fp:
            max_x_v = max(lod2_fp, key=lambda v: v[0])
            min_x_v = min(lod2_fp, key=lambda v: v[0])
            max_y_v = max(lod2_fp, key=lambda v: v[1])
            min_y_v = min(lod2_fp, key=lambda v: v[1])
            print(f"    LOD2:    {len(lod2_fp)} unique XY vertices at mid-height")
            print(f"      (+X): ({max_x_v[0]:>7.2f}, {max_x_v[1]:>7.2f})")
            print(f"      (-X): ({min_x_v[0]:>7.2f}, {min_x_v[1]:>7.2f})")
            print(f"      (+Y): ({max_y_v[0]:>7.2f}, {max_y_v[1]:>7.2f})")
            print(f"      (-Y): ({min_y_v[0]:>7.2f}, {min_y_v[1]:>7.2f})")

    # Stilt comparison
    if 'Current stilts' in data:
        print("\n" + "=" * 70)
        print("Stilt Analysis")
        print("=" * 70)
        stilts = data['Current stilts']
        print(f"  Current stilts: {len(stilts['tris'])} triangles, 4 columns")
        print(f"  NOTE: LOD2 CityGML does NOT include stilts (under building)")
        print(f"  Current stilt geometry should be improved separately:")
        print(f"    - Rotate 45 deg to match tower orientation")
        print(f"    - Place at face midpoints (not corners)")
        print(f"    - Current corners at ({stilts['bbox']['xmax']:.2f}, {stilts['bbox']['ymax']:.2f})")
        print(f"      protrude {stilts['bbox']['xmax'] - 33.835:.2f} m past tower face")

    print()


if __name__ == '__main__':
    main()
