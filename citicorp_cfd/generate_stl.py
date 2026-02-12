#!/usr/bin/env python3
"""Generate STL geometry for Citicorp Center CFD simulation.

Creates watertight box STL files for:
  1. Citicorp tower (elevated box on stilts)
  2. 4 stilt columns at face midpoints
  3. ~18 surrounding Midtown Manhattan buildings

All dimensions in meters. Building center at origin.
No external dependencies — standard library only.
"""

import os
import math

# ============================================================
# GEOMETRY PARAMETERS (meters)
# ============================================================
TOWER_W = 47.85       # 157 ft — square plan
TOWER_H_TOP = 278.9   # 915 ft — roof height
STILT_H = 34.75       # 114 ft — stilt height (tower bottom)
STILT_W = 7.32        # 24 ft — stilt cross-section
HALF_T = TOWER_W / 2  # 23.925 m
HALF_S = STILT_W / 2  # 3.66 m

# Stilt centers (at face midpoints)
STILTS = [
    (0,        -HALF_T),  # South
    (HALF_T,    0),       # East
    (0,         HALF_T),  # North
    (-HALF_T,   0),       # West
]

# Surrounding buildings: (name, cx, cy, wx, wy, height_m)
# Approximate Midtown East layout around Citicorp (601 Lexington Ave)
# Avenues ~230m apart (E-W), streets ~80m apart (N-S)
SURROUNDINGS = [
    # Immediate neighbors
    ("399_Park_Ave",       -200,  -50,  60, 45, 162),
    ("280_Park_Ave",       -200,   50,  55, 40, 135),
    ("780_Third_Ave",       200,  -25,  50, 50, 175),
    ("919_Third_Ave",       200,   60,  55, 45, 182),
    ("885_Third_Ave",       200,  180,  35, 35, 138),
    ("599_Lexington",      -100,   80,  45, 40, 199),
    ("731_Lexington",      -100, -100,  50, 45, 180),
    ("135_E_54th",            0,  100,  40, 35, 150),
    ("153_E_53rd",           80,  -80,  35, 30, 120),
    # Second ring
    ("Park_Ave_Tower_S",   -400,  -50,  55, 40, 150),
    ("Park_Ave_Tower_N",   -400,   60,  50, 35, 120),
    ("Second_Ave_S",        400,  -60,  50, 40, 130),
    ("Second_Ave_N",        400,   70,  45, 35, 110),
    ("E_55th_Lex",         -110,  170,  40, 30, 100),
    ("E_52nd_Lex",         -110, -180,  45, 35,  90),
    ("E_55th_Third",        110,  170,  35, 30,  85),
    ("E_52nd_Third",        110, -170,  40, 35,  95),
    # Infill
    ("E_53rd_mid",          -50, -100,  30, 25,  70),
]

# ============================================================
# STL WRITER
# ============================================================

def box_triangles(xmin, ymin, zmin, xmax, ymax, zmax):
    """Return 12 triangles for a watertight box with outward normals.

    Each triangle is ((v1, v2, v3), (nx, ny, nz)).
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
    # (indices, outward normal)
    faces = [
        ((0, 2, 1), (0, 0, -1)),   # bottom -z
        ((0, 3, 2), (0, 0, -1)),
        ((4, 5, 6), (0, 0, 1)),    # top +z
        ((4, 6, 7), (0, 0, 1)),
        ((0, 1, 5), (0, -1, 0)),   # front -y
        ((0, 5, 4), (0, -1, 0)),
        ((2, 3, 7), (0, 1, 0)),    # back +y
        ((2, 7, 6), (0, 1, 0)),
        ((0, 4, 7), (-1, 0, 0)),   # left -x
        ((0, 7, 3), (-1, 0, 0)),
        ((1, 2, 6), (1, 0, 0)),    # right +x
        ((1, 6, 5), (1, 0, 0)),
    ]
    return [(tuple(v[i] for i in idx), n) for idx, n in faces]


def write_stl(filepath, solid_name, triangles):
    """Write triangles to an ASCII STL file."""
    with open(filepath, 'w') as f:
        f.write(f"solid {solid_name}\n")
        for verts, normal in triangles:
            f.write(f"  facet normal {normal[0]:.6f} {normal[1]:.6f} {normal[2]:.6f}\n")
            f.write(f"    outer loop\n")
            for vx, vy, vz in verts:
                f.write(f"      vertex {vx:.4f} {vy:.4f} {vz:.4f}\n")
            f.write(f"    endloop\n")
            f.write(f"  endfacet\n")
        f.write(f"endsolid {solid_name}\n")
    return len(triangles)


# ============================================================
# MAIN
# ============================================================

def main():
    # Output directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(script_dir, "constant", "triSurface")
    os.makedirs(out_dir, exist_ok=True)

    total_tris = 0

    # --- Citicorp Tower ---
    tower_tris = box_triangles(
        -HALF_T, -HALF_T, STILT_H,
         HALF_T,  HALF_T, TOWER_H_TOP
    )
    n = write_stl(os.path.join(out_dir, "citicorp_tower.stl"),
                  "citicorp_tower", tower_tris)
    total_tris += n
    print(f"  citicorp_tower.stl: {n} triangles")
    print(f"    Bounds: ({-HALF_T:.1f}, {-HALF_T:.1f}, {STILT_H:.1f}) "
          f"to ({HALF_T:.1f}, {HALF_T:.1f}, {TOWER_H_TOP:.1f})")

    # --- Stilt Columns ---
    stilt_tris = []
    for i, (cx, cy) in enumerate(STILTS):
        tris = box_triangles(
            cx - HALF_S, cy - HALF_S, 0,
            cx + HALF_S, cy + HALF_S, STILT_H
        )
        stilt_tris.extend(tris)
    n = write_stl(os.path.join(out_dir, "citicorp_stilts.stl"),
                  "citicorp_stilts", stilt_tris)
    total_tris += n
    print(f"  citicorp_stilts.stl: {n} triangles (4 stilts)")

    # --- Surrounding Buildings ---
    surr_tris = []
    for name, cx, cy, wx, wy, h in SURROUNDINGS:
        tris = box_triangles(
            cx - wx/2, cy - wy/2, 0,
            cx + wx/2, cy + wy/2, h
        )
        surr_tris.extend(tris)
    n = write_stl(os.path.join(out_dir, "surroundings.stl"),
                  "surroundings", surr_tris)
    total_tris += n
    print(f"  surroundings.stl: {n} triangles ({len(SURROUNDINGS)} buildings)")

    print(f"\nTotal: {total_tris} triangles in {out_dir}")
    print("STL generation complete.")


if __name__ == "__main__":
    main()
