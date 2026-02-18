#!/usr/bin/env python3
"""
plot_wind_elevation.py
======================
Generate a 2D elevation view from the wind direction (west, looking east).
Projects all building geometry onto the Y-Z plane, showing silhouettes
of the city with the Citicorp tower rising above surrounding buildings.

Usage:
    python plot_wind_elevation.py [--stl-dir PATH] [--resolution 2.0] [--no-show]
"""

import os, sys, math, struct, argparse
import numpy as np

try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ─── STL READER ───────────────────────────────────────────────────────────────

def read_stl(path):
    """Read ASCII or binary STL. Returns list of ((x,y,z),(x,y,z),(x,y,z))."""
    if not os.path.exists(path):
        return []
    with open(path, 'rb') as f:
        header = f.read(5)
    if header.startswith(b'solid'):
        try:
            tris = _read_ascii(path)
            if tris:
                return tris
        except Exception:
            pass
    return _read_binary(path)

def _read_ascii(path):
    tris = []
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        verts = []
        for line in f:
            line = line.strip()
            if line.startswith('vertex'):
                parts = line.split()
                verts.append((float(parts[1]), float(parts[2]), float(parts[3])))
                if len(verts) == 3:
                    tris.append(tuple(verts))
                    verts = []
    return tris

def _read_binary(path):
    tris = []
    with open(path, 'rb') as f:
        f.read(80)
        count = struct.unpack('<I', f.read(4))[0]
        for _ in range(count):
            f.read(12)
            v = []
            for __ in range(3):
                x, y, z = struct.unpack('<fff', f.read(12))
                v.append((x, y, z))
            tris.append(tuple(v))
            f.read(2)
    return tris

# ─── RASTERIZE TRIANGLES ONTO Y-Z PLANE ──────────────────────────────────────

def rasterize_yz(tris, y_arr, z_arr, grid, value=1.0, x_values=None):
    """
    Project triangles onto Y-Z plane and fill grid cells.

    tris     : list of ((x0,y0,z0), (x1,y1,z1), (x2,y2,z2))
    y_arr    : 1D array of Y bin edges
    z_arr    : 1D array of Z bin edges
    grid     : 2D array [nz, ny] to fill
    value    : fill value (or 'x_distance' to use per-triangle X centroid)
    x_values : if not None, fill grid with X-centroid for depth ordering
    """
    dy = y_arr[1] - y_arr[0]
    dz = z_arr[1] - z_arr[0]
    ny = len(y_arr)
    nz = len(z_arr)

    for tri in tris:
        # Project to Y-Z: use (y, z) of each vertex
        y0, z0 = tri[0][1], tri[0][2]
        y1, z1 = tri[1][1], tri[1][2]
        y2, z2 = tri[2][1], tri[2][2]

        # Bounding box in grid coords
        ymin_t = min(y0, y1, y2)
        ymax_t = max(y0, y1, y2)
        zmin_t = min(z0, z1, z2)
        zmax_t = max(z0, z1, z2)

        # Skip degenerate or ground-only triangles
        if zmax_t < 0.5:
            continue

        yi0 = max(0, int((ymin_t - y_arr[0]) / dy))
        yi1 = min(ny - 1, int((ymax_t - y_arr[0]) / dy) + 1)
        zi0 = max(0, int((zmin_t - z_arr[0]) / dz))
        zi1 = min(nz - 1, int((zmax_t - z_arr[0]) / dz) + 1)

        # X centroid for depth value
        if x_values is not None:
            x_cent = (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0

        for zi in range(zi0, zi1 + 1):
            for yi in range(yi0, yi1 + 1):
                py = y_arr[0] + yi * dy
                pz = z_arr[0] + zi * dz
                if _pt_in_tri(py, pz, y0, z0, y1, z1, y2, z2):
                    if x_values is not None:
                        # Store X distance (for depth coloring)
                        grid[zi, yi] = x_cent
                    else:
                        grid[zi, yi] = value


def _pt_in_tri(px, py, x0, y0, x1, y1, x2, y2):
    d1 = (px - x1) * (y0 - y1) - (x0 - x1) * (py - y1)
    d2 = (px - x2) * (y1 - y2) - (x1 - x2) * (py - y2)
    d3 = (px - x0) * (y2 - y0) - (x2 - x0) * (py - y0)
    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (has_neg and has_pos)

# ─── MAIN ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='2D elevation view from wind direction (west-facing)')
    parser.add_argument('--stl-dir', default=None)
    parser.add_argument('--resolution', type=float, default=2.0,
                        help='Grid cell size in meters (default: 2.0)')
    parser.add_argument('--no-show', action='store_true')
    args = parser.parse_args()

    if not HAS_MPL:
        print("ERROR: pip install matplotlib numpy")
        sys.exit(1)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    stl_dir = args.stl_dir or os.path.join(script_dir, 'constant', 'triSurface')

    print("\n2D Wind-Direction Elevation View")
    print("=" * 55)
    print(f"STL directory: {stl_dir}")
    print(f"View: looking from west (upwind) toward east")
    print(f"Resolution: {args.resolution} m")

    # Load STLs
    tower_tris = read_stl(os.path.join(stl_dir, 'citicorp_tower.stl'))
    stilt_tris = read_stl(os.path.join(stl_dir, 'citicorp_stilts.stl'))
    surr_tris  = read_stl(os.path.join(stl_dir, 'surroundings.stl'))

    print(f"\n  Tower:        {len(tower_tris):,} triangles")
    print(f"  Stilts:       {len(stilt_tris):,} triangles")
    print(f"  Surroundings: {len(surr_tris):,} triangles")

    all_tris = tower_tris + stilt_tris + surr_tris
    if not all_tris:
        print("No triangles found!")
        sys.exit(1)

    # Find Y-Z extent
    ys = [v[1] for t in all_tris for v in t]
    zs = [v[2] for t in all_tris for v in t]
    xs = [v[0] for t in all_tris for v in t]

    pad = 30.0
    ymin, ymax = min(ys) - pad, max(ys) + pad
    zmin, zmax = 0.0, max(zs) + 20.0
    xmin, xmax = min(xs), max(xs)

    res = args.resolution
    ny = int((ymax - ymin) / res) + 1
    nz = int((zmax - zmin) / res) + 1
    y_arr = np.linspace(ymin, ymax, ny)
    z_arr = np.linspace(zmin, zmax, nz)

    print(f"\n  Y range: [{ymin:.0f}, {ymax:.0f}] m  ({ny} bins)")
    print(f"  Z range: [{zmin:.0f}, {zmax:.0f}] m  ({nz} bins)")
    print(f"  X depth: [{xmin:.0f}, {xmax:.0f}] m")

    # ── Painters algorithm: all triangles sorted back-to-front ─────────
    # Closer buildings (smaller X = more west = closer to viewer) drawn last,
    # so they naturally occlude Citicorp and buildings behind them.
    print("\nBuilding tagged triangle list...")
    citi_set = set(id(t) for t in tower_tris + stilt_tris)
    tagged = []  # (tri, x_centroid, is_citicorp)
    for tri in surr_tris:
        x_c = (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0
        tagged.append((tri, x_c, False))
    for tri in (tower_tris + stilt_tris):
        x_c = (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0
        tagged.append((tri, x_c, True))

    # Sort: largest X first (furthest from viewer = drawn first, overwritten by closer)
    tagged.sort(key=lambda t: t[1], reverse=True)
    print(f"  {len(tagged):,} triangles sorted by depth")

    # Two grids: category (0=surr, 1=citi) and depth value
    cat_grid   = np.full((nz, ny), np.nan, dtype=np.float32)  # NaN=empty
    depth_grid = np.full((nz, ny), np.nan, dtype=np.float32)

    x_range = xmax - xmin + 1e-6
    print("Rasterizing (back-to-front)...")

    for tri, x_cent, is_citi in tagged:
        y0, z0 = tri[0][1], tri[0][2]
        y1, z1 = tri[1][1], tri[1][2]
        y2, z2 = tri[2][1], tri[2][2]

        zmax_t = max(z0, z1, z2)
        if zmax_t < 0.5:
            continue

        ymin_t = min(y0, y1, y2)
        ymax_t = max(y0, y1, y2)
        zmin_t = min(z0, z1, z2)

        yi0 = max(0, int((ymin_t - ymin) / res))
        yi1 = min(ny - 1, int((ymax_t - ymin) / res) + 1)
        zi0 = max(0, int((zmin_t - zmin) / res))
        zi1 = min(nz - 1, int((zmax_t - zmin) / res) + 1)

        # Depth: 0 = far (east/downwind), 1 = close (west/upwind)
        x_norm = 1.0 - (x_cent - xmin) / x_range

        for zi in range(zi0, zi1 + 1):
            for yi in range(yi0, yi1 + 1):
                py = ymin + yi * res
                pz = zmin + zi * res
                if _pt_in_tri(py, pz, y0, z0, y1, z1, y2, z2):
                    cat_grid[zi, yi] = 1.0 if is_citi else 0.0
                    depth_grid[zi, yi] = x_norm

    # ── Plot ─────────────────────────────────────────────────────────────
    print("Rendering elevation view...")

    fig, ax = plt.subplots(figsize=(16, 7))
    fig.patch.set_facecolor('#0d1117')
    ax.set_facecolor('#0f1318')

    extent = [ymin, ymax, zmin, zmax]

    # Build RGBA image from composite grids
    # Surroundings: copper/brown shaded by depth
    # Citicorp: cyan
    rgba = np.zeros((nz, ny, 4), dtype=np.float32)

    # Surroundings pixels
    surr_mask = cat_grid == 0.0
    cmap_surr = plt.cm.copper
    if np.any(surr_mask):
        depth_surr = np.where(surr_mask, depth_grid, 0.0)
        surr_colors = cmap_surr(depth_surr)
        rgba[surr_mask] = surr_colors[surr_mask]

    # Citicorp pixels (overwrite surroundings where Citicorp is in front)
    citi_mask = cat_grid == 1.0
    if np.any(citi_mask):
        rgba[citi_mask] = [0.0, 0.90, 1.0, 1.0]  # cyan

    ax.imshow(rgba, extent=extent, origin='lower', aspect='auto',
              interpolation='nearest', zorder=1)

    # ── Annotations ──────────────────────────────────────────────────────

    # Ground line
    ax.axhline(y=0, color='#555555', linewidth=1.5, zorder=3)

    # Citicorp label with key heights
    ax.annotate('Citicorp Center\n278.9 m', xy=(0, 279), xytext=(0, 305),
                color='#00e5ff', fontsize=13, fontweight='bold',
                ha='center', va='bottom',
                arrowprops=dict(arrowstyle='->', color='#00e5ff', lw=2),
                zorder=5)

    # Stilt height line
    stilt_h = 34.75
    ax.plot([-50, 50], [stilt_h, stilt_h], '--', color='#00e5ff',
            alpha=0.4, lw=1, zorder=4)
    ax.text(55, stilt_h, f'Stilts: {stilt_h} m', color='#00e5ff',
            fontsize=9, alpha=0.7, va='center', zorder=4)

    # Height reference lines
    for h, label in [(100, '100 m'), (200, '200 m')]:
        if h < zmax:
            ax.axhline(y=h, color='#333333', linewidth=0.5, linestyle=':', zorder=0)
            ax.text(ymax - 10, h + 3, label, color='#666666', fontsize=8,
                    ha='right', zorder=0)

    # ── Wind indicator: into the page ──────────────────────────────────
    # Wind is +X (west to east). Viewer stands upwind (west) looking east.
    # Wind blows past the viewer INTO the scene = into the page.
    # Standard symbol: circle with X = vector going away from viewer.
    from matplotlib.patches import Circle
    wind_symbols_y = np.linspace(ymin + 120, ymax - 120, 5)
    wind_symbols_z = [zmax - 25] * 5
    for wy, wz in zip(wind_symbols_y, wind_symbols_z):
        r = 8.0
        circ = Circle((wy, wz), r, fill=False, edgecolor='#00e5ff',
                       linewidth=1.5, zorder=7)
        ax.add_patch(circ)
        # X inside the circle (vector into page)
        d = r * 0.6
        ax.plot([wy - d, wy + d], [wz - d, wz + d], color='#00e5ff',
                lw=1.5, zorder=7)
        ax.plot([wy - d, wy + d], [wz + d, wz - d], color='#00e5ff',
                lw=1.5, zorder=7)

    ax.text(0, zmax - 8,
            'Wind: 44.7 m/s  West -> East  (into page)',
            color='#00e5ff', fontsize=12, fontweight='bold',
            ha='center', va='top', zorder=7,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#0d1117',
                      edgecolor='#00e5ff', alpha=0.85))

    # Depth legend
    leg_y = ymax - 220
    leg_z = 12
    ax.text(leg_y, leg_z, 'Bright = closer (upwind)',
            color='#dda050', fontsize=9, ha='left', zorder=6)
    ax.text(leg_y, leg_z - 10, 'Dark = farther (downwind)',
            color='#4a2800', fontsize=9, ha='left', zorder=6)

    # Axis labels (viewer faces east from the west: left=north, right=south)
    ax.invert_xaxis()  # flip so +Y (North) is on left, -Y (South) on right
    ax.set_xlabel('North  <--  Position (m)  -->  South', fontsize=13, color='white')
    ax.set_ylabel('Height (m)', fontsize=13, color='white')
    ax.tick_params(colors='#cccccc', labelsize=11)

    # Detect dataset
    n_surr = len(surr_tris)
    if n_surr > 50000:
        ds_label = "1978 Historical Skyline"
    elif n_surr > 25000:
        ds_label = "2024 Present-Day Skyline"
    elif n_surr == 0:
        ds_label = "Tower Only"
    else:
        ds_label = "City Model"

    ax.set_title(
        f'Citicorp CFD — View from Upwind (Standing West, Looking East)\n'
        f'{ds_label}  |  {len(all_tris):,} triangles  |  {res}m grid',
        fontsize=14, fontweight='bold', color='white', pad=14)

    for spine in ax.spines.values():
        spine.set_color('#444444')

    ax.set_xlim(ymax, ymin)  # inverted: North (+Y) on left, South (-Y) on right
    ax.set_ylim(zmin, min(zmax, 330))

    plt.tight_layout()
    out_path = os.path.join(script_dir, 'wind_elevation.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='#0d1117')
    print(f"\n  Saved: {out_path}")

    if not args.no_show:
        plt.show()

    print("Done.")


if __name__ == '__main__':
    main()
