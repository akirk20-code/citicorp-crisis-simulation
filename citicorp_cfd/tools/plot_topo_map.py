#!/usr/bin/env python3
"""
plot_topo_map.py
================
Generate a topographic height map (urban canopy map) from city STL files.
Reads surroundings.stl + citicorp_tower.stl, projects building footprints
onto a 2D grid, and plots maximum building height per cell as a color map.

Usage:
    python plot_topo_map.py [--stl-dir PATH] [--resolution 5.0] [--no-show]

If --stl-dir is not specified, looks in ./constant/triSurface/ (default CFD output).
"""

import os
import sys
import struct
import math
import argparse
import numpy as np

# ── optional matplotlib ────────────────────────────────────────────────────────
try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.collections import PatchCollection
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ─── STL READER ───────────────────────────────────────────────────────────────

def read_stl(path):
    """
    Read ASCII or binary STL file.
    Returns list of triangles: each is ((x0,y0,z0),(x1,y1,z1),(x2,y2,z2))
    """
    if not os.path.exists(path):
        return []
    tris = []
    with open(path, 'rb') as f:
        header = f.read(5)
    if header.startswith(b'solid'):
        # Try ASCII first
        try:
            tris = _read_stl_ascii(path)
            if tris:
                return tris
        except Exception:
            pass
    return _read_stl_binary(path)


def _read_stl_ascii(path):
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


def _read_stl_binary(path):
    tris = []
    with open(path, 'rb') as f:
        f.read(80)  # header
        count = struct.unpack('<I', f.read(4))[0]
        for _ in range(count):
            f.read(12)  # normal
            v = []
            for __ in range(3):
                x, y, z = struct.unpack('<fff', f.read(12))
                v.append((x, y, z))
            tris.append(tuple(v))
            f.read(2)  # attribute byte count
    return tris

# ─── HEIGHT GRID ─────────────────────────────────────────────────────────────

def build_height_grid(tris_list, resolution=5.0, padding=50.0):
    """
    Given list of triangle arrays (combined from multiple STLs),
    build a 2D grid of maximum building heights.

    resolution : grid cell size in meters
    padding    : extra margin beyond the data extent, in meters
    """
    all_tris = []
    for tris in tris_list:
        all_tris.extend(tris)

    if not all_tris:
        print("  No triangles found.")
        return None

    # Find extent
    xs = [v[0] for tri in all_tris for v in tri]
    ys = [v[1] for tri in all_tris for v in tri]
    zs = [v[2] for tri in all_tris for v in tri]

    xmin, xmax = min(xs) - padding, max(xs) + padding
    ymin, ymax = min(ys) - padding, max(ys) + padding
    zmax_global = max(zs)

    print(f"  Domain extent:  X [{xmin:.0f}, {xmax:.0f}] m,  Y [{ymin:.0f}, {ymax:.0f}] m")
    print(f"  Max height:     {zmax_global:.1f} m")
    print(f"  Grid resolution: {resolution} m")

    nx = int((xmax - xmin) / resolution) + 1
    ny = int((ymax - ymin) / resolution) + 1
    print(f"  Grid size:      {nx} x {ny} = {nx*ny:,} cells")

    grid = np.zeros((ny, nx), dtype=np.float32)

    # Rasterize each triangle onto the grid
    # For each triangle, find bounding box, then check which cells are inside
    for tri in all_tris:
        x0, y0, z0 = tri[0]
        x1, y1, z1 = tri[1]
        x2, y2, z2 = tri[2]

        # Triangle bounding box in grid coords
        xi0 = max(0, int((min(x0, x1, x2) - xmin) / resolution))
        xi1 = min(nx - 1, int((max(x0, x1, x2) - xmin) / resolution) + 1)
        yi0 = max(0, int((min(y0, y1, y2) - ymin) / resolution))
        yi1 = min(ny - 1, int((max(y0, y1, y2) - ymin) / resolution) + 1)

        zmax_tri = max(z0, z1, z2)
        if zmax_tri < 1.0:  # skip ground-level triangles
            continue

        for yi in range(yi0, yi1 + 1):
            for xi in range(xi0, xi1 + 1):
                px = xmin + xi * resolution
                py = ymin + yi * resolution
                if _point_in_triangle_2d(px, py, x0, y0, x1, y1, x2, y2):
                    if zmax_tri > grid[yi, xi]:
                        grid[yi, xi] = zmax_tri

    # Build coordinate arrays
    x_arr = np.linspace(xmin, xmax, nx)
    y_arr = np.linspace(ymin, ymax, ny)

    return grid, x_arr, y_arr, xmin, xmax, ymin, ymax


def _point_in_triangle_2d(px, py, x0, y0, x1, y1, x2, y2):
    """Barycentric test for point in 2D triangle."""
    d1 = (px - x1) * (y0 - y1) - (x0 - x1) * (py - y1)
    d2 = (px - x2) * (y1 - y2) - (x1 - x2) * (py - y2)
    d3 = (px - x0) * (y2 - y0) - (x2 - x0) * (py - y0)
    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (has_neg and has_pos)

# ─── PLOT ────────────────────────────────────────────────────────────────────

def plot_topo_map(grid, x_arr, y_arr, title, out_path):
    """Generate and save the topographic height map."""
    fig, ax = plt.subplots(figsize=(14, 12))

    # Dark background for contrast
    fig.patch.set_facecolor('#0d1117')
    ax.set_facecolor('#161b22')

    # Mask zero cells (no building)
    masked = np.ma.masked_where(grid < 1.0, grid)

    # Draw ground as dark
    ax.imshow(np.ones_like(grid) * 0.08,
              extent=[x_arr[0], x_arr[-1], y_arr[0], y_arr[-1]],
              origin='lower', aspect='equal', cmap='Greys', vmin=0, vmax=1,
              zorder=0)

    # Draw building heights — high contrast colormap
    vmax = min(grid.max(), 320.0)
    im = ax.imshow(masked,
                   extent=[x_arr[0], x_arr[-1], y_arr[0], y_arr[-1]],
                   origin='lower', aspect='equal',
                   cmap='inferno', vmin=5.0, vmax=vmax,
                   zorder=1)

    # Contour lines at key height intervals
    try:
        levels = [30, 60, 100, 150, 200, 250]
        levels = [l for l in levels if l < vmax]
        X, Y = np.meshgrid(x_arr, y_arr)
        cs = ax.contour(X, Y, grid, levels=levels,
                        colors='white', linewidths=0.3, alpha=0.35, zorder=2)
        ax.clabel(cs, inline=True, fontsize=7, fmt='%dm', colors='#aaaaaa')
    except Exception:
        pass

    # Mark Citicorp at origin
    ax.plot(0, 0, 'c^', markersize=16, zorder=5, label='Citicorp Center (278.9 m)',
            markeredgecolor='white', markeredgewidth=1.0)
    ax.annotate('Citicorp\n278.9 m', xy=(0, 0), xytext=(50, 50),
                color='cyan', fontsize=11, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='cyan', lw=2), zorder=6)

    # ── Wind direction arrow ───────────────────────────────────────────────
    # CFD wind is +X (west to east). Place arrow in upper-left area.
    arrow_y = y_arr[-1] - 100
    arrow_x = x_arr[0] + 80
    arrow_len = 180
    ax.annotate('', xy=(arrow_x + arrow_len, arrow_y), xytext=(arrow_x, arrow_y),
                arrowprops=dict(arrowstyle='->', color='#00e5ff', lw=3.0,
                                mutation_scale=25),
                zorder=7)
    ax.text(arrow_x + arrow_len / 2, arrow_y + 35,
            'Wind (44.7 m/s)', ha='center', va='bottom',
            color='#00e5ff', fontsize=12, fontweight='bold', zorder=7)
    ax.text(arrow_x - 10, arrow_y, 'W', ha='right', va='center',
            color='#00e5ff', fontsize=11, fontweight='bold', zorder=7)
    ax.text(arrow_x + arrow_len + 10, arrow_y, 'E', ha='left', va='center',
            color='#00e5ff', fontsize=11, fontweight='bold', zorder=7)

    # ── Cardinal direction labels on axes ──────────────────────────────────
    ax.set_xlabel('West  <--  East-West (m)  -->  East', fontsize=13, color='white')
    ax.set_ylabel('South  <--  North-South (m)  -->  North', fontsize=13, color='white')
    ax.tick_params(colors='#cccccc', labelsize=11)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.82)
    cbar.set_label('Building Height (m)', fontsize=13, color='white')
    cbar.ax.tick_params(colors='#cccccc', labelsize=10)

    # Scale bar (200m)
    xbar = x_arr[-1] - 280
    ybar = y_arr[0] + 40
    ax.plot([xbar, xbar + 200], [ybar, ybar], 'w-', lw=3, zorder=6)
    ax.plot([xbar, xbar], [ybar - 8, ybar + 8], 'w-', lw=1.5, zorder=6)
    ax.plot([xbar + 200, xbar + 200], [ybar - 8, ybar + 8], 'w-', lw=1.5, zorder=6)
    ax.text(xbar + 100, ybar + 18, '200 m', ha='center', color='white',
            fontsize=10, fontweight='bold', zorder=6)

    # North arrow (compass) in lower-left
    nx_pos = x_arr[0] + 80
    ny_pos = y_arr[0] + 120
    ax.annotate('', xy=(nx_pos, ny_pos + 60), xytext=(nx_pos, ny_pos),
                arrowprops=dict(arrowstyle='->', color='white', lw=2.0,
                                mutation_scale=20),
                zorder=7)
    ax.text(nx_pos, ny_pos + 70, 'N', ha='center', va='bottom',
            color='white', fontsize=13, fontweight='bold', zorder=7)

    ax.set_title(title, fontsize=15, fontweight='bold', pad=16, color='white')
    ax.legend(loc='upper right', fontsize=11, facecolor='#161b22', edgecolor='#444',
              labelcolor='white')
    ax.grid(False)

    for spine in ax.spines.values():
        spine.set_color('#444444')

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='#0d1117')
    print(f"  Saved: {out_path}")
    return fig

# ─── MAIN ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Topographic height map from CFD STL geometry')
    parser.add_argument('--stl-dir', default=None,
                        help='Directory containing STL files (default: ./constant/triSurface)')
    parser.add_argument('--resolution', type=float, default=5.0,
                        help='Grid cell size in meters (default: 5.0)')
    parser.add_argument('--no-show', action='store_true',
                        help='Save figure but do not open display window')
    args = parser.parse_args()

    if not HAS_MPL:
        print("ERROR: matplotlib not installed. Run: pip install matplotlib numpy")
        sys.exit(1)

    # Locate STL directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if args.stl_dir:
        stl_dir = args.stl_dir
    else:
        stl_dir = os.path.join(script_dir, 'constant', 'triSurface')

    print(f"\nTopographic Height Map Generator")
    print("=" * 55)
    print(f"STL directory: {stl_dir}")
    print(f"Resolution:    {args.resolution} m")

    # Find all STL files
    stl_files = {}
    for name in ['citicorp_tower', 'citicorp_stilts', 'surroundings']:
        path = os.path.join(stl_dir, f'{name}.stl')
        if os.path.exists(path):
            stl_files[name] = path
            size_kb = os.path.getsize(path) / 1024
            print(f"  Found: {name}.stl ({size_kb:.0f} KB)")
        else:
            print(f"  Missing: {name}.stl")

    if not stl_files:
        print("\nERROR: No STL files found. Run a geometry script first:")
        print(f"  python generate_stl_nyc3d.py")
        sys.exit(1)

    # Load triangles
    print("\nLoading triangles...")
    tris_all = []
    n_total = 0
    for name, path in stl_files.items():
        tris = read_stl(path)
        n_total += len(tris)
        tris_all.append(tris)
        print(f"  {name}: {len(tris):,} triangles")

    print(f"  Total: {n_total:,} triangles")

    # Build height grid
    print("\nBuilding height grid...")
    result = build_height_grid(tris_all, resolution=args.resolution)
    if result is None:
        sys.exit(1)
    grid, x_arr, y_arr, xmin, xmax, ymin, ymax = result

    # Statistics
    nonzero = grid[grid > 1.0]
    if len(nonzero):
        print(f"\nHeight statistics (buildings only):")
        print(f"  Min  : {nonzero.min():.1f} m")
        print(f"  Mean : {nonzero.mean():.1f} m")
        print(f"  Max  : {nonzero.max():.1f} m")
        print(f"  p75  : {np.percentile(nonzero, 75):.1f} m")
        print(f"  p95  : {np.percentile(nonzero, 95):.1f} m")
        footprint_frac = len(nonzero) / grid.size
        print(f"  Building coverage: {footprint_frac*100:.1f}% of domain")

    # Detect which dataset is active (by triangle count signature)
    n_surr = len(tris_all[-1]) if tris_all else 0
    if n_surr > 50000:
        label = "1978 Historical Skyline"
    elif n_surr > 25000:
        label = "2024 Present-Day Skyline"
    elif n_surr == 0:
        label = "Tower Only (no surroundings)"
    else:
        label = "City Model"

    title = f"Citicorp CFD — Urban Canopy Height Map\n{label}  |  {n_total:,} triangles  |  {args.resolution}m grid"

    # Save map
    out_path = os.path.join(script_dir, 'topo_map.png')
    print("\nRendering map...")
    fig = plot_topo_map(grid, x_arr, y_arr, title, out_path)

    # Also save CSV of height grid for further analysis
    csv_path = os.path.join(script_dir, 'topo_grid.csv')
    print(f"\nSaving height grid to CSV: {csv_path}")
    # Save compact: x, y, height for nonzero cells only
    with open(csv_path, 'w') as f:
        f.write("x_m,y_m,height_m\n")
        for yi in range(grid.shape[0]):
            for xi in range(grid.shape[1]):
                h = grid[yi, xi]
                if h > 1.0:
                    f.write(f"{x_arr[xi]:.1f},{y_arr[yi]:.1f},{h:.1f}\n")
    print(f"  Saved: {csv_path}")

    if not args.no_show:
        plt.show()

    print("\nDone.")


if __name__ == '__main__':
    main()
