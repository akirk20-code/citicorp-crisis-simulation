#!/usr/bin/env python3
"""
analyze_exposure.py
===================
Quantify how much of the Citicorp tower is visible vs. sheltered by
surrounding buildings, from the upwind direction (270 = West).

Computes exposure ratio = visible_cyan_pixels / total_tower_pixels
for the CFD wind direction and optionally all compass bearings.

Usage:
    python analyze_exposure.py [--all-bearings]
"""

import os, sys, math, struct, argparse
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
except ImportError:
    print("ERROR: pip install matplotlib numpy"); sys.exit(1)

# ─── STL READER ───────────────────────────────────────────────────────────────

def read_stl(path):
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

def _pt_in_tri(px, py, x0, y0, x1, y1, x2, y2):
    d1 = (px - x1) * (y0 - y1) - (x0 - x1) * (py - y1)
    d2 = (px - x2) * (y1 - y2) - (x1 - x2) * (py - y2)
    d3 = (px - x0) * (y2 - y0) - (x2 - x0) * (py - y0)
    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (has_neg and has_pos)

# ─── RASTERIZE WITH PROJECTION ──────────────────────────────────────────────

def rasterize_projected(tris_tagged, bearing_deg, lat_min, lat_max,
                        z_min, z_max, res):
    """
    Project triangles onto a plane perpendicular to the viewing direction
    and rasterize with painters algorithm.

    tris_tagged: list of (tri, xc, yc, is_citicorp)
    Returns: cat_grid (nz, nlat) where 0.0=surr, 1.0=citi, NaN=empty
    """
    theta = math.radians(bearing_deg)
    sin_t, cos_t = math.sin(theta), math.cos(theta)

    n_lat = int((lat_max - lat_min) / res) + 1
    n_z   = int((z_max - z_min) / res) + 1
    cat_grid = np.full((n_z, n_lat), np.nan, dtype=np.float32)

    # Compute depth for sorting
    projected = []
    for tri, xc, yc, is_citi in tris_tagged:
        depth = xc * sin_t + yc * cos_t
        verts_2d = []
        for vx, vy, vz in tri:
            lat = vx * (-cos_t) + vy * sin_t
            verts_2d.append((lat, vz))
        projected.append((depth, verts_2d, is_citi))

    # Back-to-front
    projected.sort(key=lambda p: p[0], reverse=True)

    for depth, verts_2d, is_citi in projected:
        l0, z0 = verts_2d[0]
        l1, z1 = verts_2d[1]
        l2, z2 = verts_2d[2]

        if max(z0, z1, z2) < 0.5:
            continue

        li0 = max(0, int((min(l0, l1, l2) - lat_min) / res))
        li1 = min(n_lat - 1, int((max(l0, l1, l2) - lat_min) / res) + 1)
        zi0 = max(0, int((min(z0, z1, z2) - z_min) / res))
        zi1 = min(n_z - 1, int((max(z0, z1, z2) - z_min) / res) + 1)

        for zi in range(zi0, zi1 + 1):
            for li in range(li0, li1 + 1):
                pl = lat_min + li * res
                pz = z_min + zi * res
                if _pt_in_tri(pl, pz, l0, z0, l1, z1, l2, z2):
                    cat_grid[zi, li] = 1.0 if is_citi else 0.0

    return cat_grid

# ─── MAIN ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--all-bearings', action='store_true',
                        help='Compute exposure for all 24 compass bearings')
    parser.add_argument('--resolution', type=float, default=2.0)
    parser.add_argument('--stl-dir', default=None)
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    stl_dir = args.stl_dir or os.path.join(script_dir, 'constant', 'triSurface')

    print("\nCiticorp Tower Exposure Analysis")
    print("=" * 60)

    tower_tris = read_stl(os.path.join(stl_dir, 'citicorp_tower.stl'))
    stilt_tris = read_stl(os.path.join(stl_dir, 'citicorp_stilts.stl'))
    surr_tris  = read_stl(os.path.join(stl_dir, 'surroundings.stl'))

    citi_tris = tower_tris + stilt_tris
    print(f"  Tower+stilts: {len(citi_tris):,} triangles")
    print(f"  Surroundings: {len(surr_tris):,} triangles")

    # Pre-compute centroids
    tagged_citi_only = []
    for tri in citi_tris:
        xc = (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0
        yc = (tri[0][1] + tri[1][1] + tri[2][1]) / 3.0
        tagged_citi_only.append((tri, xc, yc, True))

    tagged_all = []
    for tri in surr_tris:
        xc = (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0
        yc = (tri[0][1] + tri[1][1] + tri[2][1]) / 3.0
        tagged_all.append((tri, xc, yc, False))
    tagged_all.extend(tagged_citi_only)

    # Extents
    all_tris = citi_tris + surr_tris
    xs = [v[0] for t in all_tris for v in t]
    ys = [v[1] for t in all_tris for v in t]
    zs = [v[2] for t in all_tris for v in t]
    max_r = max(max(abs(min(xs)), abs(max(xs))),
                max(abs(min(ys)), abs(max(ys))))
    pad = 30
    lat_ext = (-max_r - pad, max_r + pad)
    z_ext = (0.0, max(zs) + 20.0)
    res = args.resolution

    # ── Compute for wind direction (270 = West) ─────────────────────────
    bearing = 270
    print(f"\n  Computing for bearing {bearing} (West = CFD wind direction)...")
    print(f"  Resolution: {res} m")

    # Pass 1: Tower alone (no surroundings) = full silhouette
    print("  [1/2] Rasterizing tower alone (unsheltered)...")
    grid_alone = rasterize_projected(tagged_citi_only, bearing,
                                     lat_ext[0], lat_ext[1],
                                     z_ext[0], z_ext[1], res)
    total_pixels = np.sum(grid_alone == 1.0)

    # Pass 2: Tower + surroundings (with occlusion)
    print("  [2/2] Rasterizing with city (sheltered)...")
    grid_city = rasterize_projected(tagged_all, bearing,
                                    lat_ext[0], lat_ext[1],
                                    z_ext[0], z_ext[1], res)
    visible_pixels = np.sum(grid_city == 1.0)
    blocked_pixels = total_pixels - visible_pixels

    exposure = visible_pixels / total_pixels * 100 if total_pixels > 0 else 0
    sheltered = 100 - exposure

    pixel_area = res * res  # m^2 per pixel
    total_area = total_pixels * pixel_area
    visible_area = visible_pixels * pixel_area
    blocked_area = blocked_pixels * pixel_area

    print(f"\n{'='*60}")
    print(f"EXPOSURE ANALYSIS — Wind from {bearing} (West)")
    print(f"{'='*60}")
    print(f"  Tower silhouette (alone):    {total_pixels:,} pixels  ({total_area:,.0f} m2)")
    print(f"  Visible (with city):         {visible_pixels:,} pixels  ({visible_area:,.0f} m2)")
    print(f"  Blocked by upwind bldgs:     {blocked_pixels:,} pixels  ({blocked_area:,.0f} m2)")
    print(f"")
    print(f"  EXPOSURE RATIO:              {exposure:.1f}%")
    print(f"  SHELTERED:                   {sheltered:.1f}%")
    print(f"{'='*60}")

    # ── Height-band breakdown ────────────────────────────────────────────
    n_z = grid_alone.shape[0]
    n_lat = grid_alone.shape[1]
    z_arr = np.linspace(z_ext[0], z_ext[1], n_z)

    print(f"\n  Height-band exposure:")
    print(f"  {'Band':20s} {'Total px':>10s} {'Visible':>10s} {'Exposure':>10s}")
    print(f"  {'-'*52}")

    bands = [(0, 34.75, 'Stilts (0-35m)'),
             (34.75, 100, 'Lower tower (35-100m)'),
             (100, 180, 'Mid tower (100-180m)'),
             (180, 248, 'Upper tower (180-248m)'),
             (248, 280, 'Crown (248-280m)')]

    band_data = []
    for zlo, zhi, name in bands:
        zi_lo = max(0, int((zlo - z_ext[0]) / res))
        zi_hi = min(n_z - 1, int((zhi - z_ext[0]) / res))
        alone_band = np.sum(grid_alone[zi_lo:zi_hi+1, :] == 1.0)
        city_band  = np.sum(grid_city[zi_lo:zi_hi+1, :] == 1.0)
        exp_band = city_band / alone_band * 100 if alone_band > 0 else 0
        band_data.append((name, alone_band, city_band, exp_band))
        print(f"  {name:20s} {alone_band:10,} {city_band:10,} {exp_band:8.1f}%")

    # ── Side-by-side visualization ───────────────────────────────────────
    print("\nRendering comparison figure...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor('#0d1117')

    extent = [lat_ext[0], lat_ext[1], z_ext[0], z_ext[1]]

    # Panel 1: Tower alone
    ax = axes[0]
    ax.set_facecolor('#0f1318')
    rgba1 = np.zeros((n_z, n_lat, 4), dtype=np.float32)
    mask1 = grid_alone == 1.0
    rgba1[mask1] = [0.0, 0.9, 1.0, 1.0]
    ax.imshow(rgba1, extent=extent, origin='lower', aspect='auto',
              interpolation='nearest')
    ax.set_title(f'Unsheltered Tower\n{total_pixels:,} px ({total_area:,.0f} m2)',
                 fontsize=12, fontweight='bold', color='white', pad=8)
    ax.set_xlim(lat_ext[1], lat_ext[0])
    ax.set_ylim(z_ext[0], 310)
    ax.set_ylabel('Height (m)', fontsize=11, color='white')
    ax.tick_params(colors='#ccc', labelsize=9)
    for s in ax.spines.values(): s.set_color('#444')

    # Panel 2: Tower with city
    ax = axes[1]
    ax.set_facecolor('#0f1318')
    rgba2 = np.zeros((n_z, n_lat, 4), dtype=np.float32)
    # Surroundings in brown
    smask = grid_city == 0.0
    rgba2[smask] = [0.7, 0.5, 0.3, 1.0]
    # Citicorp in cyan
    cmask = grid_city == 1.0
    rgba2[cmask] = [0.0, 0.9, 1.0, 1.0]
    ax.imshow(rgba2, extent=extent, origin='lower', aspect='auto',
              interpolation='nearest')
    ax.set_title(f'With City (from West)\n{visible_pixels:,} px visible ({visible_area:,.0f} m2)',
                 fontsize=12, fontweight='bold', color='white', pad=8)
    ax.set_xlim(lat_ext[1], lat_ext[0])
    ax.set_ylim(z_ext[0], 310)
    ax.tick_params(colors='#ccc', labelsize=9)
    for s in ax.spines.values(): s.set_color('#444')

    # Panel 3: Blocked area highlighted
    ax = axes[2]
    ax.set_facecolor('#0f1318')
    rgba3 = np.zeros((n_z, n_lat, 4), dtype=np.float32)
    # Visible Citicorp = cyan
    rgba3[cmask] = [0.0, 0.9, 1.0, 1.0]
    # Blocked area = red (was tower alone but now hidden)
    blocked_mask = (grid_alone == 1.0) & (grid_city != 1.0)
    rgba3[blocked_mask] = [1.0, 0.2, 0.2, 1.0]
    ax.imshow(rgba3, extent=extent, origin='lower', aspect='auto',
              interpolation='nearest')
    ax.set_title(f'Exposure Map\nCyan={exposure:.1f}% visible, Red={sheltered:.1f}% blocked',
                 fontsize=12, fontweight='bold', color='white', pad=8)
    ax.set_xlim(lat_ext[1], lat_ext[0])
    ax.set_ylim(z_ext[0], 310)
    ax.tick_params(colors='#ccc', labelsize=9)
    for s in ax.spines.values(): s.set_color('#444')

    # Height band dividers on all panels
    for ax in axes:
        for zlo, zhi, name in bands:
            ax.axhline(y=zhi, color='#444', lw=0.5, ls=':', zorder=3)
        ax.axhline(y=0, color='#555', lw=1, zorder=3)

    fig.suptitle('Citicorp Tower Exposure — Wind from West (270)',
                 fontsize=14, fontweight='bold', color='#00e5ff', y=1.02)
    plt.tight_layout()

    out_path = os.path.join(script_dir, 'exposure_analysis.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='#0d1117')
    plt.close(fig)
    print(f"  Saved: {out_path}")

    # ── All bearings ─────────────────────────────────────────────────────
    if args.all_bearings:
        print(f"\nComputing exposure for all compass bearings...")
        print(f"  {'Bearing':>8s} {'Dir':>5s} {'Total':>8s} {'Visible':>8s} {'Exposure':>9s}")
        print(f"  {'-'*42}")

        results = []
        for b in range(0, 360, 15):
            g_alone = rasterize_projected(tagged_citi_only, b,
                                          lat_ext[0], lat_ext[1],
                                          z_ext[0], z_ext[1], res)
            g_city = rasterize_projected(tagged_all, b,
                                         lat_ext[0], lat_ext[1],
                                         z_ext[0], z_ext[1], res)
            tp = int(np.sum(g_alone == 1.0))
            vp = int(np.sum(g_city == 1.0))
            exp_pct = vp / tp * 100 if tp > 0 else 0
            results.append((b, tp, vp, exp_pct))

            compass = {0:'N',45:'NE',90:'E',135:'SE',180:'S',225:'SW',270:'W',315:'NW'}
            lbl = compass.get(b, '')
            print(f"  {b:>5d}    {lbl:>4s} {tp:>8,} {vp:>8,} {exp_pct:>8.1f}%")

        # Plot polar chart of exposure vs bearing
        fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(projection='polar'))
        fig.patch.set_facecolor('#0d1117')
        ax.set_facecolor('#0f1318')

        bearings_rad = [math.radians(90 - r[0]) for r in results]  # convert to math angle
        exposures = [r[3] for r in results]
        # Close the loop
        bearings_rad.append(bearings_rad[0])
        exposures.append(exposures[0])

        ax.plot(bearings_rad, exposures, 'o-', color='#00e5ff', lw=2, markersize=5)
        ax.fill(bearings_rad, exposures, alpha=0.15, color='#00e5ff')

        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_rlabel_position(45)
        ax.tick_params(colors='#ccc', labelsize=9)
        ax.set_rticks([20, 40, 60, 80, 100])
        ax.set_rlim(0, 100)
        ax.grid(color='#333', linewidth=0.5)

        # Mark wind direction
        wind_rad = math.radians(90 - 270)
        ax.annotate('', xy=(wind_rad, 95), xytext=(wind_rad, 110),
                    arrowprops=dict(arrowstyle='->', color='red', lw=2))
        ax.text(wind_rad, 112, 'Wind', color='red', fontsize=10,
                fontweight='bold', ha='center')

        ax.set_title('Tower Exposure by Viewing Direction (%)\n'
                     'Higher = more visible, Lower = more sheltered',
                     fontsize=12, fontweight='bold', color='white', pad=20)

        polar_path = os.path.join(script_dir, 'exposure_polar.png')
        fig.savefig(polar_path, dpi=150, bbox_inches='tight', facecolor='#0d1117')
        plt.close(fig)
        print(f"\n  Saved: {polar_path}")

    print("\nDone.")


if __name__ == '__main__':
    main()
