#!/usr/bin/env python3
"""
render_compass_views.py
=======================
Pre-render elevation views of the Citicorp city model from 24 compass bearings
(every 15 degrees). Outputs PNG frames + an interactive HTML viewer.

Usage:
    python render_compass_views.py [--stl-dir PATH] [--resolution 4.0]
"""

import os, sys, math, struct, argparse, base64, io
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.patches import Circle
except ImportError:
    print("ERROR: pip install matplotlib numpy")
    sys.exit(1)

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

# ─── COMPASS DIRECTIONS ──────────────────────────────────────────────────────

COMPASS_LABELS = {
    0: 'N', 15: 'NNE', 30: 'NNE', 45: 'NE', 60: 'ENE', 75: 'ENE',
    90: 'E', 105: 'ESE', 120: 'ESE', 135: 'SE', 150: 'SSE', 165: 'SSE',
    180: 'S', 195: 'SSW', 210: 'SSW', 225: 'SW', 240: 'WSW', 255: 'WSW',
    270: 'W', 285: 'WNW', 300: 'WNW', 315: 'NW', 330: 'NNW', 345: 'NNW',
}

def bearing_label(deg):
    """Short compass label for a bearing."""
    return COMPASS_LABELS.get(deg % 360, f'{deg}')

# ─── RENDER ONE FRAME ────────────────────────────────────────────────────────

def render_frame(bearing_deg, tagged_tris, res, extent_lat, extent_z,
                 n_surr, n_total):
    """
    Render elevation view from a given compass bearing.

    bearing_deg : compass bearing of the viewer (0=N, 90=E, 180=S, 270=W)
    tagged_tris : list of (tri, x_cent, y_cent, is_citicorp)
    res         : grid resolution in meters
    extent_lat  : (lat_min, lat_max) in meters along lateral axis
    extent_z    : (z_min, z_max) in meters

    Returns: matplotlib figure
    """
    theta = math.radians(bearing_deg)

    # Viewing direction (from viewer toward origin):
    #   u = (-sin(theta), -cos(theta))
    # Lateral direction (viewer's right):
    #   r = (-cos(theta), sin(theta))
    sin_t = math.sin(theta)
    cos_t = math.cos(theta)

    lat_min, lat_max = extent_lat
    z_min, z_max = extent_z

    n_lat = int((lat_max - lat_min) / res) + 1
    n_z   = int((z_max - z_min) / res) + 1

    # Compute projection for all triangles
    projected = []
    for tri, xc, yc, is_citi in tagged_tris:
        # Depth along viewing direction (larger = further from viewer)
        depth = xc * sin_t + yc * cos_t

        # Project each vertex to (lateral, z)
        verts_2d = []
        for vx, vy, vz in tri:
            lat = vx * (-cos_t) + vy * sin_t
            verts_2d.append((lat, vz))

        projected.append((depth, verts_2d, is_citi))

    # Sort back-to-front: largest depth first (furthest from viewer drawn first)
    projected.sort(key=lambda p: p[0], reverse=True)

    # Rasterize
    cat_grid   = np.full((n_z, n_lat), np.nan, dtype=np.float32)
    depth_grid = np.full((n_z, n_lat), np.nan, dtype=np.float32)

    depth_vals = [p[0] for p in projected]
    d_min = min(depth_vals) if depth_vals else 0
    d_max = max(depth_vals) if depth_vals else 1
    d_range = d_max - d_min + 1e-6

    for depth, verts_2d, is_citi in projected:
        l0, z0 = verts_2d[0]
        l1, z1 = verts_2d[1]
        l2, z2 = verts_2d[2]

        z_max_t = max(z0, z1, z2)
        if z_max_t < 0.5:
            continue

        l_min_t = min(l0, l1, l2)
        l_max_t = max(l0, l1, l2)
        z_min_t = min(z0, z1, z2)

        li0 = max(0, int((l_min_t - lat_min) / res))
        li1 = min(n_lat - 1, int((l_max_t - lat_min) / res) + 1)
        zi0 = max(0, int((z_min_t - z_min) / res))
        zi1 = min(n_z - 1, int((z_max_t - z_min) / res) + 1)

        # Normalized depth: 0 = far, 1 = close to viewer
        d_norm = 1.0 - (depth - d_min) / d_range

        for zi in range(zi0, zi1 + 1):
            for li in range(li0, li1 + 1):
                pl = lat_min + li * res
                pz = z_min + zi * res
                if _pt_in_tri(pl, pz, l0, z0, l1, z1, l2, z2):
                    cat_grid[zi, li] = 1.0 if is_citi else 0.0
                    depth_grid[zi, li] = d_norm

    # ── Build figure ─────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(14, 5.5))
    fig.patch.set_facecolor('#0d1117')
    ax.set_facecolor('#0f1318')

    extent = [lat_min, lat_max, z_min, z_max]

    # RGBA image
    rgba = np.zeros((n_z, n_lat, 4), dtype=np.float32)

    surr_mask = cat_grid == 0.0
    cmap_surr = plt.cm.copper
    if np.any(surr_mask):
        d_s = np.where(surr_mask, depth_grid, 0.0)
        rgba[surr_mask] = cmap_surr(d_s)[surr_mask]

    citi_mask = cat_grid == 1.0
    if np.any(citi_mask):
        rgba[citi_mask] = [0.0, 0.90, 1.0, 1.0]

    ax.imshow(rgba, extent=extent, origin='lower', aspect='auto',
              interpolation='nearest', zorder=1)

    ax.axhline(y=0, color='#555555', linewidth=1, zorder=3)

    # Height references
    for h in [100, 200]:
        if h < z_max:
            ax.axhline(y=h, color='#333333', linewidth=0.5, linestyle=':', zorder=0)
            ax.text(lat_max - 10, h + 3, f'{h}m', color='#555', fontsize=7,
                    ha='right', zorder=0)

    # Citicorp label
    ax.annotate('Citicorp\n278.9m', xy=(0, 279), xytext=(0, 300),
                color='#00e5ff', fontsize=10, fontweight='bold',
                ha='center', va='bottom',
                arrowprops=dict(arrowstyle='->', color='#00e5ff', lw=1.5),
                zorder=5)

    # Wind into-page symbols
    for wy in np.linspace(lat_min + 80, lat_max - 80, 4):
        wz = z_max - 18
        r = 6.0
        circ = Circle((wy, wz), r, fill=False, edgecolor='#00e5ff',
                       linewidth=1.2, zorder=7)
        ax.add_patch(circ)
        d = r * 0.55
        ax.plot([wy-d, wy+d], [wz-d, wz+d], color='#00e5ff', lw=1.2, zorder=7)
        ax.plot([wy-d, wy+d], [wz+d, wz-d], color='#00e5ff', lw=1.2, zorder=7)

    # Determine left/right labels based on bearing
    # Viewer's left = one compass direction, right = another
    left_bearing = (bearing_deg + 90) % 360
    right_bearing = (bearing_deg - 90) % 360
    left_label = bearing_label(int(left_bearing))
    right_label = bearing_label(int(right_bearing))

    wind_label = f"Wind from {bearing_label(bearing_deg)} ({bearing_deg}) into page"
    ax.text(0, z_max - 6, wind_label,
            color='#00e5ff', fontsize=10, fontweight='bold',
            ha='center', va='top', zorder=7,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#0d1117',
                      edgecolor='#00e5ff', alpha=0.85))

    ax.set_xlabel(f'{left_label}  <--  Position (m)  -->  {right_label}',
                  fontsize=11, color='white')
    ax.set_ylabel('Height (m)', fontsize=11, color='white')
    ax.tick_params(colors='#cccccc', labelsize=9)

    # Dataset info
    if n_surr > 50000:
        ds = "1978 Historical"
    elif n_surr > 25000:
        ds = "2024 Present-Day"
    else:
        ds = "City Model"

    ax.set_title(
        f'Citicorp — View from {bearing_label(bearing_deg)} ({bearing_deg})\n'
        f'{ds}  |  {n_total:,} triangles',
        fontsize=12, fontweight='bold', color='white', pad=10)

    for spine in ax.spines.values():
        spine.set_color('#444444')

    # Invert so viewer's left is on the left of the image
    ax.set_xlim(lat_max, lat_min)
    ax.set_ylim(z_min, min(z_max, 330))

    plt.tight_layout()
    return fig

# ─── MAIN ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Render compass elevation views')
    parser.add_argument('--stl-dir', default=None)
    parser.add_argument('--resolution', type=float, default=4.0,
                        help='Grid cell size in meters (default: 4.0)')
    parser.add_argument('--step', type=int, default=15,
                        help='Degrees between frames (default: 15)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    stl_dir = args.stl_dir or os.path.join(script_dir, 'constant', 'triSurface')

    print("\nCompass Elevation Viewer Generator")
    print("=" * 55)

    # Load STLs once
    tower_tris = read_stl(os.path.join(stl_dir, 'citicorp_tower.stl'))
    stilt_tris = read_stl(os.path.join(stl_dir, 'citicorp_stilts.stl'))
    surr_tris  = read_stl(os.path.join(stl_dir, 'surroundings.stl'))

    n_surr = len(surr_tris)
    all_tris = tower_tris + stilt_tris + surr_tris
    n_total = len(all_tris)
    print(f"  Loaded {n_total:,} triangles")

    if not all_tris:
        print("No triangles!"); sys.exit(1)

    # Pre-compute centroids
    tagged = []
    for tri in surr_tris:
        xc = (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0
        yc = (tri[0][1] + tri[1][1] + tri[2][1]) / 3.0
        tagged.append((tri, xc, yc, False))
    for tri in tower_tris + stilt_tris:
        xc = (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0
        yc = (tri[0][1] + tri[1][1] + tri[2][1]) / 3.0
        tagged.append((tri, xc, yc, True))

    # Find lateral extent (max possible across all bearings)
    xs = [v[0] for t in all_tris for v in t]
    ys = [v[1] for t in all_tris for v in t]
    zs = [v[2] for t in all_tris for v in t]
    max_r = max(max(abs(min(xs)), abs(max(xs))),
                max(abs(min(ys)), abs(max(ys))))
    pad = 30
    lat_ext = (-max_r - pad, max_r + pad)
    z_ext = (0.0, max(zs) + 20.0)

    print(f"  Lateral extent: [{lat_ext[0]:.0f}, {lat_ext[1]:.0f}] m")
    print(f"  Height extent:  [{z_ext[0]:.0f}, {z_ext[1]:.0f}] m")
    print(f"  Resolution:     {args.resolution} m")

    # Render frames
    bearings = list(range(0, 360, args.step))
    n_frames = len(bearings)
    print(f"\nRendering {n_frames} frames (every {args.step} degrees)...")

    frames_b64 = {}
    for i, bearing in enumerate(bearings):
        print(f"  [{i+1:2d}/{n_frames}] {bearing:3d} ({bearing_label(bearing)})...",
              end='', flush=True)
        fig = render_frame(bearing, tagged, args.resolution,
                           lat_ext, z_ext, n_surr, n_total)
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=100, bbox_inches='tight',
                    facecolor='#0d1117')
        plt.close(fig)
        b64 = base64.b64encode(buf.getvalue()).decode('ascii')
        frames_b64[bearing] = b64
        sz_kb = len(b64) * 3 / 4 / 1024
        print(f" {sz_kb:.0f} KB")

    total_mb = sum(len(v) for v in frames_b64.values()) * 3 / 4 / 1024 / 1024
    print(f"\n  Total image data: {total_mb:.1f} MB")

    # ── Generate HTML ────────────────────────────────────────────────────
    print("\nBuilding HTML viewer...")

    # JSON-encode the frame data
    frames_js = "{\n"
    for b in bearings:
        frames_js += f'  {b}: "data:image/png;base64,{frames_b64[b]}",\n'
    frames_js += "}"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Citicorp CFD -- Compass Elevation Viewer</title>
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ background: #0d1117; color: #e6e6e6; font-family: 'Segoe UI', sans-serif;
         display: flex; flex-direction: column; align-items: center; padding: 20px; }}
  h1 {{ font-size: 1.3em; margin-bottom: 10px; color: #00e5ff; }}
  .controls {{ display: flex; align-items: center; gap: 30px; margin: 15px 0; }}
  .compass {{ position: relative; width: 180px; height: 180px; }}
  .compass canvas {{ width: 180px; height: 180px; }}
  .slider-group {{ display: flex; flex-direction: column; align-items: center; gap: 6px; }}
  .slider-group input {{ width: 300px; accent-color: #00e5ff; }}
  .bearing-label {{ font-size: 1.8em; font-weight: bold; color: #00e5ff; min-width: 120px; text-align: center; }}
  .frame {{ border: 1px solid #333; margin-top: 10px; max-width: 95vw; }}
  .frame img {{ display: block; max-width: 95vw; height: auto; }}
  .info {{ font-size: 0.85em; color: #888; margin-top: 8px; }}
  .play-btn {{ padding: 6px 18px; background: #00e5ff; color: #0d1117; border: none;
               border-radius: 4px; font-size: 1em; cursor: pointer; font-weight: bold; }}
  .play-btn:hover {{ background: #33eeff; }}
</style>
</head>
<body>

<h1>Citicorp Center -- Compass Elevation Viewer</h1>
<p class="info">Drag the slider or click the compass to change viewing direction.
   Buildings between you and the tower occlude its lower portion.</p>

<div class="controls">
  <div class="compass">
    <canvas id="compassCanvas" width="180" height="180"></canvas>
  </div>
  <div class="slider-group">
    <div class="bearing-label" id="bearingLabel">270 (W)</div>
    <input type="range" id="bearingSlider" min="0" max="345" step="{args.step}" value="270">
    <div style="display:flex; gap:10px; margin-top:4px;">
      <button class="play-btn" id="playBtn" onclick="togglePlay()">Play</button>
      <button class="play-btn" onclick="setBearing(270)">Wind Dir (W)</button>
    </div>
  </div>
</div>

<div class="frame">
  <img id="frameImg" src="" alt="elevation view">
</div>

<script>
const FRAMES = {frames_js};
const STEP = {args.step};
const BEARINGS = Object.keys(FRAMES).map(Number).sort((a,b)=>a-b);
const LABELS = {{0:'N',45:'NE',90:'E',135:'SE',180:'S',225:'SW',270:'W',315:'NW'}};

let currentBearing = 270;
let playing = false;
let playTimer = null;

const slider = document.getElementById('bearingSlider');
const label = document.getElementById('bearingLabel');
const img = document.getElementById('frameImg');
const canvas = document.getElementById('compassCanvas');
const ctx = canvas.getContext('2d');

function bearingName(deg) {{
  deg = ((deg % 360) + 360) % 360;
  const names = ['N','NNE','NE','ENE','E','ESE','SE','SSE',
                 'S','SSW','SW','WSW','W','WNW','NW','NNW'];
  return names[Math.round(deg / 22.5) % 16];
}}

function setBearing(deg) {{
  deg = ((deg % 360) + 360) % 360;
  // Snap to nearest step
  deg = Math.round(deg / STEP) * STEP % 360;
  currentBearing = deg;
  slider.value = deg;
  label.textContent = deg + ' (' + bearingName(deg) + ')';
  img.src = FRAMES[deg];
  drawCompass(deg);
}}

function drawCompass(bearing) {{
  const cx = 90, cy = 90, r = 75;
  ctx.clearRect(0, 0, 180, 180);

  // Outer ring
  ctx.beginPath();
  ctx.arc(cx, cy, r, 0, Math.PI * 2);
  ctx.strokeStyle = '#444';
  ctx.lineWidth = 2;
  ctx.stroke();

  // Tick marks and labels
  for (let deg = 0; deg < 360; deg += 15) {{
    const rad = (deg - 90) * Math.PI / 180;
    const inner = deg % 45 === 0 ? r - 12 : r - 6;
    ctx.beginPath();
    ctx.moveTo(cx + Math.cos(rad) * inner, cy + Math.sin(rad) * inner);
    ctx.lineTo(cx + Math.cos(rad) * r, cy + Math.sin(rad) * r);
    ctx.strokeStyle = deg % 90 === 0 ? '#aaa' : '#555';
    ctx.lineWidth = deg % 45 === 0 ? 2 : 1;
    ctx.stroke();
  }}

  // Cardinal labels
  ctx.font = 'bold 13px sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  const cardinals = [['N', 0],['E', 90],['S', 180],['W', 270]];
  for (const [lbl, deg] of cardinals) {{
    const rad = (deg - 90) * Math.PI / 180;
    const lr = r + 12;
    ctx.fillStyle = deg === 270 ? '#00e5ff' : '#aaa';
    ctx.fillText(lbl, cx + Math.cos(rad) * lr, cy + Math.sin(rad) * lr);
  }}

  // Viewer direction arrow
  const bRad = (bearing - 90) * Math.PI / 180;
  ctx.beginPath();
  ctx.moveTo(cx, cy);
  ctx.lineTo(cx + Math.cos(bRad) * (r - 15), cy + Math.sin(bRad) * (r - 15));
  ctx.strokeStyle = '#00e5ff';
  ctx.lineWidth = 3;
  ctx.stroke();

  // Arrowhead
  const ax = cx + Math.cos(bRad) * (r - 15);
  const ay = cy + Math.sin(bRad) * (r - 15);
  const aSize = 8;
  ctx.beginPath();
  ctx.moveTo(ax, ay);
  ctx.lineTo(ax - Math.cos(bRad - 0.4) * aSize, ay - Math.sin(bRad - 0.4) * aSize);
  ctx.moveTo(ax, ay);
  ctx.lineTo(ax - Math.cos(bRad + 0.4) * aSize, ay - Math.sin(bRad + 0.4) * aSize);
  ctx.stroke();

  // Center dot (Citicorp)
  ctx.beginPath();
  ctx.arc(cx, cy, 4, 0, Math.PI * 2);
  ctx.fillStyle = '#00e5ff';
  ctx.fill();

  // "Wind from" label
  ctx.font = '9px sans-serif';
  ctx.fillStyle = '#888';
  ctx.fillText('Viewer', cx, cy + r + 12);
}}

// Slider event
slider.addEventListener('input', () => setBearing(parseInt(slider.value)));

// Compass click
canvas.addEventListener('click', function(e) {{
  const rect = canvas.getBoundingClientRect();
  const x = e.clientX - rect.left - 90;
  const y = e.clientY - rect.top - 90;
  let deg = Math.atan2(x, -y) * 180 / Math.PI;
  deg = ((deg % 360) + 360) % 360;
  setBearing(deg);
}});

// Play/pause animation
function togglePlay() {{
  playing = !playing;
  document.getElementById('playBtn').textContent = playing ? 'Pause' : 'Play';
  if (playing) {{
    playTimer = setInterval(() => {{
      setBearing((currentBearing + STEP) % 360);
    }}, 400);
  }} else {{
    clearInterval(playTimer);
  }}
}}

// Keyboard arrows
document.addEventListener('keydown', (e) => {{
  if (e.key === 'ArrowRight') setBearing(currentBearing + STEP);
  else if (e.key === 'ArrowLeft') setBearing(currentBearing - STEP);
  else if (e.key === ' ') {{ e.preventDefault(); togglePlay(); }}
}});

// Initialize
setBearing(270);
</script>
</body>
</html>"""

    html_path = os.path.join(script_dir, 'compass_viewer.html')
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html)
    html_mb = os.path.getsize(html_path) / 1024 / 1024
    print(f"  Saved: {html_path} ({html_mb:.1f} MB)")
    print(f"\nOpen compass_viewer.html in a browser.")
    print("Controls: slider, compass click, arrow keys, spacebar for play/pause.")


if __name__ == '__main__':
    main()
