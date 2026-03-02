#!/usr/bin/env python3
"""ParaView comparison: old (Socrata LoD1) vs new (CityGML hybrid) cityscape.

Run with:
    "C:/Program Files/ParaView 6.1.0/bin/pvpython.exe" _paraview_comparison.py

Produces PNG screenshots in the parent directory (Citi Sim/).
"""

from paraview.simple import *
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TRI_DIR = os.path.join(SCRIPT_DIR, 'constant', 'triSurface')
OUT_DIR = os.path.dirname(SCRIPT_DIR)  # Citi Sim/

# Color scheme (R, G, B) in 0-1 range
COLOR_TOWER = (0.85, 0.20, 0.20)    # red
COLOR_STILTS = (0.20, 0.40, 0.85)   # blue
COLOR_SURR = (0.70, 0.70, 0.70)     # light gray
COLOR_LOD2_SURR = (0.55, 0.75, 0.55) # green tint for CityGML
BG_COLOR = (1.0, 1.0, 1.0)          # white background

# Camera presets: (position, focal_point, view_up)
CAMERAS = {
    'aerial': {
        'pos': (150, -300, 500),
        'focus': (50, 0, 50),
        'up': (0, 0, 1),
        'label': 'Aerial View',
    },
    'street': {
        'pos': (-80, -60, 20),
        'focus': (0, 0, 40),
        'up': (0, 0, 1),
        'label': 'Street Level (Stilt Zone)',
    },
    'quarter': {
        'pos': (250, -250, 200),
        'focus': (0, 0, 80),
        'up': (0, 0, 1),
        'label': '3/4 View',
    },
    'top_down': {
        'pos': (0, 0, 600),
        'focus': (0, 0, 0),
        'up': (0, 1, 0),
        'label': 'Plan View (Top Down)',
    },
}

IMG_W = 1920
IMG_H = 1080


def load_stl(filepath, color, opacity=1.0):
    """Load an STL file and set display properties."""
    if not os.path.exists(filepath):
        print(f"  WARNING: {filepath} not found, skipping")
        return None
    reader = STLReader(FileNames=[filepath])
    display = Show(reader)
    display.Representation = 'Surface With Edges'
    display.AmbientColor = color
    display.DiffuseColor = color
    display.Opacity = opacity
    display.EdgeColor = (0.2, 0.2, 0.2)
    display.LineWidth = 0.5
    name = os.path.basename(filepath)
    n_cells = reader.GetDataInformation().GetNumberOfCells()
    print(f"  Loaded {name}: {n_cells} triangles")
    return reader


def set_camera(view, cam_dict):
    """Set camera position and orientation."""
    view.CameraPosition = cam_dict['pos']
    view.CameraFocalPoint = cam_dict['focus']
    view.CameraViewUp = cam_dict['up']
    view.CameraParallelProjection = 0
    Render()


def add_text(view, text, position=(0.02, 0.92), fontsize=18, color=(0, 0, 0)):
    """Add text annotation to the view."""
    txt = Text(Text=text)
    txt_display = Show(txt, view)
    txt_display.WindowLocation = 'Any Location'
    txt_display.Position = list(position)
    txt_display.FontSize = fontsize
    txt_display.Color = color
    txt_display.Bold = 1
    return txt


def screenshot(view, filename, width=IMG_W, height=IMG_H):
    """Save screenshot."""
    filepath = os.path.join(OUT_DIR, filename)
    layout = GetLayout()
    SaveScreenshot(filepath, layout,
                   ImageResolution=[width, height],
                   TransparentBackground=0)
    print(f"  Saved: {filepath}")


def create_scene(file_set, label, color_surr):
    """Load a set of STL files into the current view."""
    sources = []
    tower_path = os.path.join(TRI_DIR, file_set['tower'])
    stilts_path = os.path.join(TRI_DIR, file_set['stilts'])
    surr_path = os.path.join(TRI_DIR, file_set['surroundings'])

    s = load_stl(surr_path, color_surr, opacity=0.85)
    if s:
        sources.append(s)
    s = load_stl(tower_path, COLOR_TOWER)
    if s:
        sources.append(s)
    s = load_stl(stilts_path, COLOR_STILTS)
    if s:
        sources.append(s)

    return sources


def main():
    print("=" * 60)
    print("ParaView Cityscape Comparison")
    print("=" * 60)

    # --- Old geometry (Socrata LoD1) ---
    old_files = {
        'tower': 'citicorp_tower_bin.stl',
        'stilts': 'citicorp_stilts_bin.stl',
        'surroundings': 'surroundings_bin.stl',
    }

    # --- New geometry (CityGML hybrid) ---
    new_files = {
        'tower': 'citicorp_tower.stl',
        'stilts': 'citicorp_stilts.stl',
        'surroundings': 'surroundings.stl',
    }

    # ============================================================
    # INDIVIDUAL VIEWS (full-screen, one geometry set at a time)
    # ============================================================
    for cam_name, cam in CAMERAS.items():
        # --- OLD ---
        print(f"\n--- Old geometry: {cam['label']} ---")
        ResetSession()
        view = GetActiveViewOrCreate('RenderView')
        view.Background = list(BG_COLOR)
        view.OrientationAxesVisibility = 0

        create_scene(old_files, 'Socrata LoD1', COLOR_SURR)
        add_text(view, f"OLD: Socrata LoD1 (flat-top extrusions) - {cam['label']}")
        set_camera(view, cam)
        screenshot(view, f"compare_old_{cam_name}.png")

        # --- NEW ---
        print(f"\n--- New geometry: {cam['label']} ---")
        ResetSession()
        view = GetActiveViewOrCreate('RenderView')
        view.Background = list(BG_COLOR)
        view.OrientationAxesVisibility = 0

        create_scene(new_files, 'CityGML Hybrid', COLOR_LOD2_SURR)
        add_text(view, f"NEW: CityGML Hybrid (actual surfaces) - {cam['label']}")
        set_camera(view, cam)
        screenshot(view, f"compare_new_{cam_name}.png")

    # ============================================================
    # CITICORP DETAIL VIEW (stilts + core close-up)
    # ============================================================
    print(f"\n--- Citicorp stilt detail ---")
    ResetSession()
    view = GetActiveViewOrCreate('RenderView')
    view.Background = list(BG_COLOR)
    view.OrientationAxesVisibility = 0

    # Load just tower + stilts (new geometry)
    tower = load_stl(os.path.join(TRI_DIR, new_files['tower']), COLOR_TOWER, 0.4)
    stilts = load_stl(os.path.join(TRI_DIR, new_files['stilts']), COLOR_STILTS)
    add_text(view, "Citicorp Stilts + Center Core (photo-refined)")

    # Close-up camera for stilt zone
    stilt_cam = {
        'pos': (-60, -50, 25),
        'focus': (0, 0, 17),
        'up': (0, 0, 1),
    }
    set_camera(view, stilt_cam)
    screenshot(view, "compare_stilt_detail.png")

    # ============================================================
    # LOD2 TOWER DETAIL (if available)
    # ============================================================
    lod2_tower = os.path.join(TRI_DIR, 'citicorp_tower_lod2_bin.stl')
    if os.path.exists(lod2_tower):
        print(f"\n--- LoD2 tower comparison ---")
        ResetSession()
        view = GetActiveViewOrCreate('RenderView')
        view.Background = list(BG_COLOR)
        view.OrientationAxesVisibility = 0

        load_stl(os.path.join(TRI_DIR, 'citicorp_tower.stl'),
                 (0.9, 0.3, 0.3), 0.5)  # simplified = translucent red
        load_stl(lod2_tower, (0.3, 0.3, 0.9), 0.7)  # LoD2 = blue
        add_text(view, "Tower: Simplified (red) vs CityGML LoD2 (blue)")

        tower_cam = {
            'pos': (80, -80, 160),
            'focus': (0, 0, 140),
            'up': (0, 0, 1),
        }
        set_camera(view, tower_cam)
        screenshot(view, "compare_tower_lod2.png")

    print(f"\n{'=' * 60}")
    print(f"Done! Screenshots saved to: {OUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
