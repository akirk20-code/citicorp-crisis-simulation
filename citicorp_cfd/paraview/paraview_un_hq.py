#!/usr/bin/env python3
"""ParaView visualization of UN HQ surroundings mesh.

Run with:
    "C:/Program Files/ParaView 6.1.0/bin/pvpython.exe" _paraview_un_hq.py

Produces PNG screenshots in the parent directory (Citi Sim/).
"""

from paraview.simple import *
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TRI_DIR = os.path.join(SCRIPT_DIR, 'constant', 'triSurface')
OUT_DIR = os.path.dirname(SCRIPT_DIR)  # Citi Sim/

# Color scheme
COLOR_SURR = (0.55, 0.75, 0.55)      # green tint for CityGML LoD2
BG_COLOR = (1.0, 1.0, 1.0)

# Camera presets for UN HQ area (buildings ~170m max, 600m domain)
CAMERAS = {
    'aerial': {
        'pos': (100, -250, 400),
        'focus': (0, 0, 40),
        'up': (0, 0, 1),
        'label': 'Aerial View',
    },
    'quarter': {
        'pos': (200, -200, 180),
        'focus': (0, 0, 50),
        'up': (0, 0, 1),
        'label': '3/4 View',
    },
    'street': {
        'pos': (-60, -50, 15),
        'focus': (0, 0, 50),
        'up': (0, 0, 1),
        'label': 'Street Level',
    },
    'top_down': {
        'pos': (0, 0, 500),
        'focus': (0, 0, 0),
        'up': (0, 1, 0),
        'label': 'Plan View (Top Down)',
    },
    'east_river': {
        'pos': (250, 50, 60),
        'focus': (-20, 0, 50),
        'up': (0, 0, 1),
        'label': 'East River Approach',
    },
}

IMG_W = 1920
IMG_H = 1080


def load_stl(filepath, color, opacity=1.0):
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
    view.CameraPosition = cam_dict['pos']
    view.CameraFocalPoint = cam_dict['focus']
    view.CameraViewUp = cam_dict['up']
    view.CameraParallelProjection = 0
    Render()


def add_text(view, text, position=(0.02, 0.92), fontsize=18, color=(0, 0, 0)):
    txt = Text(Text=text)
    txt_display = Show(txt, view)
    txt_display.WindowLocation = 'Any Location'
    txt_display.Position = list(position)
    txt_display.FontSize = fontsize
    txt_display.Color = color
    txt_display.Bold = 1
    return txt


def screenshot(view, filename, width=IMG_W, height=IMG_H):
    filepath = os.path.join(OUT_DIR, filename)
    layout = GetLayout()
    SaveScreenshot(filepath, layout,
                   ImageResolution=[width, height],
                   TransparentBackground=0)
    print(f"  Saved: {filepath}")


def main():
    print("=" * 60)
    print("ParaView UN HQ Visualization")
    print("=" * 60)

    stl_path = os.path.join(TRI_DIR, 'un_hq_surroundings.stl')

    for cam_name, cam in CAMERAS.items():
        print(f"\n--- {cam['label']} ---")
        ResetSession()
        view = GetActiveViewOrCreate('RenderView')
        view.Background = list(BG_COLOR)
        view.OrientationAxesVisibility = 0

        load_stl(stl_path, COLOR_SURR, opacity=0.9)
        add_text(view, f"UN HQ -- CityGML Hybrid ({cam['label']})")
        set_camera(view, cam)
        screenshot(view, f"un_hq_{cam_name}.png")

    print(f"\n{'=' * 60}")
    print(f"Done! Screenshots saved to: {OUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
