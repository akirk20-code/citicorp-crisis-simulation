#!/usr/bin/env pvpython
"""Render CFD results v2 — better views with building surfaces + flow slices."""
import os
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

case_dir = '/home/kirka/citicorp_cfd'
output_dir = os.path.join(case_dir, 'figures')
os.makedirs(output_dir, exist_ok=True)

# Open case
reader = OpenFOAMReader(FileName=os.path.join(case_dir, 'citicorp_cfd.foam'))
reader.MeshRegions = ['internalMesh']
reader.CellArrays = ['U', 'p', 'k']
latest_time = max(reader.TimestepValues)
print(f"Timestep: {latest_time}")

view = CreateRenderView()
view.ViewSize = [1920, 1080]
view.Background = [0.08, 0.08, 0.12]
view.ViewTime = latest_time

# ---- Helper: set camera and render ----
def set_camera(view, name):
    cams = {
        'perspective': ([-400, -350, 300], [0, 0, 80], [0, 0, 1]),
        'top':         ([0, 0, 700],       [0, 0, 0],  [0, 1, 0]),
        'front':       ([-500, 0, 120],    [0, 0, 120], [0, 0, 1]),
        'side':        ([0, -500, 120],    [0, 0, 120], [0, 0, 1]),
    }
    pos, foc, up = cams[name]
    view.CameraPosition = pos
    view.CameraFocalPoint = foc
    view.CameraViewUp = up

# ======================================
# 1. HORIZONTAL SLICES — velocity
# ======================================
for z, label in [(10, 'ground_10m'), (30, 'stilt_30m'), (80, 'pedestrian_80m'),
                  (124, 'mid_124m'), (200, 'upper_200m'), (248, 'roof_248m')]:
    print(f"Slice z={z}m ({label})...")
    sl = Slice(Input=reader)
    sl.SliceType = 'Plane'
    sl.SliceType.Origin = [0, 0, z]
    sl.SliceType.Normal = [0, 0, 1]

    d = Show(sl, view)
    ColorBy(d, ('CELLS', 'U', 'Magnitude'))
    d.SetScalarBarVisibility(view, True)
    lut = GetColorTransferFunction('U')
    lut.RescaleTransferFunction(0, 70)
    lut.ApplyPreset('Cool to Warm (Extended)', True)

    bar = GetScalarBar(lut, view)
    bar.Title = f'|U| at z={z}m (m/s)'
    bar.ComponentTitle = ''
    bar.TitleFontSize = 18
    bar.LabelFontSize = 14
    bar.ScalarBarLength = 0.35
    bar.Position = [0.87, 0.3]

    set_camera(view, 'top')
    Render()
    SaveScreenshot(os.path.join(output_dir, f'vel_hz_{label}.png'), view,
                   ImageResolution=[1920, 1080])

    Delete(sl)
    del sl

# ======================================
# 2. VERTICAL SLICES — velocity
# ======================================
for axis, normal, cam in [('xz_y0', [0,1,0], 'front'), ('yz_x0', [1,0,0], 'side')]:
    print(f"Vertical slice {axis}...")
    sl = Slice(Input=reader)
    sl.SliceType = 'Plane'
    sl.SliceType.Origin = [0, 0, 0]
    sl.SliceType.Normal = normal

    d = Show(sl, view)
    ColorBy(d, ('CELLS', 'U', 'Magnitude'))
    d.SetScalarBarVisibility(view, True)
    lut = GetColorTransferFunction('U')
    lut.RescaleTransferFunction(0, 80)
    lut.ApplyPreset('Cool to Warm (Extended)', True)

    bar = GetScalarBar(lut, view)
    bar.Title = '|U| (m/s)'
    bar.ComponentTitle = ''
    bar.TitleFontSize = 18
    bar.LabelFontSize = 14

    set_camera(view, cam)
    Render()
    SaveScreenshot(os.path.join(output_dir, f'vel_vert_{axis}.png'), view,
                   ImageResolution=[1920, 1080])

    Delete(sl)
    del sl

# ======================================
# 3. VERTICAL SLICES — pressure
# ======================================
for axis, normal, cam in [('xz_y0', [0,1,0], 'front'), ('yz_x0', [1,0,0], 'side')]:
    print(f"Pressure vertical slice {axis}...")
    sl = Slice(Input=reader)
    sl.SliceType = 'Plane'
    sl.SliceType.Origin = [0, 0, 0]
    sl.SliceType.Normal = normal

    d = Show(sl, view)
    ColorBy(d, ('CELLS', 'p'))
    d.SetScalarBarVisibility(view, True)
    plut = GetColorTransferFunction('p')
    plut.RescaleTransferFunction(-1500, 500)
    plut.ApplyPreset('Cool to Warm', True)

    bar = GetScalarBar(plut, view)
    bar.Title = 'p (m²/s²)'
    bar.ComponentTitle = ''
    bar.TitleFontSize = 18
    bar.LabelFontSize = 14

    set_camera(view, cam)
    Render()
    SaveScreenshot(os.path.join(output_dir, f'pres_vert_{axis}.png'), view,
                   ImageResolution=[1920, 1080])

    Delete(sl)
    del sl

# ======================================
# 4. COMPOSITE: slice + building surfaces
# ======================================
print("Composite view: horizontal slice + building clip...")

# Horizontal velocity slice at z=30m
sl = Slice(Input=reader)
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0, 0, 30]
sl.SliceType.Normal = [0, 0, 1]
d1 = Show(sl, view)
ColorBy(d1, ('CELLS', 'U', 'Magnitude'))
lut = GetColorTransferFunction('U')
lut.RescaleTransferFunction(0, 70)
lut.ApplyPreset('Cool to Warm (Extended)', True)

# Clip to show only lower portion of mesh (buildings visible)
clip = Clip(Input=reader)
clip.ClipType = 'Plane'
clip.ClipType.Origin = [0, 0, 60]
clip.ClipType.Normal = [0, 0, 1]
d2 = Show(clip, view)
d2.Representation = 'Surface'
d2.Opacity = 0.3
ColorBy(d2, ('CELLS', 'U', 'Magnitude'))

set_camera(view, 'perspective')
d1.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = '|U| at stilt level (m/s)'
bar.TitleFontSize = 18
Render()
SaveScreenshot(os.path.join(output_dir, 'composite_stilt_perspective.png'), view,
               ImageResolution=[1920, 1080])

# Top view of composite
set_camera(view, 'top')
Render()
SaveScreenshot(os.path.join(output_dir, 'composite_stilt_top.png'), view,
               ImageResolution=[1920, 1080])

Delete(clip)
Delete(sl)

# ======================================
# 5. TURBULENT KINETIC ENERGY — horizontal
# ======================================
print("TKE at stilt level...")
sl = Slice(Input=reader)
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0, 0, 30]
sl.SliceType.Normal = [0, 0, 1]
d = Show(sl, view)
ColorBy(d, ('CELLS', 'k'))
d.SetScalarBarVisibility(view, True)
klut = GetColorTransferFunction('k')
klut.RescaleTransferFunction(0, 100)
klut.ApplyPreset('Cool to Warm (Extended)', True)
bar = GetScalarBar(klut, view)
bar.Title = 'TKE k (m²/s²)'
bar.TitleFontSize = 18

set_camera(view, 'top')
Render()
SaveScreenshot(os.path.join(output_dir, 'tke_stilt_top.png'), view,
               ImageResolution=[1920, 1080])

Delete(sl)

print(f"\nAll renders saved to: {output_dir}")
