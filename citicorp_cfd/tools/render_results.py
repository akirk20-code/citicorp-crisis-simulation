#!/usr/bin/env pvpython
"""Render CFD results as PNG screenshots using ParaView's offscreen rendering."""
import os
import sys

from paraview.simple import *

# Setup offscreen rendering
paraview.simple._DisableFirstRenderCameraReset()

case_dir = '/home/kirka/citicorp_cfd'
output_dir = os.path.join(case_dir, 'figures')
os.makedirs(output_dir, exist_ok=True)

# Open the OpenFOAM case
foam_file = os.path.join(case_dir, 'citicorp_cfd.foam')
reader = OpenFOAMReader(FileName=foam_file)
reader.MeshRegions = ['internalMesh']
reader.CellArrays = ['U', 'p', 'k']

# Get available time steps and use the latest
timesteps = reader.TimestepValues
print(f"Available timesteps: {timesteps}")
latest_time = max(timesteps) if timesteps else 0
print(f"Using timestep: {latest_time}")

# Create render view
view = CreateRenderView()
view.ViewSize = [1920, 1080]
view.Background = [0.05, 0.05, 0.1]  # Dark background

# Camera positions (building centered at origin, Z-up in OpenFOAM)
cameras = {
    'perspective': {
        'position': [-500, -400, 350],
        'focal': [0, 0, 120],
        'up': [0, 0, 1],
    },
    'top': {
        'position': [0, 0, 800],
        'focal': [0, 0, 0],
        'up': [0, 1, 0],
    },
    'front': {
        'position': [-600, 0, 150],
        'focal': [0, 0, 150],
        'up': [0, 0, 1],
    },
    'side': {
        'position': [0, -600, 150],
        'focal': [0, 0, 150],
        'up': [0, 0, 1],
    },
}

# ===== Velocity magnitude =====
print("Rendering velocity magnitude...")
display = Show(reader, view)
ColorBy(display, ('CELLS', 'U', 'Magnitude'))
display.SetScalarBarVisibility(view, True)

# Get color transfer function and set range
uLUT = GetColorTransferFunction('U')
uLUT.RescaleTransferFunction(0, 80)  # 0 to 80 m/s
uLUT.ApplyPreset('Cool to Warm (Extended)', True)

# Get scalar bar and style it
uBar = GetScalarBar(uLUT, view)
uBar.Title = 'Velocity Magnitude (m/s)'
uBar.ComponentTitle = ''
uBar.TitleFontSize = 16
uBar.LabelFontSize = 14
uBar.ScalarBarLength = 0.4
uBar.Position = [0.85, 0.25]

# Add time annotation
annotate = AnnotateTimeFilter(reader)
annotate_display = Show(annotate, view)
annotate_display.FontSize = 14

# Set to latest timestep
view.ViewTime = latest_time
reader.UpdatePipeline(latest_time)

for cam_name, cam in cameras.items():
    view.CameraPosition = cam['position']
    view.CameraFocalPoint = cam['focal']
    view.CameraViewUp = cam['up']
    view.ResetCamera()
    # Zoom in a bit
    view.CameraPosition = cam['position']
    view.CameraFocalPoint = cam['focal']
    Render()

    fname = os.path.join(output_dir, f'velocity_{cam_name}.png')
    SaveScreenshot(fname, view, ImageResolution=[1920, 1080])
    print(f"  Saved: {fname}")

# ===== Pressure =====
print("Rendering pressure...")
ColorBy(display, ('CELLS', 'p'))
display.SetScalarBarVisibility(view, True)

pLUT = GetColorTransferFunction('p')
pLUT.RescaleTransferFunction(-2000, 500)  # Pressure range (Pa, kinematic)
pLUT.ApplyPreset('Blue to Red Rainbow', True)

pBar = GetScalarBar(pLUT, view)
pBar.Title = 'Pressure (m²/s²)'
pBar.ComponentTitle = ''
pBar.TitleFontSize = 16
pBar.LabelFontSize = 14
pBar.ScalarBarLength = 0.4
pBar.Position = [0.85, 0.25]

for cam_name, cam in cameras.items():
    view.CameraPosition = cam['position']
    view.CameraFocalPoint = cam['focal']
    view.CameraViewUp = cam['up']
    Render()

    fname = os.path.join(output_dir, f'pressure_{cam_name}.png')
    SaveScreenshot(fname, view, ImageResolution=[1920, 1080])
    print(f"  Saved: {fname}")

# ===== Slice at stilt level (z=30m) for velocity =====
print("Rendering horizontal slice at stilt level (z=30m)...")
Hide(reader, view)

sliceFilter = Slice(Input=reader)
sliceFilter.SliceType = 'Plane'
sliceFilter.SliceType.Origin = [0, 0, 30]
sliceFilter.SliceType.Normal = [0, 0, 1]

sliceDisplay = Show(sliceFilter, view)
ColorBy(sliceDisplay, ('CELLS', 'U', 'Magnitude'))
sliceDisplay.SetScalarBarVisibility(view, True)

uLUT.RescaleTransferFunction(0, 60)

# Top view for slice
view.CameraPosition = [0, 0, 800]
view.CameraFocalPoint = [0, 0, 30]
view.CameraViewUp = [0, 1, 0]
Render()

fname = os.path.join(output_dir, 'velocity_slice_stilt_level.png')
SaveScreenshot(fname, view, ImageResolution=[1920, 1080])
print(f"  Saved: {fname}")

# ===== Slice at mid-height (z=124m) =====
print("Rendering horizontal slice at mid-height (z=124m)...")
sliceFilter.SliceType.Origin = [0, 0, 124]
sliceFilter.UpdatePipeline(latest_time)
Render()

fname = os.path.join(output_dir, 'velocity_slice_midheight.png')
SaveScreenshot(fname, view, ImageResolution=[1920, 1080])
print(f"  Saved: {fname}")

# ===== Vertical slice through building center (y=0) =====
print("Rendering vertical slice through building center...")
sliceFilter.SliceType.Origin = [0, 0, 0]
sliceFilter.SliceType.Normal = [0, 1, 0]
sliceFilter.UpdatePipeline(latest_time)

# Side view for vertical slice
view.CameraPosition = [0, -600, 150]
view.CameraFocalPoint = [0, 0, 150]
view.CameraViewUp = [0, 0, 1]
Render()

fname = os.path.join(output_dir, 'velocity_slice_vertical.png')
SaveScreenshot(fname, view, ImageResolution=[1920, 1080])
print(f"  Saved: {fname}")

# Front view of vertical slice
view.CameraPosition = [-600, 0, 150]
view.CameraFocalPoint = [0, 0, 150]
view.CameraViewUp = [0, 0, 1]
Render()

fname = os.path.join(output_dir, 'velocity_slice_vertical_front.png')
SaveScreenshot(fname, view, ImageResolution=[1920, 1080])
print(f"  Saved: {fname}")

# ===== Streamlines =====
print("Rendering streamlines...")
Hide(sliceFilter, view)
Show(reader, view)
ColorBy(GetDisplayProperties(reader, view), ('CELLS', 'U', 'Magnitude'))

streamTracer = StreamTracerWithCustomSource(Input=reader,
    SeedSource='Point Cloud')
streamTracer.Vectors = ['CELLS', 'U']
streamTracer.MaximumStreamlineLength = 1000
streamTracer.SeedSource.NumberOfPoints = 200
streamTracer.SeedSource.Center = [-300, 0, 100]
streamTracer.SeedSource.Radius = 150

streamDisplay = Show(streamTracer, view)
ColorBy(streamDisplay, ('POINTS', 'U', 'Magnitude'))
streamDisplay.LineWidth = 1.5

# Perspective view
view.CameraPosition = [-500, -400, 350]
view.CameraFocalPoint = [0, 0, 120]
view.CameraViewUp = [0, 0, 1]
Render()

fname = os.path.join(output_dir, 'streamlines_perspective.png')
SaveScreenshot(fname, view, ImageResolution=[1920, 1080])
print(f"  Saved: {fname}")

print(f"\nAll renders complete! Files saved to: {output_dir}")
print("Copy to Windows: cp -r figures/ /mnt/c/.../citicorp_cfd/figures/")
