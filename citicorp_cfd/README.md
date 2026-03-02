# Citicorp Center CFD Simulation

Computational Fluid Dynamics simulation of atmospheric boundary layer flow around the Citicorp Center in Midtown Manhattan. Uses OpenFOAM for the simulation and Python for procedural geometry generation from NYC Open Data.

## Quick Start

```bash
# Install Python dependencies
pip install -r requirements.txt

# Run full simulation (geometry + mesh + solve)
./Allrun

# Clean all generated files
./Allclean
```

**Important:** OpenFOAM cannot handle spaces in paths. Copy the case to a space-free location before running:
```bash
cp -r "<windows path>" ~/citicorp_cfd
cd ~/citicorp_cfd && bash Allrun
```

## Prerequisites

- **OpenFOAM:** Foundation v13 (`foamRun -solver incompressibleFluid`) or ESI v2512 (`simpleFoam`)
- **Python:** 3.8+ with `requests` module
- **System:** WSL2 Ubuntu (tested on 24.04)

## Directory Structure

```
Allrun, Allclean          Run/clean scripts
0/                        Initial and boundary conditions (ABL inlet)
constant/                 Physical properties and geometry
  triSurface/             STL files for snappyHexMesh
system/                   Solver control dictionaries

generators/               STL geometry generation scripts
  generate_stl.py         Original (hardcoded fallback)
  generate_stl_nyc3d.py   NYC Open Data API (RECOMMENDED)
  generate_stl_hybrid.py  CityGML + API (most flexible, supports --year)
  README.md               Comparison and selection guide

tools/                    Analysis and utility scripts
paraview/                 ParaView visualization + presentation generators
un_hq/                    UN Headquarters extraction demo
dashboards/               Interactive HTML viewers
docs/                     Technical documentation
```

## Simulation Details

- **Solver:** `foamRun -solver incompressibleFluid` (steady-state RANS)
- **Turbulence:** k-omega SST
- **Inlet:** Atmospheric boundary layer (`atmBoundaryLayerInletVelocity`, z0=0.1m, Uref=44.7 m/s)
- **Mesh:** 3.16M cells (snappyHexMesh, 3 refinement levels)
- **Hardware:** 8 cores meshing, 10 cores solving (~4 hours on Intel Ultra 9 185H)
- **Validation:** Cd = 3.24, within 8% of empirical correlations

## Documentation

See [`docs/`](docs/) for detailed technical notes:
- [ASSUMPTIONS.md](docs/ASSUMPTIONS.md) — Simulation assumptions and limitations
- [CD_VALIDATION.md](docs/CD_VALIDATION.md) — Drag coefficient hand calculation
- [VELOCITY_ANALYSIS.md](docs/VELOCITY_ANALYSIS.md) — ABL inlet profile analysis
- [IMPROVEMENT_ROADMAP.md](docs/IMPROVEMENT_ROADMAP.md) — Future enhancement plan
- [BUILDING_HEIGHTS.md](docs/BUILDING_HEIGHTS.md) — Surrounding building height analysis
- [MESH_SIZING_GUIDE.md](docs/MESH_SIZING_GUIDE.md) — Mesh refinement guidance
- [HOPPER_SETUP.md](docs/HOPPER_SETUP.md) — GMU HPC cluster setup
