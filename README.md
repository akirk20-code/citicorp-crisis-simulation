# Citicorp Center Wind Crisis Simulation

**Course:** OR 750 — Reliability, Safety, and Risk (George Mason University)
**Date:** February 2026

A multi-fidelity simulation of the 1978 Citicorp Center structural crisis, combining MATLAB structural analysis with OpenFOAM computational fluid dynamics.

## Overview

In 1978, structural engineer William LeMessurier discovered that the Citicorp Center's chevron bracing system — with bolted connections substituted for the designed welded connections — was vulnerable to collapse under quartering winds with a 1-in-16-year return period. This simulation recreates the engineering analysis that led to a secret emergency repair.

## Analysis Modules

| # | Module | Script | Description |
|---|--------|--------|-------------|
| 1 | 3D Building Geometry | `visualize_building_3d.m` | Structural system visualization |
| 2 | Wind Pressure | `analyze_wind_pressure.m` | ASCE 7-22 wind load analysis |
| 3 | Chevron Bracing FEA | `analyze_chevron_fea.m` | Direct stiffness method FEA |
| 4 | Connection Analysis | `analyze_connections.m` | AISC 360-22 bolt/weld checks |
| 5 | Monte Carlo Reliability | `monte_carlo_reliability.m` | Gumbel wind + structural uncertainty |
| 6 | Validation | `visualize_validation.m` | Results summary and validation |
| 7 | Interactive Structure | `interactive_structure.m` | Interactive 3D structural editor |
| 8-9 | CFD Post-Processing | `process_cfd_results.m` | OpenFOAM integration |
| 10 | 1978 Comparison | `analyze_1978_comparison.m` | Historical vs modern code analysis |
| 11 | Data Provenance | `visualize_provenance.m` | Source verification dashboard |
| 12 | Structural Elevation | `draw_structural_elevation.m` | Engineering drawings |

**Run all modules:** `run_citicorp_simulation.m`

## CFD Simulation

See [`citicorp_cfd/README.md`](citicorp_cfd/README.md) for the OpenFOAM case, including:
- Atmospheric boundary layer simulation of Manhattan
- 3.16M cell mesh with surrounding building context
- Drag coefficient validation (Cd = 3.24, within 8% of empirical)
- [UN Headquarters extraction demo](citicorp_cfd/un_hq/) (CityGML pipeline)

## Directory Structure

```
*.m                     MATLAB analysis modules (run from root)
citicorp_cfd/           OpenFOAM CFD case
  generators/           STL geometry generation scripts
  tools/                Analysis and utility scripts
  paraview/             ParaView visualization and presentation generators
  un_hq/                UN Headquarters building extraction demo
  dashboards/           Interactive HTML viewers
  docs/                 Technical documentation
presentations/          PowerPoint summary presentations
screenshots/            Development screenshots and reference photos
archive/                Preserved development iterations
```

## Prerequisites

- MATLAB R2020b+ (no toolbox dependencies)
- Python 3.8+ with `requests` (for STL generation)
- OpenFOAM Foundation v13 or ESI v2512 (for CFD, via WSL2)

## References

- Morgenstern, J. (1995). "The Fifty-Nine-Story Crisis." *The New Yorker*.
- NIST (2021). "Modern Reassessment of the Citicorp Building Design Wind Loads."
- ASCE 7-22: Minimum Design Loads and Associated Criteria for Buildings
- AISC 360-22: Specification for Structural Steel Buildings

## License

MIT License — see [LICENSE](LICENSE).
