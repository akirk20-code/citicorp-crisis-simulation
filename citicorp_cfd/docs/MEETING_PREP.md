# CFD Meeting Prep — Files to Show/Share

## Quick Demo Path

1. **Open the dashboard** (main showpiece):
   - File: `cfd_dashboard.html` (open in Chrome/Edge)
   - Shows: 3D building mesh, live convergence plot, all solver parameters
   - Interactive: drag to orbit, scroll to zoom, camera presets

2. **Key talking points from the dashboard**:
   - 707 real NYC buildings from Open Data API (photogrammetric heights)
   - 3.16M cell mesh, k-omega SST, SIMPLE steady-state
   - 44.7 m/s (100 mph) atmospheric boundary layer at z=10m
   - ~75.8 m/s at roof height (248m) via log-law profile

---

## Files to Share (by topic)

### Geometry & Data Pipeline
| File | What it shows |
|------|--------------|
| `generate_stl.py` | NYC Open Data API query, coord transform, STL extrusion |
| `constant/triSurface/*.stl` | Binary STL files (tower, stilts, 707 surroundings) |

### OpenFOAM Configuration
| File | What it shows |
|------|--------------|
| `system/controlDict` | Solver setup: foamRun, incompressibleFluid, 3000 iters |
| `system/fvSchemes` | Discretization: upwind (stability), limited corrected |
| `system/fvSolution` | GAMG pressure, smoothSolver U/k/omega, relaxation factors |
| `system/snappyHexMeshDict` | 5-level refinement, Foundation syntax |
| `system/blockMeshDict` | 48x48x40 background mesh (15m base cells) |
| `0/U`, `0/k`, `0/omega`, `0/p` | Boundary conditions (ABL inlet, wall functions) |

### Run Script
| File | What it shows |
|------|--------------|
| `Allrun` | Full pipeline: STL gen -> meshing -> decompose -> solve |

---

## Key Stability Story (good for discussion)

1. **Initial attempt**: Uniform 44.7 m/s field -> SIGFPE crash at iter 157
   - Continuity errors: 0.03 -> 70 -> 7,225 -> 1.8M -> CRASH
   - Root cause: GAMG pressure solver hitting 408 highly-skew cells

2. **Fix applied**: potentialFoam warm start
   - Creates divergence-free velocity field before RANS solve
   - Combined with: lower relaxation (U=0.2, p=0.15), nNonOrthCorrectors=2
   - Result: stable convergence, Ux residual 0.014 -> 0.0005 by iter 300

3. **Mesh quality**:
   - Max non-orthogonality: 64.75 deg (OK, threshold is 70)
   - Max skewness: 7.45 (408 faces, 0.01% of mesh)
   - renumberMesh applied for cache optimization

---

## Solver Current Status

- Running on 10 cores (Intel Ultra 9 185H, WSL2 Ubuntu)
- ~5s per iteration
- Target: 3000 iterations (~4 hours total)
- At meeting time: expect ~800-900 iterations completed
- Residuals: monotonically decreasing (Ux ~ 5e-4, p ~ 0.04)

---

## GitHub Repo
https://github.com/akirk20-code/citicorp-crisis-simulation

## Project Context
- OR 750: Reliability, Safety, and Risk Analysis (GMU)
- Citicorp Center 1978 wind crisis: bolted vs welded chevron bracing
- 6 MATLAB modules (Monte Carlo, FEA, 3D viz) + OpenFOAM CFD
