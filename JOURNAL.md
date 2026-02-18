# Project Journal — Citicorp Center CFD Simulation
**Course:** OR 750 — Reliability, Safety, and Risk
**Repository:** https://github.com/akirk20-code/citicorp-crisis-simulation
**Last updated:** 2025-02-15

---

## Project Overview

Simulation of wind loads on the Citicorp Center (now 601 Lexington Ave, NYC) to study the
1978 structural crisis. The building is famous for an engineering near-disaster: after
construction, student Diane Hartman discovered that a design change (bolted instead of welded
Chevron bracing joints) made the building dangerously vulnerable to quartering winds. Emergency
repairs were secretly performed overnight to prevent collapse.

**Six modules in MATLAB + OpenFOAM CFD:**
1. Monte Carlo structural reliability (16 scenarios)
2. Transfer zone finite element analysis (FEA)
3. Chevron brace analysis
4. Structural elevation visualization
5. Interactive 3D model
6. CFD wind simulation pipeline

---

## Session Log

### Session 1 — Initial Setup (pre-journal)

**Accomplished:**
- Created 6 MATLAB modules for structural analysis
- Set up OpenFOAM CFD case (`citicorp_cfd/`) with Foundation v13 on WSL2 Ubuntu 24.04
- Implemented ATM boundary layer inlet (log-law ABL profile, z₀ = 0.1 m, U_ref = 44.7 m/s at z=10m)
- First CFD run with simplified 18-building hand-drawn surroundings (LOD0)
- Computed drag coefficient: Cd ≈ 3.24 (vs hand-calc Cd ≈ 3.03, ~6.5% difference)
- Created ParaView post-processing pipeline and CFD dashboard (HTML)

**Key parameters set:**
- Domain: 2200 × 1800 × 420 m (streamwise × transverse × vertical)
- Wind direction: +X (west-to-east, matching prevailing NYC wind)
- Reference velocity: 44.7 m/s (100 mph, design wind speed for 1978 era)
- Turbulence model: k-ω SST (steady-state RANS)
- Mesh: ~3.16M cells

---

### Session 2 — Documentation, Data Sources, and 6-Case Mesh Setup
**Date:** ~2025-02-14 to 2025-02-15

#### What was done

**Documentation created:**
| File | Contents |
|------|----------|
| `PROJECT_ASSUMPTIONS.md` | All simplifications in the CFD model |
| `IMPROVEMENT_ROADMAP.md` | Technical path to fix each assumption (LES, FSI, etc.) |
| `HOPPER_SETUP.md` | GMU HPC (Hopper) Slurm job scripts for RANS/LES/GPU runs |
| `MESH_SIZING_GUIDE.md` | Comparison to Camelli et al. (2006) MSG study — 50-80M cells |
| `CD_VALIDATION_HAND_CALC.md` | Hand calculation of Cd = 1.35 × 1.7 × 1.2 × 1.1 = 3.03 |
| `VELOCITY_ANALYSIS.md` | Discovery that ABL inlet was already implemented |
| `ACTUAL_BUILDING_HEIGHTS.md` | Real NYC API data: 1,665 buildings, median 18 m, max 425.5 m |
| `HISTORICAL_1978_DATA.md` | 1978 filter: 1,472 buildings, mean 29.4 m, max 288.6 m (Chrysler) |
| `BUILDING_HEIGHT_ANALYSIS.md` | Summary of height data with context |

**Key technical discovery:**
> The 0/U boundary condition already uses `atmBoundaryLayerInletVelocity` (log-law ABL profile
> with z₀ = 0.1 m). The "uniform inlet" assumption in PROJECT_ASSUMPTIONS.md was incorrect.
> At roof height (278.9 m), U ≈ 77.9 m/s — much higher than the 44.7 m/s reference.
> However, both CFD and hand-calc normalize by U_ref = 44.7 m/s (wind engineering convention),
> so the Cd comparison remains valid.

**Professor's question about Cd interference:**
> "Should the Cd not be lower due to interference/blockage effects?"
>
> Answer: For typical clustered buildings, sheltering lowers Cd. But Citicorp is atypical:
> the median neighbor height is only 18 m (6.5% of Citicorp's 278.9 m). The upper 230 m is
> fully exposed — no sheltering. The dense low-rise cluster actually increases turbulence
> intensity in the lower wake, slightly increasing Cd. The domain blockage ratio (1.3H
> clearance vs. recommended 5-10H) adds an artificial +5-15%. Corrected Cd ≈ 3.24/1.10 ≈
> 2.95, which agrees with the hand-calc of 3.03.

**Data sources compared:**
| Source | Buildings | Mean height | Notes |
|--------|-----------|-------------|-------|
| Original (hand-drawn) | 18 | ~25m est. | LOD0, no real data |
| NYC Footprints API 2024 | 1,665 | 36.0 m | Used in nyc3d, lod2 cases |
| NYC Footprints API 1978 | 1,472 | 29.4 m | Uses `construction_year <= 1978` filter |
| OSM Overpass | N/A | N/A | 504 timeout — surroundings empty |

**Six geometry scripts created:**
| Script | Surroundings | Triangles |
|--------|-------------|-----------|
| `generate_stl_osm.py` | OSM Overpass (timed out → 0 buildings) | 20+48+0 |
| `generate_stl_nyc3d.py` | NYC API 2024 | 22+48+28,376 |
| `generate_stl_lod2.py` | NYC API 2024 + roof inference | 12+48+32,080 |
| `generate_stl_osm_1978.py` | NYC API pre-1978 | 19+40+57,948 |
| `generate_stl_nyc3d_1978.py` | NYC API pre-1978 | 19+40+57,948 |
| `generate_stl_lod2_1978.py` | NYC API pre-1978 + roof types | 19+40+69,001 |

The LOD2 script infers roof types by construction era:
- mansard roofs: buildings pre-1900 (beaux-arts, Gilded Age)
- hipped roofs: buildings pre-1940 (early 20c brick)
- gabled roofs: buildings pre-1960
- flat roofs: modern high-rises

**WSL case directories created (all in /home/kirka/):**
- `citicorp_cfd_osm`
- `citicorp_cfd_nyc3d`
- `citicorp_cfd_lod2`
- `citicorp_cfd_osm_1978`
- `citicorp_cfd_nyc3d_1978`
- `citicorp_cfd_lod2_1978`

Each has full OpenFOAM structure: `0/`, `constant/`, `system/`, `.foam` file.

**Master mesh script:** `mesh_all_cases.sh`
Runs all 6: STL generation → copy → surfaceFeatureExtract → blockMesh → snappyHexMesh (8 cores)

**Topographic height map tool:** `plot_topo_map.py`
Rasterizes STL triangles onto a 5m grid and plots urban canopy height as a colored map.
Output: `topo_map.png` + `topo_grid.csv`

**Test run on 1978 dataset:**
```
Domain:    X [-656, 657] m,  Y [-844, 849] m
Max height: 288.6 m (Chrysler Building)
Buildings:  30.1% of domain coverage
Mean height (buildings): 69.5 m
p95:        192.1 m
```

#### Issues encountered and fixed
| Issue | Fix |
|-------|-----|
| Bash tool runs in Git Bash (MINGW64), NOT WSL | Use `/c/` paths, not `/mnt/c/` |
| Windows Python not found as `python3` | Use full path `/c/Users/kirka/.../python.exe` |
| Unicode chars (`→`, `✓`) crash cp1252 terminal | Replace with ASCII equivalents (`->`, `+`) |
| OSM Overpass 504 timeout | Surroundings STL is empty; use NYC API instead |
| `mesh_all_cases.sh` had `/c/` paths (Git Bash) | Fixed to `/mnt/c/` for WSL |
| `2>nul` creates literal file on Windows (OneDrive renames to `_ul`) | Use `2>/dev/null` |

---

### Session 3 — [Next session will be added here]

**Planned:**
- [ ] Run `mesh_all_cases.sh` in WSL to generate 6 OpenFOAM meshes
- [ ] Open `.foam` files in ParaView to compare geometries
- [ ] Run simpleFoam for 2024 and 1978 cases, compare Cd
- [ ] Generate topographic maps for all 6 cases
- [ ] Compare with literature Cd values

---

## Current State

### What exists
- 6 MATLAB structural analysis modules (complete, working)
- 1 OpenFOAM CFD case `citicorp_cfd/` (original, meshed, solved)
- 6 new case directories in WSL (structured, not yet meshed)
- 6 STL geometry scripts (all tested, all produce output)
- Current STLs in `constant/triSurface/` = last run = **1978 nyc3d/osm dataset** (57,948 surrounding triangles)

### What's next
1. **Run meshes** — WSL terminal: `bash "/mnt/c/.../mesh_all_cases.sh"`
   Expected: 8-10M cells per case, ~4-8h per wind direction on this hardware
2. **Visualize** — Open `.foam` files in ParaView (via `\\wsl$\Ubuntu\home\kirka\...`)
3. **Compare Cd** — 2024 vs. 1978 skyline, 3 data sources each

---

## Key Physical Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| Citicorp height | 278.9 m | As-built drawings |
| Tower plan width | 47.85 m | As-built drawings |
| Stilt height | 34.75 m | Morgenstern 1995 |
| Stilt width | 7.32 m | Morgenstern 1995 |
| Reference wind speed | 44.7 m/s (100 mph) | 1978 NYC design code |
| ABL roughness z₀ | 0.1 m (current) / 1.0 m (urban target) | ASCE 7 terrain Category D |
| Air density | 1.225 kg/m³ | ISA sea level |
| Wind direction | +X = westerly (prevailing NYC wind) | NYC Wind Atlas |

---

## Drag Coefficient Results

| Case | Cd | Notes |
|------|----|-------|
| Hand calculation | 3.03 | Cd_shape × aspect × gust × ABL factor |
| CFD (original LOD0) | 3.24 | 3.16M cells, k-ω SST, ABL inlet, z₀=0.1m |
| Difference | +6.5% | Partly explained by domain blockage (+5-10%) |
| Corrected CFD estimate | ~2.95 | After 1.10× blockage correction |

> The hand-calc and CFD agree within ~3% after blockage correction — good validation.

---

## Architecture Decisions

### Why NYC Open Data API (not TUM CityGML)
- TUM download: 2.4 GB vs. NYC API: ~2 MB per query
- TUM requires `cjio`, `lxml`, `citygml2stl` (complex setup)
- NYC API: `requests` only (stdlib-friendly)
- LOD2 roof detail adds 1-3 days parsing effort for marginal RANS benefit
- Roof shapes negligible for mean wind loads in steady-state RANS

### Why Foundation v13 (not ESI v2512)
- Foundation v13: simpler setup, well-documented for incompressible flows
- ESI v2512 with GPU (AmgX/PETSc): available but adds complexity
- GPU (RTX 4070, 8GB VRAM) only beneficial for meshes ≥ 500k cells — current mesh qualifies
- Decision: use Foundation v13 for quick iterations, switch to ESI+GPU for production runs

### Why k-ω SST (not LES)
- LES requires 50-80M cells and 128+ cores for urban flows (Camelli et al., 2006)
- k-ω SST: good wall-bounded flow prediction, robust convergence
- For a class project, RANS is appropriate; LES documented as future work in IMPROVEMENT_ROADMAP.md

---

## References

- Morgenstern, J. (1995). The fifty-nine-story crisis. *The New Yorker*.
- Camelli, F. et al. (2006). VLES study of the Madison Square Garden area winds. AIAA Paper 2006-1384.
- ASCE 7-22 (Minimum Design Loads for Buildings and Other Structures)
- Citicorp Center Wikipedia / SkyscraperPage for architectural dimensions
- NYC Open Data Building Footprints API: https://data.cityofnewyork.us/resource/5zhs-2jue.geojson
