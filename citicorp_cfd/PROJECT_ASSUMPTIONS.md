# Citicorp CFD Simulation — Assumptions & Limitations

**Project:** Wind load analysis for Citicorp Center structural crisis (1978)
**Date:** 2026-02-15
**Purpose:** Educational simulation for OR 750 (Reliability, Safety, and Risk)

---

## 1. Geometry Simplifications

### 1.1 Building Detail Level
- **LOD1 (Level of Detail 1)** buildings from NYC Open Data Footprints API
- Simple extruded footprints to `height_roof` elevation
- **No roof geometry:** All buildings have flat tops (no sloped roofs, crowns, setbacks)
- **No architectural details:** No cornices, parapets, balconies, or facade features

### 1.2 Citicorp Tower Specific
- ✗ **Missing iconic 45° slanted roof** (actual: 248m to 278.9m sloped top)
- ✓ Current model: Rectangular box with flat top at 278.9m
- ✓ 4-column stilt geometry preserved (correct topology)
- ✓ Building footprint and overall dimensions accurate (59.44m × 59.44m)
- ✗ No detailed facade (actual: reflective aluminum panels)
- ✗ No 45° rotation relative to Manhattan street grid (oriented N-S-E-W instead of NE-SW-NW-SE)

### 1.3 Surrounding Context
- 707 buildings from NYC Open Data within ~360m radius
- Accurate footprint polygons (L-shapes, setbacks preserved in plan view)
- **No street-level detail:** No trees, street furniture, awnings
- **No terrain elevation:** Ground assumed perfectly flat at z=0
- **Simplified domain:** 720m × 720m × 450m box (real urban context extends further)

---

## 2. Fluid Dynamics Modeling

### 2.1 Turbulence Model
- **RANS (Reynolds-Averaged Navier-Stokes)** with k-ω SST
- **Steady-state:** No time-accurate simulation of vortex shedding
- ✗ **Cannot capture unsteady effects:**
  - Vortex-Induced Vibration (VIV)
  - Cross-wind oscillations
  - Time-varying pressure fluctuations
- ✓ Good for: Mean wind loads, drag coefficient, time-averaged pressure field
- ⚠ **For unsteady effects, need LES/VLES** (Large Eddy Simulation / Very Large Eddy Simulation)

### 2.2 Computational Mesh
- **3.16 million cells** (snappyHexMesh with 3 refinement levels)
- **Near-surface resolution:** ~2m cell size at building walls
- ⚠ **Adequate but not highly refined** (research-grade: 10-30M cells)
- Mesh quality: 408 highly-skewed cells (max skew 3.65, mostly in stilt region)
- **No boundary layer mesh** (addLayers disabled for stability)
- ⚠ **Wall functions used** instead of resolving viscous sublayer (y+ ~30-300)

### 2.3 Boundary Conditions
- **Inlet:** Uniform profile (44.7 m/s constant), no atmospheric boundary layer profile
- **Outlet:** Zero gradient (zeroGradient for all fields)
- **Ground:** Slip wall (no ground friction) — real case would use roughWallFunction
- **Top/Sides:** Slip walls (symmetryPlane approximation)
- **Buildings:** No-slip walls (kqRWallFunction for k, omegaWallFunction for ω)

---

## 3. Structural Modeling

### 3.1 Fluid-Structure Interaction (FSI)
- ✗ **No FSI coupling** — structure is completely rigid
- ✗ **No building deformation** under wind load
- ✗ **No Tuned Mass Damper (TMD)** effects
- ⚠ Real Citicorp: 400-ton TMD on 63rd floor to reduce vibration
- ⚠ **One-way uncoupled:** CFD provides loads to MATLAB structural modules, but structure doesn't feed back to flow

### 3.2 What This Means
- CFD results show wind loads on **rigid** building
- Actual building would deflect (deflection → changes pressure field → changes loads)
- For most buildings: uncoupled approach is acceptable (deflections << building size)
- For slender/flexible towers: fully coupled FSI may be needed for accuracy

---

## 4. Computational Constraints

### 4.1 Hardware Limitations
- **Target:** 3-10M cells on laptop (Intel Ultra 9 185H, 64GB RAM, WSL2)
- **Current run:** 3.16M cells, 10 cores, ~4 hours runtime
- ⚠ **Research-grade simulations:** 15-100M cells, HPC clusters, days-to-weeks

### 4.2 Single Wind Direction
- ✓ **Simulated:** East wind (+X direction, 44.7 m/s = 100 mph = Cat 2 hurricane)
- ✗ **Not simulated:** Multiple wind angles (real DBD studies: 0°-360° in 10° increments)
- ✗ **Not simulated:** Wind speed sensitivity (real: 10-100 m/s parametric study)
- ⚠ Real DBD (Database-Assisted Design): ~1000+ wind tunnel tests, interpolated database

---

## 5. Physical Simplifications

### 5.1 Atmospheric Conditions
- **Uniform wind:** Constant 44.7 m/s at all heights (no wind shear)
- ⚠ Real atmosphere: Logarithmic or power-law velocity profile
- **No thermal effects:** Isothermal flow (no buoyancy, no temperature stratification)
- **No humidity, precipitation, or phase change**
- **Air properties:** ρ = 1.225 kg/m³, ν = 1.5e-5 m²/s (sea level, 15°C)

### 5.2 Time Scales
- **Steady-state RANS:** No temporal resolution
- ✗ **Cannot predict:** Frequency of vortex shedding, dynamic amplification, resonance
- ⚠ Real Citicorp crisis: Across-wind forces at 0.15-0.2 Hz matched building natural frequency
- ⚠ Strouhal number St ≈ 0.15 for square cylinders → shedding at f = St × U / D ≈ 0.11 Hz (close to resonance)

---

## 6. Validation & Verification

### 6.1 What Was Validated
- ✓ **Drag coefficient Cd = 3.24** (tower+stilts), 3.57 (tower only)
- ✓ **Hand calculation:** Empirical Cd ≈ 1.4 (isolated) × 1.7 (urban) × 1.25 (stilt) ≈ 3.0
- ✓ **Agreement within 8%** — reasonable for educational simulation

### 6.2 What Was NOT Validated
- ✗ No comparison to wind tunnel data (none publicly available for Citicorp)
- ✗ No mesh independence study (would need 6M, 12M, 24M cell runs)
- ✗ No turbulence model sensitivity (k-ε, Realizable k-ε, LES comparison)
- ✗ No experimental pressure coefficient (Cp) distribution on building faces

---

## 7. Data Sources & Provenance

### 7.1 Building Geometry
- **Source:** NYC Open Data — Building Footprints (DoITT, 2023 update)
- **API:** Socrata GeoJSON endpoint, 707 buildings within 40.7507-40.7651°N, 73.9618-73.9762°W
- **Accuracy:** Footprints ±1m, heights from photogrammetric survey (±1-3m)
- ⚠ **No semantic labels:** Cannot distinguish residential/commercial/historic

### 7.2 Historical Context
- **Citicorp dimensions:** Morgenstern (1995), LeMessurier (1995)
- **Wind speed:** 100 mph ≈ 44.7 m/s (Cat 2 hurricane, 1-in-16-year NYC event)
- **Structural crisis:** Quartering wind (45°) on bolted chevron bracing
- ⚠ **Our simulation:** Cardinal wind (90°), no 45° rotation

---

## 8. Known Discrepancies from Reality

| Aspect | Real Citicorp (1978) | Our Simulation |
|--------|---------------------|----------------|
| **Roof geometry** | 45° sloped wedge (iconic) | Flat rectangular box |
| **Orientation** | 45° rotated (NE corner facing) | Aligned to N-S-E-W grid |
| **Wind direction** | Quartering (45°) most critical | East wind (90°) only |
| **Bracing** | Chevron V-bracing (apex down) | Not modeled (rigid walls) |
| **Connections** | Welded (repaired) | N/A |
| **TMD** | 400-ton concrete block, 63rd floor | Not modeled |
| **Surrounding buildings** | 1978 skyline (shorter neighbors) | 2023 skyline (some taller) |
| **Dynamic response** | Across-wind flutter, resonance | Static rigid structure |
| **Simulation type** | Wind tunnel (database-assisted) | CFD RANS steady-state |

---

## 9. Future Improvements (Priority Order)

### High Priority
1. **Add slanted roof to tower** (TUM NYC LOD2 CityGML or manual geometry)
2. **Rotate building 45°** to match historical orientation
3. **Simulate 45° quartering wind** (critical case from Morgenstern 1995)
4. **Mesh refinement study** (6M, 12M cells on Hopper HPC)

### Medium Priority
5. **Atmospheric boundary layer profile** (log-law inlet BC)
6. **Ground roughness** (urban terrain, z₀ ≈ 1m)
7. **Multiple wind directions** (0°-360° in 30° steps)
8. **Boundary layer mesh** (addLayers, resolve y+ < 5)

### Low Priority (Research Scope)
9. **LES/VLES for unsteady loads** (24M+ cells, weeks runtime)
10. **Coupled FSI** (OpenFOAM + preCICE + CalculiX or MATLAB)
11. **TMD implementation** (mass-spring-damper in structural model)
12. **Full parametric study** (wind speed, direction, building stiffness)

---

## 10. Summary: What This Simulation IS and ISN'T

### ✓ This Simulation IS:
- Educational demonstration of CFD workflow for urban wind engineering
- Reasonable estimate of **mean drag force** on Citicorp tower
- Valid for **qualitative flow visualization** (wake patterns, stagnation regions)
- Proof-of-concept for open-source CFD tools (OpenFOAM, ParaView, Python)
- **Correct order of magnitude** for wind loads (Cd ≈ 3.2 validated)

### ✗ This Simulation IS NOT:
- Research-grade wind engineering analysis (would need 10-100M cells, LES, mesh study)
- Replacement for wind tunnel testing (no dynamic effects, no unsteady loads)
- Accurate representation of 1978 historical crisis (wrong orientation, static structure)
- Suitable for design decisions (no validation against experiment or high-fidelity data)
- Capable of predicting across-wind vibration or resonance (steady RANS limitation)

---

## Notes for Meeting Discussion

**Strengths to emphasize:**
- Comprehensive workflow (geometry → meshing → solving → post-processing)
- Validated Cd against empirical correlations (3.24 vs 3.0, within 8%)
- Stable 3000-iteration convergence after debugging initial crash
- Integration with MATLAB structural modules (uncoupled one-way FSI)

**Limitations to acknowledge:**
- LOD1 geometry (no slanted roof, simplified buildings)
- Steady RANS (cannot capture VIV or unsteady effects)
- Single wind direction (not parametric study)
- No experimental validation data available for Citicorp

**If asked about improvements:**
- Refer to Section 9 (Future Improvements)
- TUM LOD2 dataset for better geometry (+1-3 days setup)
- GMU Hopper HPC for larger mesh (8-15M cells feasible)
- LES would require 10-100× computational cost (weeks on HPC)

**Key message:**
"This is an educational simulation demonstrating CFD methodology for urban wind loads. The mean drag coefficient matches empirical estimates, but the model simplifies geometry (no slanted roof), uses steady-state turbulence (no vortex shedding), and assumes rigid structure (no FSI). For research-grade results, we'd need LOD2 geometry, mesh refinement study, LES turbulence model, and experimental validation — all feasible as Phase 2 on GMU's Hopper cluster."
