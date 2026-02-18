# CFD Improvement Roadmap — Technical Deep Dive

**Purpose:** Detailed technical plan for addressing each limitation in PROJECT_ASSUMPTIONS.md
**Audience:** Graduate-level understanding of CFD and structural dynamics
**Scope:** From quick fixes (hours) to research-grade enhancements (weeks)

---

## Table of Contents
1. [Geometry Improvements](#1-geometry-improvements)
2. [Turbulence Model Upgrades](#2-turbulence-model-upgrades)
3. [Mesh Quality & Resolution](#3-mesh-quality--resolution)
4. [Boundary Condition Realism](#4-boundary-condition-realism)
5. [Fluid-Structure Interaction](#5-fluid-structure-interaction)
6. [Computational Infrastructure](#6-computational-infrastructure)
7. [Validation & Verification](#7-validation--verification)

---

## 1. Geometry Improvements

### 1.1 Add Citicorp Slanted Roof (Priority: HIGH)

**Current Issue:** Tower is a rectangular box with flat top at 278.9m
**Target:** 45° sloped wedge from 248m (floor 60) to 278.9m peak

#### Option A: Manual Geometry in generate_stl.py (Recommended)
**Effort:** 2-4 hours
**Approach:**
```python
# In generate_stl.py, replace box_triangles() for tower with:

def tower_with_sloped_roof(half_width, base_z, roof_start_z, peak_z):
    """Generate Citicorp tower: box base + 45° sloped wedge top."""
    hw = half_width
    tris = []

    # 1. Box portion (base to roof_start_z = 248m)
    # Bottom face (z = base_z = 24.4m stilt top)
    tris.extend([
        [0,0,-1, -hw,-hw,base_z, hw,-hw,base_z, hw,hw,base_z],
        [0,0,-1, -hw,-hw,base_z, hw,hw,base_z, -hw,hw,base_z],
    ])

    # Four vertical walls
    for (x1,y1,x2,y2) in [(-hw,-hw, hw,-hw), (hw,-hw, hw,hw),
                           (hw,hw, -hw,hw), (-hw,hw, -hw,-hw)]:
        nx, ny = -(y2-y1), (x2-x1)  # Outward normal
        norm = (nx**2 + ny**2)**0.5
        nx, ny = nx/norm, ny/norm
        tris.extend([
            [nx,ny,0, x1,y1,base_z, x2,y2,base_z, x2,y2,roof_start_z],
            [nx,ny,0, x1,y1,base_z, x2,y2,roof_start_z, x1,y1,roof_start_z],
        ])

    # 2. Sloped roof (248m to 278.9m peak)
    # Peak point at center
    peak = [0, 0, peak_z]

    # Four sloped triangular faces
    corners = [(-hw,-hw,roof_start_z), (hw,-hw,roof_start_z),
               (hw,hw,roof_start_z), (-hw,hw,roof_start_z)]

    for i in range(4):
        c1 = corners[i]
        c2 = corners[(i+1)%4]
        # Two triangles per face (split quad)
        v1 = [c2[j]-c1[j] for j in range(3)]
        v2 = [peak[j]-c1[j] for j in range(3)]
        # Normal = v1 × v2
        nx = v1[1]*v2[2] - v1[2]*v2[1]
        ny = v1[2]*v2[0] - v1[0]*v2[2]
        nz = v1[0]*v2[1] - v1[1]*v2[0]
        norm = (nx**2+ny**2+nz**2)**0.5
        nx,ny,nz = nx/norm, ny/norm, nz/norm

        tris.append([nx,ny,nz, *c1, *c2, *peak])

    return tris

# Replace line ~430:
tower_tris = tower_with_sloped_roof(HALF_T, STILT_H, 248.0, TOWER_H_TOP)
```

**Testing:**
1. Run `python3 generate_stl.py`
2. Open `citicorp_tower.stl` in ParaView or MeshLab
3. Verify sloped roof appears correctly
4. Check triangle count increased by ~8-12 triangles (negligible mesh impact)

**Mesh Implications:**
- Surface feature refinement will capture slope automatically
- May need `resolveFeatureAngle 30;` in snappyHexMeshDict (currently 30°, adequate)
- No significant cell count increase (slope is smooth feature)

---

### 1.2 Rotate Building 45° to Historical Orientation (Priority: HIGH)

**Current Issue:** Tower aligned N-S-E-W (cardinal directions)
**Historical Reality:** Rotated 45° to street grid (NE corner prominent)
**Why It Matters:** Quartering wind (45°) was the critical load case in 1978 crisis

#### Implementation in generate_stl.py:
```python
import numpy as np

def rotate_z(points, angle_deg):
    """Rotate points around Z-axis by angle_deg."""
    theta = np.radians(angle_deg)
    R = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta),  np.cos(theta), 0],
        [0, 0, 1]
    ])
    return points @ R.T

# After generating tower geometry:
tower_verts = np.array([[tri[3+i*3:6+i*3] for i in range(3)]
                        for tri in tower_tris]).reshape(-1, 3)
tower_verts = rotate_z(tower_verts, 45.0)  # Rotate 45°

# Rebuild triangles with rotated vertices and recalculate normals
```

**Better Approach:** Rotate in OpenFOAM `surfaceTransformPoints`:
```bash
surfaceTransformPoints -rotateAlongVector '((0 0 0)(0 0 1) 45)' \
    constant/triSurface/citicorp_tower.stl \
    constant/triSurface/citicorp_tower_rotated.stl
```

**Wind Direction Setup:**
After rotating building 45°, simulate multiple wind directions:
- **0° (East):** Faces flat face (current case)
- **45° (NE):** Quartering wind, hits corner (critical case)
- **90° (North):** Faces flat face
- **135° (NW):** Quartering wind, opposite corner

**Effort:** 1 hour for geometry rotation, 4-6 hours per wind direction simulation

---

### 1.3 Upgrade to LOD2 Geometry (Priority: MEDIUM)

**Current Issue:** Simple extruded footprints, no roof detail
**Target:** TUM NYC LOD2 CityGML (1M+ buildings with roof shapes)

#### Data Source:
- **URL:** https://3d.bk.tudelft.nl/opendata/3dbag/download.html (NYC subset)
- **Format:** CityGML 2.0 or CityJSON
- **Size:** ~2.4 GB for all NYC (can download by borough)
- **Content:** LoD2.2 — roof shapes, dormers, setbacks, chimneys

#### Conversion Pipeline:
```bash
# 1. Download Manhattan buildings (CityJSON format)
wget https://3d.bk.tudelft.nl/opendata/3dbag/v2/manhattan.city.json

# 2. Install cjio (CityJSON tools)
pip install cjio

# 3. Extract bounding box around Citicorp
cjio manhattan.city.json subset --bbox 40.7507 -73.9762 40.7651 -73.9618 \
    save manhattan_citicorp_area.city.json

# 4. Convert to OBJ (then to STL)
cjio manhattan_citicorp_area.city.json export \
    --format obj \
    --split-buildings \
    manhattan_buildings/

# 5. Merge OBJ files to STL using meshio
python3 << 'EOF'
import meshio
import glob
import numpy as np

all_points = []
all_cells = []
offset = 0

for obj_file in glob.glob('manhattan_buildings/*.obj'):
    mesh = meshio.read(obj_file)
    all_points.append(mesh.points)
    cells = mesh.cells_dict['triangle'] + offset
    all_cells.append(cells)
    offset += len(mesh.points)

merged = meshio.Mesh(
    points=np.vstack(all_points),
    cells=[('triangle', np.vstack(all_cells))]
)
meshio.write('surroundings_lod2.stl', merged)
EOF
```

**Challenges:**
- **Coordinate transform:** CityJSON uses RD New (EPSG:28992) → need pyproj for WGS84 → local meters
- **Semantic surfaces:** CityGML labels walls/roofs/ground — need to filter only walls/roofs
- **Data size:** 2.4 GB → 500+ MB STL (manageable but large)
- **Parsing complexity:** CityJSON is nested JSON, requires careful schema navigation

**Effort:** 1-3 days (data download, parsing, conversion, validation)

**Payoff:**
- Realistic roof shapes (gabled, hipped, flat with parapets)
- Building setbacks (step pyramids like Empire State Building)
- Better wake flow prediction (roof geometry affects separation)

---

## 2. Turbulence Model Upgrades

### 2.1 Large Eddy Simulation (LES) for Unsteady Effects (Priority: MEDIUM-HIGH)

**Current Issue:** RANS k-ω SST time-averages turbulence, cannot predict:
- Vortex shedding frequency
- Cross-wind oscillating loads
- Pressure fluctuations (RMS Cp)
- Dynamic amplification near resonance

**Target:** LES with Smagorinsky or WALE subgrid model

#### Key Differences: RANS vs LES

| Aspect | RANS (k-ω SST) | LES (Smagorinsky) |
|--------|---------------|-------------------|
| **Resolved scales** | None (all turbulence modeled) | Large eddies resolved, small modeled |
| **Time** | Steady-state | Transient (10+ flow-through times) |
| **Mesh** | 3-10M cells | 20-100M cells (Δx < 0.1 × building width) |
| **Time step** | N/A | Δt < 0.01 sec (CFL < 1) |
| **Runtime** | 4 hours (3000 SIMPLE iterations) | 200-1000 hours (10,000-100,000 timesteps) |
| **Output** | Mean fields (U, p, k, ω) | Time series (U(t), p(t), forces(t)) |

#### OpenFOAM LES Setup:

**1. Change solver:** `foamRun` → `pimpleFoam` (transient incompressible)

**2. Update turbulence model in `constant/turbulenceProperties`:**
```cpp
simulationType      LES;

LES
{
    model           Smagorinsky;
    turbulence      on;
    printCoeffs     on;

    SmagorinskyCoeffs
    {
        Ck              0.094;
        Ce              1.048;
    }
}
```

**3. Refine mesh for LES resolution:**
```
Target: Δx ≈ 0.1 × building width = 0.1 × 59.44m ≈ 6m
Current: Δx ≈ 2m near surfaces → adequate!
But need: Δx ≈ 1-2m throughout wake region → 20-30M cells

snappyHexMeshDict:
- Increase refinementRegions/wakeRegion from level 1 to level 2-3
- Add refinementSurfaces for building to level 4-5 (0.5-1m)
```

**4. Set transient time controls in `system/controlDict`:**
```cpp
application     pimpleFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         100;        // 100 seconds = ~8 flow-through times at U=44.7 m/s, L=720m
deltaT          0.005;      // 5 ms timestep (CFL ≈ 0.5-1.0)
writeControl    adjustableRunTime;
writeInterval   1.0;        // Write every 1 second (200 timesteps)
```

**5. Add PIMPLE controls in `system/fvSolution`:**
```cpp
PIMPLE
{
    nOuterCorrectors    2;
    nCorrectors         3;
    nNonOrthogonalCorrectors 2;

    residualControl
    {
        U      1e-6;
        p      1e-5;
    }
}
```

**6. Force output for time-history:**
```cpp
functions
{
    forces
    {
        type            forces;
        libs            ("libforces.so");
        patches         (citicorp_tower citicorp_stilts);
        rho             rhoInf;
        rhoInf          1.225;
        CofR            (0 0 124);
        writeControl    timeStep;
        writeInterval   1;      // Write every timestep
    }
}
```

#### Post-Processing LES Results:

**Time-averaged fields:**
```bash
# Average last 50 seconds (quasi-steady regime)
postProcess -func 'timeAverage(U,p)' -time '50:100'
```

**FFT for vortex shedding frequency:**
```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch

# Load force time-history
t, Fx, Fy, Fz = np.loadtxt('postProcessing/forces/0/force.dat', unpack=True)

# FFT of cross-wind force (Fy)
freqs, psd = welch(Fy, fs=1/0.005, nperseg=2048)

# Find dominant frequency
f_peak = freqs[np.argmax(psd)]
St = f_peak * 59.44 / 44.7  # Strouhal number = f × D / U

print(f"Shedding frequency: {f_peak:.3f} Hz")
print(f"Strouhal number: {St:.3f}")
print(f"Expected for square cylinder: St ≈ 0.12-0.15")
```

**Effort:** 1 week setup + 1-4 weeks runtime on Hopper HPC
**Hardware:** 128-core node, 30M cells, 1-2 weeks walltime

---

### 2.2 VLES (Very Large Eddy Simulation) — Hybrid Approach

**Compromise:** RANS near walls, LES in free shear layers
**Model:** `kOmegaSSTDES` (Detached Eddy Simulation)

```cpp
// In constant/turbulenceProperties
simulationType      LES;

LES
{
    model           kOmegaSSTDES;
    turbulence      on;

    kOmegaSSTDESCoeffs
    {
        CDES            0.61;   // DES constant
    }
}
```

**Advantage:** Coarser mesh OK (10-15M cells), shorter runtime (2-5 days)
**Tradeoff:** Less accurate than pure LES in boundary layer, but 5-10× cheaper

---

## 3. Mesh Quality & Resolution

### 3.1 Boundary Layer Mesh with addLayers (Priority: MEDIUM)

**Current Issue:** Wall functions (y+ = 30-300), no viscous sublayer resolution
**Target:** Resolve boundary layer (y+ < 1) or buffer layer (y+ = 1-5)

#### Enable addLayers in snappyHexMeshDict:
```cpp
addLayers
{
    relativeSizes true;

    layers
    {
        "citicorp_.*|surroundings"
        {
            nSurfaceLayers 5;           // 5 prism layers
            expansionRatio 1.2;          // 20% growth per layer
            finalLayerThickness 0.3;     // 30% of base cell
            minThickness 0.1;            // Don't go below 10% of cell
        }
    }

    // Quality controls (critical for stability)
    featureAngle 180;                    // Don't break at sharp edges
    nRelaxIter 5;                        // Iterations to fit layers
    nSmoothSurfaceNormals 1;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}
```

#### Why We Disabled It Initially:
- **Stability:** Prism layers at sharp corners (stilt edges) can create high-aspect-ratio cells → solver crash
- **Complexity:** Requires careful tuning of quality parameters
- **RANS:** Wall functions are adequate for mean loads in RANS

#### When You NEED It:
- **LES:** Must resolve or nearly resolve boundary layer (y+ < 5)
- **Heat transfer:** Thermal boundary layer resolution
- **Separation accuracy:** Better prediction of separation point on rounded edges

**Testing Strategy:**
1. Enable only on tower (smooth surface), disable on stilts
2. Run checkMesh — verify maxAspectRatio < 100, maxSkewness < 4
3. Start with 3 layers before going to 5
4. Use potentialFoam initialization (mandatory with layers)

**Effort:** 1-2 days tuning, frequent checkMesh failures expected

---

### 3.2 Mesh Independence Study (Priority: HIGH for publication)

**Current Mesh:** 3.16M cells
**Target:** Prove results converge with refinement

#### Systematic Refinement:
| Mesh | Base Res | Refinement Levels | Cell Count | Δx (tower) | Runtime |
|------|----------|-------------------|------------|------------|---------|
| Coarse | 20m | 2 | 1.5M | 5m | 2 hrs |
| **Medium (current)** | 15m | 3 | 3.16M | 2m | 4 hrs |
| Fine | 12m | 4 | 8M | 1m | 12 hrs |
| Extra-fine | 10m | 5 | 20M | 0.5m | 40 hrs |

#### Convergence Metrics:
```python
# Calculate Grid Convergence Index (GCI)
Cd_coarse = 3.45
Cd_medium = 3.24
Cd_fine   = 3.18

r = 2.0  # Refinement ratio (Δx_coarse / Δx_medium)
p = np.log((Cd_fine - Cd_medium) / (Cd_medium - Cd_coarse)) / np.log(r)
phi_ext = (r**p * Cd_fine - Cd_medium) / (r**p - 1)  # Richardson extrapolation

GCI = 1.25 * abs((Cd_fine - Cd_medium) / Cd_medium) / (r**p - 1)
print(f"Order of convergence p = {p:.2f}")
print(f"Extrapolated Cd = {phi_ext:.3f}")
print(f"GCI (fine) = {GCI*100:.2f}%")
```

**Acceptance:** GCI < 3% for "grid-converged" results

**Effort:** 4 mesh runs × 4-40 hours = 2-7 days on Hopper

---

## 4. Boundary Condition Realism

### 4.1 Atmospheric Boundary Layer Inlet (Priority: HIGH)

**Current Issue:** Uniform inlet (U = 44.7 m/s constant)
**Reality:** Logarithmic or power-law velocity profile

#### Power-Law Profile (Engineering Standard):
```
U(z) = U_ref × (z / z_ref)^α

Where:
- U_ref = 44.7 m/s at z_ref = 10m (standard reference height)
- α = 0.25 for urban terrain (Manhattan)
- Result: U(248m) ≈ 80 m/s at Citicorp roof
```

#### Implementation with groovyBC (swak4foam):
```cpp
// In 0/U
inlet
{
    type            codedFixedValue;
    value           uniform (44.7 0 0);

    code
    #{
        const vectorField& Cf = patch().Cf();  // Face centers
        vectorField& U = *this;

        scalar U_ref = 44.7;
        scalar z_ref = 10.0;
        scalar alpha = 0.25;
        scalar z_min = 0.5;  // Avoid log(0)

        forAll(Cf, i)
        {
            scalar z = max(Cf[i].z(), z_min);
            scalar U_mag = U_ref * pow(z / z_ref, alpha);
            U[i] = vector(U_mag, 0, 0);
        }
    #};
}
```

#### With Turbulence Profiles:
```cpp
// TKE: k(z) = (u*)^2 / sqrt(Cmu)
// where u* = friction velocity ≈ U_ref / 20 for urban

// Omega: omega(z) = u* / (kappa × z × sqrt(Cmu))
// where kappa = 0.41 (von Karman constant)
```

**Impact:** 40-80% higher velocities at roof → 2-3× higher loads (U²)

**Effort:** 2-4 hours coding + validation

---

### 4.2 Ground Roughness (Priority: MEDIUM)

**Current Issue:** Slip ground (frictionless)
**Target:** Rough wall with z₀ = 0.5-1.0m (urban terrain)

```cpp
// In 0/U
ground
{
    type            atmBoundaryLayerInletVelocity;
    Uref            44.7;
    Href            248;        // Building height
    z0              uniform 1.0;  // Roughness length (urban)
    zGround         uniform 0.0;
}

// In 0/k
ground
{
    type            kqRWallFunction;
    value           uniform 1;
}

// In 0/omega
ground
{
    type            omegaWallFunction;
    Cmu             0.09;
    kappa           0.41;
    E               9.8;
    beta1           0.075;
    value           uniform 1;
}
```

**Effect:** Reduced velocities near ground (wake region slower), more realistic wake structure

---

## 5. Fluid-Structure Interaction (FSI)

### 5.1 One-Way Coupled (Current → Enhanced)

**Current:** CFD pressure → MATLAB (manual copy)
**Target:** Automated one-way coupling

#### Export Pressure Distribution:
```cpp
// In system/controlDict functions
surfacePressure
{
    type            surfaceFieldValue;
    libs            ("libfieldFunctionObjects.so");
    writeControl    writeTime;
    writeFields     true;
    surfaceFormat   vtk;
    regionType      patch;
    name            citicorp_tower;
    operation       none;
    fields          (p);
}
```

#### Parse in MATLAB:
```matlab
% Read VTK pressure from OpenFOAM
vtkFile = 'postProcessing/surfacePressure/3000/citicorp_tower.vtk';
fid = fopen(vtkFile);
% Parse POINT_DATA section...
% Map pressures to structural nodes (nearest-neighbor or RBF interpolation)
% Apply to FEA model in analyze_chevron_fea.m
```

**Effort:** 4-8 hours for automation script

---

### 5.2 Two-Way Coupled FSI (Priority: LOW — research scope)

**Target:** Flow ↔ structure bidirectional coupling

#### Option A: preCICE (Recommended)
**Architecture:** OpenFOAM + MATLAB/CalculiX + preCICE coupler

```
┌────────────┐  Pressure   ┌──────────┐  Displacement  ┌────────────┐
│ OpenFOAM   │────────────▶│ preCICE  │────────────────▶│  MATLAB or │
│ (fluid)    │◀────────────│ (coupler)│◀────────────────│ CalculiX   │
└────────────┘  Mesh motion └──────────┘  Forces        └────────────┘
```

**Steps:**
1. Install preCICE: `sudo apt install libprecice-dev`
2. Install OpenFOAM adapter: `git clone https://github.com/precice/openfoam-adapter`
3. Create preCICE config (`precice-config.xml`)
4. Wrap MATLAB structural solver as preCICE participant
5. Run: `mpirun -np 12 foamRun -solver incompressibleFluid : -np 1 matlab -r structural_solver`

**Challenges:**
- **Mesh mapping:** CFD surface (triangles) ↔ structural nodes (sparse)
- **Time step coordination:** CFD (Δt=0.005s) vs structural (Δt=0.001s or adaptive)
- **Stability:** Added mass effect (light structure in dense fluid) requires implicit coupling
- **Debugging:** Two solvers + coupler = 3× complexity

**Effort:** 2-4 weeks initial setup, 1-2 weeks per case

**Payoff:**
- Predict actual building deflection under wind load
- Capture aeroelastic effects (deflection changes pressure → changes load)
- Model TMD dynamics in structural solver

---

### 5.3 Simplified TMD Model (Priority: MEDIUM)

**Target:** Add 400-ton TMD to MATLAB structural model without full FSI

#### Approach in MATLAB:
```matlab
% In analyze_chevron_fea.m or new module_13_tmd_dynamics.m

% TMD parameters (from literature on Citicorp TMD)
m_TMD = 400e3;      % 400 tons = 400,000 kg
k_TMD = ?;          % Tuned to building frequency (calculate)
c_TMD = ?;          % Damping ratio ζ = 0.05-0.15

% Building natural frequency (from eigenvalue analysis)
omega_n = 2*pi * 0.15;  % 0.15 Hz (Morgenstern 1995)

% TMD tuning (optimal frequency ratio μ ≈ 0.95-1.0)
omega_TMD = 0.97 * omega_n;
k_TMD = m_TMD * omega_TMD^2;

zeta_TMD = 0.10;  % 10% critical damping
c_TMD = 2 * zeta_TMD * sqrt(k_TMD * m_TMD);

% Two-DOF system: building + TMD
M = [m_building, 0; 0, m_TMD];
K = [k_building + k_TMD, -k_TMD; -k_TMD, k_TMD];
C = [c_building + c_TMD, -c_TMD; -c_TMD, c_TMD];

% Solve for displacement under harmonic wind load
F_wind_t = @(t) F0 * sin(omega_shed * t);  % Vortex shedding force

[T, X] = ode45(@(t,x) [x(3:4); M \ (F_wind_t(t) - K*x(1:2) - C*x(3:4))], ...
               [0 100], [0 0 0 0]);

% Plot displacement with/without TMD
```

**Effort:** 1-2 days for implementation, assuming structural eigenvalues known

---

## 6. Computational Infrastructure

### 6.1 GPU Acceleration with ESI v2512 (Priority: HIGH for speed)

**Current:** Foundation v13 CPU (10 cores, 4 hours)
**Target:** ESI v2512 GPU (RTX 4070, 1-2 hours estimated)

#### Setup Steps:
```bash
# 1. Activate ESI environment in WSL
source /opt/OpenFOAM/OpenFOAM-v2512/etc/bashrc

# 2. Copy case to new directory
cp -r ~/citicorp_cfd ~/citicorp_cfd_gpu
cd ~/citicorp_cfd_gpu

# 3. Modify system/fvSolution for PETSc/AmgX
```

```cpp
solvers
{
    p
    {
        solver          petsc;
        petsc
        {
            options
            {
                ksp_type    cg;
                pc_type     amgx;
                pc_amgx_config_file "$FOAM_CASE/system/amgx_pcg.json";
            }
        }
        tolerance       1e-6;
        relTol          0.01;
    }

    U
    {
        solver          PBiCGStab;  // Keep on CPU (less benefit from GPU)
        preconditioner  DILU;
        tolerance       1e-6;
        relTol          0.1;
    }
}
```

**4. Create AmgX config:**
```json
// system/amgx_pcg.json
{
    "config_version": 2,
    "solver": {
        "solver": "AMG",
        "presweeps": 1,
        "postsweeps": 3,
        "max_iters": 100,
        "convergence": "RELATIVE_INI_CORE",
        "tolerance": 1e-6,
        "norm": "L2",
        "cycle": "V",
        "smoother": "BLOCK_JACOBI",
        "coarsest_sweeps": 2
    }
}
```

**5. Load PETSc library in controlDict:**
```cpp
libs
(
    "libpetscFoam.so"
);
```

**6. Run:**
```bash
mpirun -np 4 simpleFoam -parallel  # 4 CPU processes share 1 GPU
```

**Speedup Estimate:**
- Pressure solve: 50-70% of runtime → 5-10× faster on GPU
- Overall: 2-3× faster (1.5-2 hours vs 4 hours)

**Caveats:**
- GPU benefit plateaus below 500k cells (overhead dominates)
- 3.16M cell mesh is good candidate (2-5× expected)
- Memory: 8GB VRAM adequate for this mesh

**Effort:** 1-2 hours setup (configs already exist), 2 hours test run

---

### 6.2 Hopper HPC Migration (Priority: HIGH for large runs)

**See HOPPER_SETUP.md for full details**

**Quick Summary:**
- 128-core AMD EPYC nodes (256 GB RAM)
- GPU nodes: NVIDIA A100 (40 GB)
- Slurm batch system
- Target: 8-30M cell meshes, parametric studies (10-20 wind directions)

**Effort:** 2-4 hours initial setup, then 1 hour per new case

---

## 7. Validation & Verification

### 7.1 Wind Tunnel Data Comparison (Priority: HIGH for publication)

**Challenge:** No published wind tunnel data for Citicorp specifically

**Option A: Generic Square Cylinder Benchmarks**
- **Source:** AIJ Benchmark (Architectural Institute of Japan)
- **Geometry:** 1:1:4 aspect ratio prism (similar to Citicorp)
- **Data:** Pressure coefficients (Cp) at 5 face locations × 16 wind angles
- **Access:** https://www.aij.or.jp/jpn/publish/cfdguide/

**Validation Metrics:**
```python
# Compare Cp at building faces
Cp_CFD = (p - p_ref) / (0.5 * rho * U_ref^2)
Cp_experiment = [-1.2, -0.8, -0.5, 0.7, -0.6]  # From AIJ

error = np.abs(Cp_CFD - Cp_experiment) / np.abs(Cp_experiment)
print(f"Mean error: {error.mean()*100:.1f}%")
```

**Option B: Similar Building (Proxy)**
- Target tall buildings with published data (e.g., CAARC standard tall building)
- Run CFD on their geometry first to validate methodology
- Apply same methodology to Citicorp with confidence

**Effort:** 1 week literature search + 1 week benchmark runs

---

### 7.2 Turbulence Model Sensitivity (Priority: MEDIUM)

**Run same mesh with different models:**
1. **k-ω SST** (current, good for separation)
2. **k-ε Realizable** (better for free shear, worse at walls)
3. **Spalart-Allmaras** (one-equation, fast but less accurate)
4. **v²-f** (better near walls, expensive)

**Comparison:**
| Model | Cd (tower+stilts) | Runtime (3000 iter) | Separation Point |
|-------|-------------------|---------------------|------------------|
| k-ω SST | 3.24 | 4 hrs | z/H = 0.3 |
| k-ε RNG | ? | 3.5 hrs | ? |
| S-A | ? | 2.5 hrs | ? |

**Acceptance:** If all models agree within 10%, result is robust

**Effort:** 3 runs × 4 hours = 12 hours

---

## Summary: Effort vs Impact Matrix

| Improvement | Effort | Impact | Priority | Dependencies |
|-------------|--------|--------|----------|--------------|
| **Add slanted roof** | 2-4 hrs | High (geometry accuracy) | HIGH | None |
| **Rotate 45°** | 1 hr | High (historical accuracy) | HIGH | Slanted roof |
| **Atmospheric BL inlet** | 2-4 hrs | High (load magnitude) | HIGH | None |
| **GPU acceleration** | 1-2 hrs | Medium (2× speedup) | HIGH | ESI installed |
| **Mesh independence** | 2-7 days | High (validation) | HIGH | Hopper HPC |
| **Boundary layer mesh** | 1-2 days | Medium (separation accuracy) | MEDIUM | Stable solver |
| **LOD2 geometry** | 1-3 days | Medium (wake accuracy) | MEDIUM | CityJSON tools |
| **LES turbulence** | 1-4 weeks | Very High (unsteady loads) | MED-HIGH | Hopper HPC, fine mesh |
| **Two-way FSI** | 2-4 weeks | High (dynamic response) | LOW | preCICE, advanced |
| **Wind tunnel validation** | 2 weeks | Very High (credibility) | HIGH | Benchmark data |

---

## Recommended Sequence (Next 2 Weeks)

### Week 1: Quick Wins
**Day 1-2:**
- ✓ Add slanted roof to tower (4 hrs)
- ✓ Rotate building 45° (1 hr)
- ✓ Test new geometry (2 hrs mesh + 4 hrs solve)

**Day 3-4:**
- ✓ Implement atmospheric BL inlet (4 hrs)
- ✓ Ground roughness BC (2 hrs)
- ✓ Rerun with realistic BCs (4 hrs solve)
- ✓ Compare Cd changes (should increase 1.5-2×)

**Day 5:**
- ✓ Setup GPU run on ESI v2512 (2 hrs)
- ✓ Benchmark CPU vs GPU (4 hrs total)

### Week 2: Validation & Scale-Up
**Day 6-7:**
- ✓ Setup Hopper HPC account (see HOPPER_SETUP.md)
- ✓ Port case to Hopper (4 hrs)
- ✓ Run mesh independence study: 3M, 8M, 15M cells (48 hrs walltime)

**Day 8-10:**
- ✓ Setup LES case (1 day)
- ✓ Submit LES job on Hopper (3-7 days walltime, run in background)
- ✓ While LES runs: process RANS results, write report

---

**Next Steps:** See [HOPPER_SETUP.md](HOPPER_SETUP.md) for HPC configuration details.
