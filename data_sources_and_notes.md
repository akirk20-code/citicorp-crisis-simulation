# Citicorp Center Simulation — Data Sources & Engineering Notes

## Primary References

1. **NIST (2021)** — "Modern Reassessment of the Citicorp Building Design Wind Loads"
   - PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC8240663/
   - Key finding: Face winds govern over quartering winds (contradicts 1978 analysis)
   - Overturning moment data, ASCE 7 methodology validation

2. **NIST Blog** — "Blown Away: Revisiting a Famous Engineering Case"
   - https://www.nist.gov/blogs/taking-measure/blown-away-revisiting-famous-engineering-case

3. **Morgenstern (1995)** — "The Fifty-Nine-Story Crisis", The New Yorker
   - PDF: https://people.duke.edu/~hpgavin/cee421/citicorp1.pdf
   - LeMessurier's own account, engineering decision-making narrative

4. **Wikipedia** — Citicorp Center Engineering Crisis
   - https://en.wikipedia.org/wiki/Citicorp_Center_engineering_crisis
   - Verified building dimensions, timeline, repair details

5. **PBS Building Big** — Citicorp Center Databank
   - https://www.pbs.org/wgbh/buildingbig/wonder/structure/citicorp.html

6. **Online Ethics Center** — William LeMessurier Case Study
   - https://onlineethics.org/cases/engineers-and-scientists-behaving-well/william-lemessurier-fifty-nine-story-crisis-lesson
   - Ethics/risk framing relevant to OR 750

## Engineering Standards Used

- **ASCE 7-22**: Minimum Design Loads and Associated Criteria for Buildings and Other Structures
  - Ch. 26: Wind Loads — General Requirements
  - Ch. 27: Wind Loads — MWFRS (Directional Procedure)
  - Table 26.10-1: Velocity Pressure Exposure Coefficients (Kz)
  - Table 26.6-1: Wind Directionality Factor (Kd = 0.85)
  - Sec. 26.11.5: Gust Effect Factor for Flexible Buildings (Gf)

- **AISC 360-22**: Specification for Structural Steel Buildings
  - Chapter J: Design of Connections
  - Table J3.2: Nominal Strength of Fasteners (A325-N: Fnv=54, Fnt=90 ksi)
  - Eq. J3-3a: Combined Tension-Shear Interaction

## Building Data

| Parameter | Value | Source |
|-----------|-------|--------|
| Height | 915 ft (279 m) | Wikipedia, PBS |
| Plan dimensions | 157 × 157 ft | Wikipedia |
| Stories | 59 | Wikipedia |
| Stilt height | 114 ft (9 stories) | Wikipedia |
| Stilt cross-section | 24 × 24 ft | Wikipedia |
| Corner cantilever | 72 ft | Wikipedia |
| Total steel | 24,000 short tons | Wikipedia |
| Fundamental period | 6.5 sec | NIST paper |
| Structural damping | 1% of critical | Typical for steel |
| TMD mass | 400 short tons | Wikipedia, PBS |
| TMD damping effect | +8% (total ~9%) | Historical record |
| Chevron tiers | 5 tiers of 8 stories each | Wikipedia |
| Bolted joints | ~200 total | Wikipedia |

## Wind Data

| Parameter | Value | Source |
|-----------|-------|--------|
| 1970 NYC design wind | 100 mph (fastest-mile) | Historical NYC code |
| ASCE 7-22 Risk Cat II | 115 mph (3-sec gust) | ASCE 7-22 Fig. 26.5-1B |
| ASCE 7-22 Risk Cat III | 120 mph (3-sec gust) | ASCE 7-22 Fig. 26.5-1C |
| Quartering fail speed (no TMD) | ~70 mph | LeMessurier analysis |
| Quartering fail speed (with TMD) | ~85 mph | LeMessurier analysis |
| Return period (no TMD) | 16 years | LeMessurier / Wikipedia |
| Return period (with TMD) | 55 years | LeMessurier / Wikipedia |
| Exposure category | B (urban) | ASCE 7-22 for Manhattan |

## Tool Options Evaluated

| Tool | Type | Verdict |
|------|------|---------|
| **MATLAB R2024b** | FEA + Monte Carlo + Visualization | **Selected** — no toolboxes, full control |
| **OpenFOAM v2312+** | CFD (steady RANS, k-omega SST) | **Implemented** — validates ASCE 7 analytical loads via 3D CFD |
| CalculiX | FEA solver | Good but less visual for presentations |
| Code_Aster | Industrial FEA | Overkill for educational demo |
| FEniCS | Python FEA | Less visual than MATLAB |
| SAP2000/ETABS | Commercial FEA | Not available |

## OpenFOAM CFD Setup

| Parameter | Value |
|-----------|-------|
| Solver | simpleFoam (steady-state, incompressible RANS) |
| Turbulence model | k-omega SST |
| Domain | 720 m × 720 m × 600 m (wind in +x) |
| Base mesh | 36 × 36 × 30 cells (~20 m base resolution) |
| Final mesh target | ~2–3 M cells (snappyHexMesh) |
| Tower surface refinement | Level 3–4 (3–6 m cells) |
| Surrounding bldg refinement | Level 2–3 (6–12 m cells) |
| ABL inlet | Richards & Hoxey (1993): Uref = 44.7 m/s, Zref = 10 m, z0 = 0.1 m |
| Surrounding buildings | 18 approximate boxes (Midtown East, heights 70–199 m) |
| Parallel | 20 cores, scotch decomposition |
| Post-processing | Forces on tower/stilts, y+, field min/max, residuals |

ABL reference: Richards, P.J. and Hoxey, R.P. (1993). "Appropriate boundary conditions
for computational wind engineering models using the k-epsilon turbulence model."
*J. Wind Eng. Ind. Aerodyn.*, 46–47, 145–153.

## Key Simulation Decisions

1. **3D FEA with direct stiffness method** — no toolbox dependencies, full transparency
2. **ASCE 7-22 wind analysis** — current code, not simplified power-law
3. **AISC 360-22 bolt interaction** — Eq. J3-3a for combined tension-shear
4. **Gumbel distribution** — standard for annual extreme wind speeds
5. **100,000 Monte Carlo samples** — sufficient convergence for Pf ~ 0.06
6. **Representative W-sections** — W14x730 columns, W14x176 braces (actual sections unknown)

## NIST Reassessment Note

The 2021 NIST paper concludes that LeMessurier's 1978 analysis was conservative. Using modern Database-Assisted Design (DAD) methodology with actual wind tunnel data, face winds govern the design — not quartering winds. The original error was treating quartering wind loads as the static sum of two perpendicular components, which overestimates the actual (correlated, dynamic) load combination.

However, the bolted connections were genuinely under-designed relative to the original welded specification, regardless of wind angle. The repair was appropriate.
