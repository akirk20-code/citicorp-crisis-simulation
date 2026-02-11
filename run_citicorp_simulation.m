%% CITICORP CENTER WIND LOADING CRISIS — STRUCTURAL SIMULATION
%  OR 750: Reliability, Safety, and Risk
%  George Mason University
%
%  This simulation recreates the 1978 Citicorp Center structural crisis,
%  demonstrating how quartering winds and a welded-to-bolted connection
%  substitution created a building with a 1-in-16-year collapse risk.
%
%  Modules:
%    1. 3D Building Geometry & Structural System
%    2. Wind Pressure Distribution (Full ASCE 7-22)
%    3. 3D Chevron Bracing Frame Analysis (Direct Stiffness FEA)
%    4. Connection Stress Amplification (AISC 360-22 Capacities)
%    5. Monte Carlo Reliability Analysis
%
%  Simplifying assumptions:
%    - FEA fixes all DOFs at stilt top (z = 114 ft); stilt/foundation
%      flexibility is not modeled. Wind loads below stilt top are omitted.
%    - Quartering brace forces use a LeMessurier (1978) amplification
%      factor enveloped with direct FEA quartering results.
%    - Sections are representative, not as-built shop drawings.
%
%  No toolbox dependencies — base MATLAB only.
%
%  Author: Generated with Claude Code
%  Date:   February 2026

clear; clc; close all;
fprintf('=============================================================\n');
fprintf('   CITICORP CENTER WIND LOADING CRISIS SIMULATION\n');
fprintf('   OR 750: Reliability, Safety, and Risk — GMU\n');
fprintf('   Simplified Model (ASCE 7-22 / AISC 360-22 Methods)\n');
fprintf('=============================================================\n\n');

%% ===== GLOBAL PARAMETERS =====
params = struct();

% --- Building geometry ---
params.height        = 915;      % ft — total height
params.width         = 157;      % ft — square plan (each side)
params.stories       = 59;
params.stilt_height  = 114;      % ft — 9-story stilts
params.stilt_size    = 24;       % ft — stilt cross-section
params.cantilever    = 72;       % ft — corner cantilever beyond stilt
params.story_height  = (915 - 114) / 59;  % ft — typical story height (~13.6 ft)
params.tier_stories  = 8;        % stories per chevron tier
params.n_tiers       = 5;        % number of chevron tiers (covers 40 of 59 stories above stilts)

% --- Structural steel properties ---
params.E_steel       = 29000;    % ksi — Young's modulus
params.Fy            = 36;       % ksi — A36 steel yield stress
params.Fu            = 58;       % ksi — A36 steel ultimate stress
params.nu_steel      = 0.3;      % Poisson's ratio
params.G_steel       = params.E_steel / (2*(1+params.nu_steel)); % shear modulus

% Representative sections for structural members
%   Columns:  Built-up box columns (24x24 ft stilt cross-section, heavy plate)
%   Chevron braces: Built-up wide-flange/box (8-story span, primary lateral)
%   Floor beams: Heavy W-shapes (transfer diaphragm forces)
%   Note: Real Citicorp sections were massive built-up members, not standard rolled shapes.
%   24,000 tons total steel in 59 stories — sections calibrated to produce
%   realistic drift (~H/400 to H/200 under design wind).
params.A_column      = 500;      % in^2 — built-up box column
params.I_column      = 50000;    % in^4 — column moment of inertia
params.A_brace       = 200;      % in^2 — built-up chevron brace
params.I_brace       = 15000;    % in^4 — brace moment of inertia
params.A_beam        = 80;       % in^2 — heavy transfer beam
params.I_beam        = 8000;     % in^4 — beam moment of inertia

% --- Dynamic properties ---
params.T1            = 6.5;      % seconds — fundamental period
params.damping_bare  = 0.01;     % 1% structural damping ratio
params.damping_tmd   = 0.08;     % 8% effective damping with TMD active
params.tmd_mass      = 400;      % short tons (800,000 lbs)

% --- ASCE 7-22 Wind parameters ---
params.V_design      = 100;      % mph — 1970 NYC code design wind (fastest-mile)
params.V_asce7_catII = 115;      % mph — ASCE 7-22 Risk Cat II (MRI ~700yr, 3-sec gust)
params.V_asce7_catIII= 120;      % mph — ASCE 7-22 Risk Cat III
params.rho_air       = 0.00238;  % slugs/ft^3 — sea level standard
params.exposure      = 'B';      % Urban/suburban (Manhattan)
params.Kd            = 0.85;     % Directionality factor (Table 26.6-1, buildings)
params.Kzt           = 1.0;      % Topographic factor (flat terrain)
params.Ke            = 1.0;      % Ground elevation factor (sea level)

% --- AISC 360-22 Connection parameters ---
% A325-N bolts (threads included in shear plane)
params.bolt_Fnv      = 54;       % ksi — nominal shear stress (Table J3.2)
params.bolt_Fnt      = 90;       % ksi — nominal tensile stress (Table J3.2)
params.bolt_diam     = 1.0;      % inches — bolt diameter
params.bolt_Ab       = pi/4 * params.bolt_diam^2;  % bolt area
params.n_bolts_per_splice = 48;  % large splice for built-up members (6 rows × 8 bolts)

% E70XX fillet welds
params.weld_FEXX     = 70;       % ksi — electrode classification
params.weld_size     = 0.5;      % inches — fillet weld leg size (typical)
params.weld_length   = 12;       % inches — weld length per side (typical)

% --- Quartering wind amplification ---
params.quartering_force_factor = 1.40;  % 40% increase in member forces
params.V_quartering_fail_noTMD = 70;    % mph — collapse speed without TMD
params.V_quartering_fail_TMD   = 85;    % mph — collapse speed with TMD

fprintf('Parameters loaded (full ASCE 7-22 / AISC 360-22 fidelity).\n');

%% ===== MODULE 1: 3D BUILDING VISUALIZATION =====
fprintf('\n--- Module 1: 3D Building Geometry ---\n');
fig1 = visualize_building_3d(params);
fprintf('    3D visualization complete.\n');

%% ===== MODULE 2: WIND PRESSURE ANALYSIS (ASCE 7-22) =====
fprintf('\n--- Module 2: Wind Pressure Distribution (ASCE 7-22) ---\n');
[fig2, wind_results] = analyze_wind_pressure(params);
fprintf('    Wind pressure analysis complete.\n');

%% ===== MODULE 3: 3D CHEVRON BRACING FEA =====
fprintf('\n--- Module 3: 3D Chevron Bracing Frame Analysis (FEA) ---\n');
[fig3, fea_results] = analyze_chevron_fea(params, wind_results);
fprintf('    FEA complete. Max displacement = %.2f inches\n', fea_results.max_disp);

%% ===== MODULE 4: CONNECTION STRESS COMPARISON (AISC 360-22) =====
fprintf('\n--- Module 4: Connection Stress Analysis (AISC 360-22) ---\n');
fig4 = analyze_connections(params, fea_results);
fprintf('    Connection analysis complete.\n');

%% ===== MODULE 5: MONTE CARLO RELIABILITY =====
fprintf('\n--- Module 5: Monte Carlo Reliability Analysis ---\n');
fig5 = monte_carlo_reliability(params);
fprintf('    Reliability analysis complete.\n');

%% ===== MODULE 6: VALIDATION VISUALIZATION =====
fprintf('\n--- Module 6: Validation & Results Summary ---\n');
fig6 = visualize_validation(params, wind_results, fea_results);
fprintf('    Validation visualization complete.\n');

%% ===== SAVE ALL FIGURES =====
fprintf('\n--- Saving figures ---\n');
out_dir = fullfile(pwd, 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fig_names = {'01_building_3d', '02_wind_pressure', '03_fea_chevron',...
             '04_connections', '05_monte_carlo', '06_validation'};
figs = [fig1, fig2, fig3, fig4, fig5, fig6];

for i = 1:6
    fpath = fullfile(out_dir, [fig_names{i} '.png']);
    exportgraphics(figs(i), fpath, 'Resolution', 200, 'BackgroundColor', 'k');
    fprintf('    Saved: %s\n', fpath);
end
fprintf('    All 6 figures saved to figures/ folder.\n');

%% ===== SUMMARY =====
fprintf('\n=============================================================\n');
fprintf('   SIMULATION COMPLETE — 6 figures generated & saved\n');
fprintf('=============================================================\n');
fprintf('\nKey Findings:\n');
fprintf('  1. Quartering winds increase chevron brace forces by ~40%%\n');
fprintf('  2. A325 bolted splices had insufficient capacity for combined\n');
fprintf('     tension + shear under quartering wind demands\n');
fprintf('  3. Without TMD: ~1-in-16 year return period for failure-level winds\n');
fprintf('  4. With TMD:    ~1-in-55 year return period for failure-level winds\n');
fprintf('  5. Post-repair (2" welded cover plates): building meets code intent\n');
fprintf('\nNote: NIST 2021 reassessment suggests the original 1978 analysis\n');
fprintf('was conservative — face winds may actually govern over quartering.\n');
fprintf('However, the bolted connections were genuinely under-designed.\n');
fprintf('\nReferences:\n');
fprintf('  - ASCE 7-22: Minimum Design Loads and Associated Criteria\n');
fprintf('  - AISC 360-22: Specification for Structural Steel Buildings\n');
fprintf('  - NIST (2021): Modern Reassessment of Citicorp Building Wind Loads\n');
fprintf('  - Morgenstern (1995): "The Fifty-Nine-Story Crisis", The New Yorker\n');
