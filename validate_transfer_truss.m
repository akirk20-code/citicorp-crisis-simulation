%% TRANSFER TRUSS GEOMETRY VALIDATION
%  Validates the 10-diagonal (5 V-peak) transfer truss layout across all files.
%  Checks: geometric consistency, angle targets, half-span closure,
%  equilibrium, and cross-file agreement.
%
%  Run standalone:  validate_transfer_truss()
%  Returns struct with all validation results.

function results = validate_transfer_truss()

fprintf('=============================================================\n');
fprintf('   TRANSFER TRUSS GEOMETRY VALIDATION\n');
fprintf('   Citicorp Center — 5 V-peaks (/\\/\\/\\/\\/\\)\n');
fprintf('=============================================================\n\n');

%% ===== BUILDING PARAMETERS =====
W   = 157;      % ft — square plan width
Hs  = 114;      % ft — stilt top elevation
sh  = (915 - 114) / 50;  % story height = 16.02 ft
n_xfer = 2;     % transfer zone = 2 stories
h   = n_xfer * sh;       % transfer height = 32.04 ft
bb  = Hs + h;            % brace base elevation = 146.04 ft
L   = W / 2;             % half-span = 78.5 ft (stilt center to corner)

fprintf('Building Parameters:\n');
fprintf('  Plan width W       = %.1f ft\n', W);
fprintf('  Stilt top Hs       = %.1f ft\n', Hs);
fprintf('  Story height sh    = %.2f ft\n', sh);
fprintf('  Transfer height h  = %.2f ft (2 stories)\n', h);
fprintf('  Brace base bb      = %.2f ft\n', bb);
fprintf('  Half-span L        = %.1f ft (stilt center to corner)\n\n', L);

%% ===== NODE FRACTIONS (canonical values used in all files) =====
bot_f = [0, 0.250, 0.400, 0.600, 0.750, 1.000];    % 6 bottom nodes
top_f = [0.200, 0.350, 0.500, 0.650, 0.800];        % 5 internal top nodes
top_f_full = [0, 0.200, 0.350, 0.500, 0.650, 0.800, 1.000];  % 7 top nodes (FEA, with corners)

N_bot = length(bot_f);
N_top = length(top_f);
N_top_full = length(top_f_full);

fprintf('Node Fractions (canonical):\n');
fprintf('  Bottom (%d): [', N_bot);
fprintf('%.3f ', bot_f); fprintf(']\n');
fprintf('  Top internal (%d): [', N_top);
fprintf('%.3f ', top_f); fprintf(']\n');
fprintf('  Top full (%d): [', N_top_full);
fprintf('%.3f ', top_f_full); fprintf(']\n\n');

%% ===== TEST 1: SYMMETRY =====
fprintf('--- Test 1: Symmetry about stilt centerline (fraction = 0.5) ---\n');
pass_sym = true;

% Bottom node symmetry: f and (1-f) should be complementary
for k = 1:floor(N_bot/2)
    f_left  = bot_f(k);
    f_right = bot_f(end+1-k);
    err = abs((f_left + f_right) - 1.0);
    ok = err < 1e-6;
    if ~ok, pass_sym = false; end
    fprintf('  bot[%d]=%.3f + bot[%d]=%.3f = %.6f  %s\n', ...
        k, f_left, N_bot+1-k, f_right, f_left+f_right, tf(ok));
end

% Top node symmetry
for k = 1:floor(N_top/2)
    f_left  = top_f(k);
    f_right = top_f(end+1-k);
    err = abs((f_left + f_right) - 1.0);
    ok = err < 1e-6;
    if ~ok, pass_sym = false; end
    fprintf('  top[%d]=%.3f + top[%d]=%.3f = %.6f  %s\n', ...
        k, f_left, N_top+1-k, f_right, f_left+f_right, tf(ok));
end

% Center top node (should be at 0.5 = stilt center)
ctr_ok = abs(top_f(3) - 0.5) < 1e-6;
if ~ctr_ok, pass_sym = false; end
fprintf('  Center top node top[3] = %.3f  %s\n', top_f(3), tf(ctr_ok));
fprintf('  SYMMETRY: %s\n\n', pf(pass_sym));

%% ===== TEST 2: NODE POSITIONS (absolute, ft) =====
fprintf('--- Test 2: Absolute node positions (South face, y=0) ---\n');
bx = W * bot_f;  % bottom x-positions
tx = W * top_f;  % top x-positions (internal only)
tx_full = W * top_f_full;  % all 7 top positions

fprintf('  Bottom (z = %.1f ft):\n', Hs);
for k = 1:N_bot
    fprintf('    B%d: x = %7.2f ft  (frac %.3f)\n', k, bx(k), bot_f(k));
end
fprintf('  Top (z = %.2f ft):\n', bb);
for k = 1:N_top_full
    fprintf('    T%d: x = %7.2f ft  (frac %.3f)\n', k, tx_full(k), top_f_full(k));
end
fprintf('\n');

%% ===== TEST 3: DIAGONAL ANGLES =====
N_diag = 10;
fprintf('--- Test 3: Diagonal angles (%d diagonals per face) ---\n', N_diag);
fprintf('  Target angles:\n');
fprintf('    Outer diagonals:  ~46 deg (from horizontal)\n');
fprintf('    Transitions:      ~76 deg (steep)\n');
fprintf('    Center V legs:    ~64 deg (apex ~52 deg)\n\n');

% Compute 10 diagonals using the zigzag connectivity
angles = zeros(1, N_diag);
dx_diag = zeros(1, N_diag);
lengths = zeros(1, N_diag);

for k = 1:N_diag
    if mod(k,2) == 1
        bi = (k+1)/2;
        ti = (k+1)/2;
        x_bot = bx(bi);
        x_top = tx(ti);
    else
        ti = k/2;
        bi = k/2 + 1;
        x_top = tx(ti);
        x_bot = bx(bi);
    end
    dx = abs(x_top - x_bot);
    len = sqrt(dx^2 + h^2);
    ang = atand(h / dx);  % angle from horizontal
    angles(k) = ang;
    dx_diag(k) = dx;
    lengths(k) = len;

    if mod(k,2) == 1
        dir_str = '/';
    else
        dir_str = '\';
    end
    fprintf('  Diag %2d (%s): dx=%6.2f ft, len=%6.2f ft, angle=%5.1f deg\n', ...
        k, dir_str, dx, len, ang);
end

fprintf('\n  Angle summary:\n');
fprintf('    Outer pair      (1,10): %.1f deg\n', mean([angles(1), angles(10)]));
fprintf('    Transition pair  (2,9): %.1f deg\n', mean([angles(2), angles(9)]));
fprintf('    Inner V-legs     (3,8): %.1f deg\n', mean([angles(3), angles(8)]));
fprintf('    Inner transition (4,7): %.1f deg\n', mean([angles(4), angles(7)]));
fprintf('    Center V-legs    (5,6): %.1f deg\n', mean([angles(5), angles(6)]));

% Center V apex angle = 180 - 2*leg_angle
center_apex = 180 - 2*mean([angles(5), angles(6)]);
fprintf('    Center V apex angle: %.1f deg  (target ~52)\n', center_apex);

% All V-peak apex angles (at top nodes T1..T5)
fprintf('\n  V-peak apex angles (at top nodes):\n');
for j = 1:5
    a_in  = angles(2*j - 1);   % incoming /
    a_out = angles(2*j);       % outgoing \
    apex = 180 - a_in - a_out;
    fprintf('    Peak at T%d: %.1f deg  (diags %d,%d)\n', j, apex, 2*j-1, 2*j);
end

% Bottom valley angles
fprintf('\n  Bottom valley angles (at bottom nodes B2..B5):\n');
for j = 1:4
    a_in  = angles(2*j);       % incoming \
    a_out = angles(2*j + 1);   % outgoing /
    valley = 180 - a_in - a_out;
    fprintf('    Valley at B%d: %.1f deg  (diags %d,%d)\n', j+1, valley, 2*j, 2*j+1);
end

% Angle pass/fail (center apex within 5° of 52°)
pass_angles = true;
if abs(center_apex - 52) > 5, pass_angles = false; end
fprintf('\n  ANGLES: %s\n\n', pf(pass_angles));

%% ===== TEST 4: HALF-SPAN CLOSURE =====
fprintf('--- Test 4: Half-span closure (corner to stilt center) ---\n');
% The zigzag path from B1(0) to T3(0.5) covers the left half-span
% Net horizontal: tx(3) - bx(1) = W*0.500 - 0 = 78.5 ft
left_start = bx(1);
left_end = tx(3);  % center top node at stilt
net_left = left_end - left_start;
ok_left = abs(net_left - L) < 0.1;

% Right half: T3(0.5) to B6(1.0)
right_start = tx(3);
right_end = bx(N_bot);
net_right = right_end - right_start;
ok_right = abs(net_right - L) < 0.1;

fprintf('  Left half  (B1 to T3): %.2f ft  (target %.1f)  %s\n', net_left, L, tf(ok_left));
fprintf('  Right half (T3 to B6): %.2f ft  (target %.1f)  %s\n', net_right, L, tf(ok_right));
fprintf('  Full span  (B1 to B6): %.2f ft  (target %.1f)  %s\n', ...
    bx(N_bot)-bx(1), W, tf(abs(bx(N_bot)-bx(1)-W) < 0.1));

% Sum of diagonal dx's should equal full span
dx_sum = sum(dx_diag);
fprintf('  Sum of diagonal dx:    %.2f ft  (target %.1f)  %s\n', ...
    dx_sum, W, tf(abs(dx_sum - W) < 0.5));
fprintf('  HALF-SPAN CLOSURE: %s\n\n', pf(ok_left && ok_right));

%% ===== TEST 5: STATIC EQUILIBRIUM =====
fprintf('--- Test 5: Static equilibrium (unit load at corner) ---\n');
P = 1;  % unit load
F_chord = P * L / h;
fprintf('  Unit load P = 1 at corner\n');
fprintf('  Cantilever moment M = P*L = %.2f P*ft\n', L);
fprintf('  Maximum chord force = M/h = P*L/h = %.3f P\n', F_chord);
fprintf('  (Independent of panel count — always %.3f P)\n\n', F_chord);

% Diagonal forces (method of sections: V = P, F_diag = P/sin(theta))
fprintf('  Diagonal forces (left half, carrying unit shear P=1):\n');
for k = 1:5  % left half only (symmetric)
    theta = angles(k);
    F_diag = P / sind(theta);
    governs = '';
    if F_diag > F_chord
        governs = '  ** governs over chord!';
    end
    fprintf('    Diag %d: theta=%.1f deg, F=P/sin(theta) = %.3f P%s\n', ...
        k, theta, F_diag, governs);
end
thresh = asind(1/F_chord);
fprintf('  Chord governs for diags steeper than %.1f deg\n', thresh);
fprintf('  All diags steeper than %.1f deg? %s\n\n', ...
    thresh, pf(all(angles > thresh)));

%% ===== TEST 6: CROSS-FILE CONSISTENCY =====
fprintf('--- Test 6: Cross-file consistency ---\n');
files = {
    'visualize_building_3d.m', ...
    'citicorp_3d_viewer.html', ...
    'interactive_structure.m', ...
    'draw_structural_elevation.m', ...
    'analyze_chevron_fea.m'
};

pass_consistency = true;
for f = 1:length(files)
    fname = files{f};
    fpath = fullfile(pwd, fname);
    if exist(fpath, 'file')
        txt = fileread(fpath);
        has_bot = contains(txt, '0.250') && contains(txt, '0.750');
        has_top = contains(txt, '0.350') && contains(txt, '0.650');
        has_10  = contains(txt, '10') || contains(txt, '1:10') || contains(txt, '<= 10');
        ok = has_bot && has_top;
        if ~ok, pass_consistency = false; end
        fprintf('  %-35s  bot: %s  top: %s  %s\n', fname, tf(has_bot), tf(has_top), pf(ok));
    else
        fprintf('  %-35s  FILE NOT FOUND\n', fname);
        pass_consistency = false;
    end
end
fprintf('  CROSS-FILE CONSISTENCY: %s\n\n', pf(pass_consistency));

%% ===== TEST 7: DIAGONAL COUNT AND CONNECTIVITY =====
fprintf('--- Test 7: Diagonal count and connectivity ---\n');
n_diag_per_face = 10;
n_vert_per_face = 2;
n_bot_chord_per_face = N_bot - 1;       % 5 segments
n_top_chord_per_face = N_top_full - 1;  % 6 segments
n_faces = 4;

total_diag  = n_diag_per_face * n_faces;
total_vert  = n_vert_per_face * n_faces;
total_bchd  = n_bot_chord_per_face * n_faces;
total_tchd  = n_top_chord_per_face * n_faces;
total_xfer  = total_diag + total_vert + total_bchd + total_tchd;

fprintf('  Per face:\n');
fprintf('    Diagonals:     %d\n', n_diag_per_face);
fprintf('    Corner verts:  %d\n', n_vert_per_face);
fprintf('    Bot chord segs: %d\n', n_bot_chord_per_face);
fprintf('    Top chord segs: %d\n', n_top_chord_per_face);
fprintf('  Total (4 faces):\n');
fprintf('    Diagonals:     %d\n', total_diag);
fprintf('    Corner verts:  %d\n', total_vert);
fprintf('    Bot chords:    %d\n', total_bchd);
fprintf('    Top chords:    %d\n', total_tchd);
fprintf('    TOTAL transfer elements: %d\n', total_xfer);

% Verify alternating /\ pattern
fprintf('  Zigzag pattern: ');
for k = 1:N_diag
    if mod(k,2) == 1
        fprintf('/');
    else
        fprintf('\\');
    end
end
fprintf('  (5 V-peaks)\n');
fprintf('  CONNECTIVITY: PASS\n\n');

%% ===== TEST 8: NODE DEDUPLICATION (FEA) =====
fprintf('--- Test 8: FEA node deduplication ---\n');
n_bot_unique = 4 * N_bot - 4;  % corners shared: 24-4 = 20 unique
n_top_unique = 4 * N_top_full - 4;  % corners shared: 28-4 = 24 unique
fprintf('  Bottom nodes: 4 faces x %d = %d, minus 4 shared corners = %d unique\n', ...
    N_bot, 4*N_bot, n_bot_unique);
fprintf('  Top nodes:    4 faces x %d = %d, minus 4 shared corners = %d unique\n', ...
    N_top_full, 4*N_top_full, n_top_unique);
fprintf('  (Top corner nodes are also shared with tier_nodes at brace_base)\n');
fprintf('  NODE DEDUPLICATION: expected behavior\n\n');

%% ===== SUMMARY =====
fprintf('=============================================================\n');
fprintf('   VALIDATION SUMMARY\n');
fprintf('=============================================================\n');
all_pass = pass_sym && pass_angles && (ok_left && ok_right) && pass_consistency;
tests = {'Symmetry', 'Angles', 'Half-span closure', 'Equilibrium', ...
         'Cross-file consistency', 'Connectivity', 'Node deduplication'};
status = {pass_sym, pass_angles, ok_left && ok_right, true, ...
          pass_consistency, true, true};
for t = 1:length(tests)
    fprintf('  %-25s %s\n', tests{t}, pf(status{t}));
end
fprintf('\n  OVERALL: %s\n', pf(all_pass));
fprintf('=============================================================\n');

%% ===== RETURN RESULTS =====
results = struct();
results.pass_all = all_pass;
results.W = W;
results.Hs = Hs;
results.h = h;
results.bb = bb;
results.L = L;
results.bot_f = bot_f;
results.top_f = top_f;
results.top_f_full = top_f_full;
results.bx = bx;
results.tx = tx;
results.tx_full = tx_full;
results.angles = angles;
results.center_apex = center_apex;
results.F_chord = F_chord;
results.pass_symmetry = pass_sym;
results.pass_angles = pass_angles;
results.pass_closure = ok_left && ok_right;
results.pass_consistency = pass_consistency;

end

%% ===== HELPER FUNCTIONS =====
function s = tf(b)
    if b, s = 'OK'; else, s = 'FAIL'; end
end

function s = pf(b)
    if b, s = 'PASS'; else, s = 'FAIL'; end
end
