function [fig, fea_results] = analyze_chevron_fea(params, wind_results)
%ANALYZE_CHEVRON_FEA  3D Frame FEA of Citicorp Center lateral system
%  Direct stiffness method with 3D beam-column elements (12 DOF/element).
%  Models chevron bracing on all 4 faces with proper node deduplication.
%  Compares perpendicular (0 deg) vs quartering (45 deg) wind loading.

    W  = params.width;           % 157 ft
    Hs = params.stilt_height;    % 114 ft
    sh = params.story_height;    % ~13.6 ft
    E  = params.E_steel * 144;   % ksi -> ksf (for ft units)
    G  = params.G_steel * 144;   % ksi -> ksf
    n_tiers = params.n_tiers;    % 6 chevron tiers
    tier_h  = params.tier_stories * sh;  % tier height in ft
    H       = params.height;             % 915 ft total height

    % Transfer truss (W-pattern between stilt top and V-bracing zone)
    if isfield(params, 'transfer_stories') && params.transfer_stories > 0
        n_transfer = params.transfer_stories;
    else
        n_transfer = 0;
    end
    brace_base = Hs + n_transfer * sh;  % V-tiers start above transfer truss

    % Section properties (in^2 -> ft^2, in^4 -> ft^4)
    A_col   = params.A_column / 144;
    I_col   = params.I_column / 20736;
    A_brace = params.A_brace  / 144;
    I_brace = params.I_brace  / 20736;
    A_beam  = params.A_beam   / 144;
    I_beam  = params.I_beam   / 20736;

    %% ===== BUILD NODES (with deduplication) =====
    % Helper: add node only if it doesn't already exist
    nodes = zeros(0, 3);

    add_node = @(xyz) add_unique_node(xyz);
    node_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');

    function id = add_unique_node(xyz)
        key = sprintf('%.1f_%.1f_%.1f', xyz(1), xyz(2), xyz(3));
        if node_map.isKey(key)
            id = node_map(key);
        else
            nodes(end+1, :) = xyz;
            id = size(nodes, 1);
            node_map(key) = id;
        end
    end

    % Face definitions: [left_corner, midpoint, right_corner]
    face_pts = {[0,0], [W/2,0], [W,0];      % South (y=0)
                [W,0], [W,W/2], [W,W];      % East  (x=W)
                [W,W], [W/2,W], [0,W];      % North (y=W)
                [0,W], [0,W/2], [0,0]};     % West  (x=0)

    % Transfer truss: /\/\/\/\/\ (5 V-peaks, 10 diagonals per face)
    % 6 bottom nodes (fractions: 0, 0.250, 0.400, 0.600, 0.750, 1.000)
    % 7 top nodes (fractions: 0, 0.200, 0.350, 0.500, 0.650, 0.800, 1.000)
    bot_frac_xfer = [0, 0.250, 0.400, 0.600, 0.750, 1.000];
    top_frac_xfer = [0, 0.200, 0.350, 0.500, 0.650, 0.800, 1.000];
    N_bot_xfer = 6;  % bottom nodes per face
    N_top_xfer = 7;  % top nodes per face (including corners)

    % Create nodes at tier boundaries (V-bracing zone starts at brace_base)
    tier_nodes = zeros(4, n_tiers+1, 3);  % face, level, position -> node_id

    for f = 1:4
        for t = 0:n_tiers
            z_lev = brace_base + t * tier_h;
            for p = 1:3
                xy = face_pts{f, p};
                tier_nodes(f, t+1, p) = add_unique_node([xy(1), xy(2), z_lev]);
            end
        end
    end

    % --- Transfer zone nodes (separate bottom and top positions) ---
    % Bottom: 6 nodes per face at z=Hs
    % Top: 7 nodes per face at z=brace_base (including corners)
    xfer_bot = [];
    xfer_top = [];
    stilt_mids = zeros(4, 1);  % stilt center nodes at z=Hs (face midpoints)
    if n_transfer > 0
        xfer_bot = zeros(4, N_bot_xfer);  % face x bottom_position -> node_id
        xfer_top = zeros(4, N_top_xfer);  % face x top_position -> node_id
        for f = 1:4
            for k = 1:N_bot_xfer
                frac = bot_frac_xfer(k);
                xy = face_pts{f,1}*(1-frac) + face_pts{f,3}*frac;
                xfer_bot(f, k) = add_unique_node([xy, Hs]);
            end
            for k = 1:N_top_xfer
                frac = top_frac_xfer(k);
                xy = face_pts{f,1}*(1-frac) + face_pts{f,3}*frac;
                xfer_top(f, k) = add_unique_node([xy, brace_base]);
            end
            % Stilt center at z=Hs (face midpoint, frac 0.5)
            xy_mid = face_pts{f, 2};
            stilt_mids(f) = add_unique_node([xy_mid, Hs]);
        end
    end

    n_nodes = size(nodes, 1);
    ndof = n_nodes * 6;
    fprintf('    Nodes: %d, DOF: %d\n', n_nodes, ndof);

    %% ===== BUILD ELEMENTS (with deduplication) =====
    % Elements: [n1, n2, A, I, type]
    %   type: 1=column, 2=brace, 3=beam
    elements = zeros(0, 5);
    elem_set = containers.Map('KeyType', 'char', 'ValueType', 'logical');

    function add_elem(n1, n2, A_e, I_e, etype)
        key1 = sprintf('%d_%d', min(n1,n2), max(n1,n2));
        if ~elem_set.isKey(key1) && n1 ~= n2
            elements(end+1, :) = [n1, n2, A_e, I_e, etype];
            elem_set(key1) = true;
        end
    end

    % --- Columns (vertical members on each face) ---
    for f = 1:4
        for p = 1:3
            for t = 1:n_tiers
                add_elem(tier_nodes(f,t,p), tier_nodes(f,t+1,p), A_col, I_col, 1);
            end
        end
    end

    % --- Chevron braces on each face ---
    % 1st tier: truncated V (\_/) at transfer truss inner W-peaks (0.350, 0.650)
    %           2 brace members + 1 horizontal beam
    % Remaining tiers: normal sharp V with apex at bottom center (0.500)
    %           2 brace members, no horizontal beam
    for f = 1:4
        for t = 1:n_tiers
            tl = tier_nodes(f, t+1, 1);    % top-left corner
            tr = tier_nodes(f, t+1, 3);    % top-right corner
            z_bot_t = brace_base + (t-1) * tier_h;
            if t == 1
                % Truncated V: legs at W-peak fractions 0.350 and 0.650
                xy_L = face_pts{f,1}*(1-0.350) + face_pts{f,3}*0.350;
                xy_R = face_pts{f,1}*(1-0.650) + face_pts{f,3}*0.650;
                bm_L = add_unique_node([xy_L, z_bot_t]);
                bm_R = add_unique_node([xy_R, z_bot_t]);
                add_elem(tl, bm_L, A_brace, I_brace, 2);
                add_elem(tr, bm_R, A_brace, I_brace, 2);
                add_elem(bm_L, bm_R, A_beam, I_beam, 3);
            else
                % Normal V: apex at face midpoint (fraction 0.500)
                xy_mid = face_pts{f,1}*0.5 + face_pts{f,3}*0.5;
                apex = add_unique_node([xy_mid, z_bot_t]);
                add_elem(tl, apex, A_brace, I_brace, 2);
                add_elem(tr, apex, A_brace, I_brace, 2);
            end
        end
    end

    % --- Transfer truss: /\/\/\/\/\ (stilt top to V-brace base) ---
    % 10 diagonals per face (5 V-peaks) with bottom/top nodes at different positions
    if n_transfer > 0
        A_xfer = A_brace * 1.5;   % heavier transfer truss members
        I_xfer = I_brace * 2.0;

        for f = 1:4
            % 10 diagonals: zigzag /\/\/\/\/\
            % B=bottom(Hs), T=top(brace_base)
            add_elem(xfer_bot(f,1), xfer_top(f,2), A_xfer, I_xfer, 2);  % diag 1:  /
            add_elem(xfer_top(f,2), xfer_bot(f,2), A_xfer, I_xfer, 2);  % diag 2:  \
            add_elem(xfer_bot(f,2), xfer_top(f,3), A_xfer, I_xfer, 2);  % diag 3:  /
            add_elem(xfer_top(f,3), xfer_bot(f,3), A_xfer, I_xfer, 2);  % diag 4:  \
            add_elem(xfer_bot(f,3), xfer_top(f,4), A_xfer, I_xfer, 2);  % diag 5:  /
            add_elem(xfer_top(f,4), xfer_bot(f,4), A_xfer, I_xfer, 2);  % diag 6:  \
            add_elem(xfer_bot(f,4), xfer_top(f,5), A_xfer, I_xfer, 2);  % diag 7:  /
            add_elem(xfer_top(f,5), xfer_bot(f,5), A_xfer, I_xfer, 2);  % diag 8:  \
            add_elem(xfer_bot(f,5), xfer_top(f,6), A_xfer, I_xfer, 2);  % diag 9:  /
            add_elem(xfer_top(f,6), xfer_bot(f,6), A_xfer, I_xfer, 2);  % diag 10: \

            % Corner verticals (bot corners <-> top corners)
            add_elem(xfer_bot(f,1), xfer_top(f,1), A_col, I_col, 1);
            add_elem(xfer_bot(f,6), xfer_top(f,7), A_col, I_col, 1);

            % Bottom chord (6 nodes, 5 segments)
            for k = 1:(N_bot_xfer-1)
                add_elem(xfer_bot(f,k), xfer_bot(f,k+1), A_beam, I_beam, 3);
            end

            % Top chord (7 nodes, 6 segments)
            for k = 1:(N_top_xfer-1)
                add_elem(xfer_top(f,k), xfer_top(f,k+1), A_beam, I_beam, 3);
            end
        end
        fprintf('    Transfer truss: %d-story, 10-diagonal pattern added\n', n_transfer);
    end

    % --- Horizontal beams at each tier level ---
    for f = 1:4
        for t = 1:n_tiers+1
            add_elem(tier_nodes(f,t,1), tier_nodes(f,t,2), A_beam, I_beam, 3);
            add_elem(tier_nodes(f,t,2), tier_nodes(f,t,3), A_beam, I_beam, 3);
        end
    end

    % --- Floor diaphragm: cross-building transfer girders ---
    % The real Citicorp building has rigid floor diaphragms connecting all
    % four faces. Without these, wind load on one face can't transfer to
    % the bracing on perpendicular faces. The transfer girders run between
    % the 4 stilt locations (face midpoints), which are the primary load
    % paths in this building.
    A_diaphragm = A_beam * 4;    % very stiff transfer girder
    I_diaphragm = I_beam * 8;    % rigid diaphragm behavior

    for t = 1:n_tiers+1
        % Major transfer girders: South<->North, East<->West
        sm = tier_nodes(1, t, 2);  % South face midpoint
        nm = tier_nodes(3, t, 2);  % North face midpoint
        em = tier_nodes(2, t, 2);  % East face midpoint
        wm = tier_nodes(4, t, 2);  % West face midpoint
        add_elem(sm, nm, A_diaphragm, I_diaphragm, 3);
        add_elem(em, wm, A_diaphragm, I_diaphragm, 3);

        % Diagonal bracing in plan (torsional resistance)
        add_elem(sm, em, A_beam * 2, I_beam * 4, 3);
        add_elem(em, nm, A_beam * 2, I_beam * 4, 3);
        add_elem(nm, wm, A_beam * 2, I_beam * 4, 3);
        add_elem(wm, sm, A_beam * 2, I_beam * 4, 3);
    end

    % --- Transfer zone diaphragm (at z=Hs, connecting face midpoints) ---
    if n_transfer > 0
        xf_sm = stilt_mids(1);  % South midpoint at z=Hs (stilt center, frac 0.5)
        xf_em = stilt_mids(2);  % East midpoint
        xf_nm = stilt_mids(3);  % North midpoint
        xf_wm = stilt_mids(4);  % West midpoint
        add_elem(xf_sm, xf_nm, A_diaphragm, I_diaphragm, 3);
        add_elem(xf_em, xf_wm, A_diaphragm, I_diaphragm, 3);
        add_elem(xf_sm, xf_em, A_beam * 2, I_beam * 4, 3);
        add_elem(xf_em, xf_nm, A_beam * 2, I_beam * 4, 3);
        add_elem(xf_nm, xf_wm, A_beam * 2, I_beam * 4, 3);
        add_elem(xf_wm, xf_sm, A_beam * 2, I_beam * 4, 3);
    end

    n_elem = size(elements, 1);
    fprintf('    Elements: %d (after deduplication)\n', n_elem);

    %% ===== ASSEMBLE GLOBAL STIFFNESS MATRIX =====
    K = sparse(ndof, ndof);

    for e = 1:n_elem
        n1 = elements(e, 1);
        n2 = elements(e, 2);
        Ae = elements(e, 3);
        Ie = elements(e, 4);

        x1 = nodes(n1,:);
        x2 = nodes(n2,:);
        L = norm(x2 - x1);

        if L < 1e-3, continue; end

        % 3D beam element (12x12)
        J_torsion = 2 * Ie;  % approximate torsion constant
        ke_local = beam3d_stiffness(E, Ae, Ie, Ie, G, J_torsion, L);
        T = rotation_matrix_3d(x1, x2);

        ke_global = T' * ke_local * T;

        dofs = [6*(n1-1)+(1:6), 6*(n2-1)+(1:6)];
        K(dofs, dofs) = K(dofs, dofs) + ke_global;
    end

    %% ===== APPLY WIND LOADS =====
    F_perp = zeros(ndof, 1);
    F_quar = zeros(ndof, 1);

    z_wind = wind_results.z;
    p_perp = wind_results.p_net_perp;
    p_q1   = wind_results.p_net_face1_45;
    p_q2   = wind_results.p_net_face2_45;

    brace_top = brace_base + n_tiers * tier_h;

    for n = 1:n_nodes
        z_n = nodes(n, 3);
        if z_n <= brace_base + 1, continue; end  % skip base/transfer nodes

        % Tributary height: half-tier at top boundary only, full tier elsewhere
        if abs(z_n - brace_top) < 1
            trib_h = tier_h / 2;  % top of bracing zone: half tributary above
        else
            trib_h = tier_h;      % interior nodes: full tier tributary
        end

        p_0    = interp1(z_wind, p_perp, z_n, 'linear', 'extrap');
        p_45_1 = interp1(z_wind, p_q1, z_n, 'linear', 'extrap');
        p_45_2 = interp1(z_wind, p_q2, z_n, 'linear', 'extrap');

        % Tributary width: each face has 3 nodes spanning W, so ~W/3 per node
        % But corner nodes are shared between faces — use W/4 for corners, W/2 for mid
        x_n = nodes(n, 1);
        y_n = nodes(n, 2);

        % Determine which face(s) this node belongs to and apply loads
        % South face (y ≈ 0): wind pushes in +Y direction (perpendicular)
        if abs(y_n) < 0.1
            trib_w = W / 4;
            if abs(x_n - W/2) < 1, trib_w = W / 2; end  % midpoint gets more
            F_node = p_0 * trib_h * trib_w;  % force in lbs
            F_perp(6*(n-1)+2) = F_perp(6*(n-1)+2) + F_node / 1000;  % convert to kips
        end

        % Quartering wind: south face (+Y) and west face (+X)
        if abs(y_n) < 0.1  % south face
            trib_w = W / 4;
            if abs(x_n - W/2) < 1, trib_w = W / 2; end
            F_node = p_45_1 * trib_h * trib_w;
            F_quar(6*(n-1)+2) = F_quar(6*(n-1)+2) + F_node / 1000;
        end
        if abs(x_n) < 0.1  % west face
            trib_w = W / 4;
            if abs(y_n - W/2) < 1, trib_w = W / 2; end
            F_node = p_45_2 * trib_h * trib_w;
            F_quar(6*(n-1)+1) = F_quar(6*(n-1)+1) + F_node / 1000;
        end
    end

    % Total applied loads
    F_perp_total = sum(abs(F_perp));
    F_quar_total = sum(abs(F_quar));
    fprintf('    Total applied load (perp.):  %.0f kips\n', F_perp_total);
    fprintf('    Total applied load (quart.): %.0f kips\n', F_quar_total);

    %% ===== BOUNDARY CONDITIONS =====
    fixed_dofs = [];
    for n = 1:n_nodes
        if abs(nodes(n,3) - Hs) < 1.0
            fixed_dofs = [fixed_dofs, 6*(n-1)+(1:6)];
        end
    end
    free_dofs = setdiff(1:ndof, fixed_dofs);

    fprintf('    Fixed DOF: %d, Free DOF: %d\n', length(fixed_dofs), length(free_dofs));

    %% ===== SOLVE =====
    K_ff = K(free_dofs, free_dofs);

    % Check conditioning
    K_diag = full(diag(K_ff));
    min_diag = min(abs(K_diag));
    max_diag = max(abs(K_diag));
    diag_ratio = max_diag / max(min_diag, 1e-20);
    fprintf('    K diagonal range: [%.2e, %.2e], ratio: %.2e\n',...
            min_diag, max_diag, diag_ratio);

    % Guard: check for zero diagonal entries (mechanism/singularity)
    n_zero_diag = sum(abs(K_diag) < 1e-10 * max_diag);
    if n_zero_diag > 0
        warning('FEA:SingularSystem', ...
            '%d near-zero diagonal entries detected — model may be a mechanism.', n_zero_diag);
    end

    % Guard: check condition number (estimated via normest for large sparse)
    if length(free_dofs) < 2000
        cond_est = condest(K_ff);
    else
        cond_est = diag_ratio;  % fallback for very large systems
    end
    fprintf('    Estimated condition number: %.2e\n', cond_est);
    if cond_est > 1e14
        warning('FEA:IllConditioned', ...
            'Stiffness matrix is ill-conditioned (cond ~%.1e). Results may be unreliable.', cond_est);
    end

    % Solve
    U_perp = zeros(ndof, 1);
    U_quar = zeros(ndof, 1);

    U_perp(free_dofs) = K_ff \ F_perp(free_dofs);
    U_quar(free_dofs) = K_ff \ F_quar(free_dofs);

    % Guard: check solution for NaN/Inf
    if any(~isfinite(U_perp)) || any(~isfinite(U_quar))
        error('FEA:BadSolution', ...
            'Solver produced NaN/Inf displacements — check boundary conditions and model connectivity.');
    end

    max_u_perp = max(abs(U_perp)) * 12;  % ft -> inches
    max_u_quar = max(abs(U_quar)) * 12;

    fprintf('    Max displacement (perp.):  %.2f in (%.4f ft)\n',...
            max_u_perp, max(abs(U_perp)));
    fprintf('    Max displacement (quart.): %.2f in (%.4f ft)\n',...
            max_u_quar, max(abs(U_quar)));

    %% ===== VALIDATION: STATIC EQUILIBRIUM CHECK =====
    % The sum of applied forces must equal the sum of support reactions.
    % Reactions are extracted from K*U at the fixed DOFs.

    % Y-direction check (perpendicular wind loads in +Y)
    F_app_y = sum(F_perp(2:6:ndof));       % total applied Y-force (kips)
    R_all   = K * U_perp;                  % full force vector: applied at free, reactions at fixed
    fixed_y_mask = false(ndof, 1);
    for fd = fixed_dofs
        if mod(fd, 6) == 2  % Y-direction DOF (1-indexed: 2,8,14,...)
            fixed_y_mask(fd) = true;
        end
    end
    R_y_support = sum(R_all(fixed_y_mask));  % total Y-reaction at supports

    % Equilibrium: F_applied + R_reactions = 0  →  error = |F + R| / |F|
    eq_error_perp = abs(F_app_y + R_y_support) / max(abs(F_app_y), 1e-10);

    fprintf('\n    --- VALIDATION CHECKS ---\n');
    fprintf('    Equilibrium (perp Y): applied=%.1f kips, reactions=%.1f kips, error=%.2e\n',...
            F_app_y, -R_y_support, eq_error_perp);

    % Quartering: check both X and Y
    F_app_x_q = sum(F_quar(1:6:ndof));
    F_app_y_q = sum(F_quar(2:6:ndof));
    R_all_q   = K * U_quar;
    fixed_x_mask = false(ndof, 1);
    for fd = fixed_dofs
        if mod(fd, 6) == 1  % X-direction DOF
            fixed_x_mask(fd) = true;
        end
    end
    R_x_q = sum(R_all_q(fixed_x_mask));
    R_y_q = sum(R_all_q(fixed_y_mask));
    eq_error_qx = abs(F_app_x_q + R_x_q) / max(abs(F_app_x_q), 1e-10);
    eq_error_qy = abs(F_app_y_q + R_y_q) / max(abs(F_app_y_q), 1e-10);

    fprintf('    Equilibrium (quar X): applied=%.1f kips, reactions=%.1f kips, error=%.2e\n',...
            F_app_x_q, -R_x_q, eq_error_qx);
    fprintf('    Equilibrium (quar Y): applied=%.1f kips, reactions=%.1f kips, error=%.2e\n',...
            F_app_y_q, -R_y_q, eq_error_qy);

    if max([eq_error_perp, eq_error_qx, eq_error_qy]) < 1e-6
        fprintf('    PASS: Static equilibrium satisfied (error < 1e-6)\n');
    else
        warning('FEA:Equilibrium', 'Equilibrium error exceeds 1e-6 — check model.');
    end

    % Simplified cantilever comparison (order-of-magnitude check)
    % Treat building above stilts as a cantilever: L = H - Hs
    % With effective EI from column stiffness on two faces (2 corners per face)
    L_cant = (H - Hs);                      % ft — cantilever height
    EI_col = E * I_col;                     % ksf * ft^4
    n_cols_resisting = 4;                   % corners on 2 loaded faces
    EI_total = n_cols_resisting * EI_col;
    % Uniformly distributed load: w = V_base / L  (approximate)
    w_cant = abs(F_app_y) / L_cant;         % kips/ft
    delta_cant = w_cant * L_cant^4 / (8 * EI_total) * 12;  % inches

    fprintf('    Cantilever estimate: delta = %.1f in (FEA = %.1f in, ratio = %.2f)\n',...
            delta_cant, max_u_perp, max_u_perp / max(delta_cant, 0.01));
    fprintf('    (Ratio > 1: FEA is more flexible due to brace/diaphragm participation)\n');

    % Drift ratio
    drift_ratio = H / (max_u_perp / 12);
    fprintf('    Drift ratio: H/%.0f (typical limit: H/400 to H/200)\n', drift_ratio);

    %% ===== COMPUTE MEMBER FORCES =====
    forces_perp = zeros(n_elem, 1);
    forces_quar = zeros(n_elem, 1);

    for e = 1:n_elem
        n1 = elements(e, 1);
        n2 = elements(e, 2);
        Ae = elements(e, 3);
        Ie = elements(e, 4);

        x1 = nodes(n1,:);
        x2 = nodes(n2,:);
        L = norm(x2 - x1);
        if L < 1e-3, continue; end

        J_t = 2 * Ie;
        ke_local = beam3d_stiffness(E, Ae, Ie, Ie, G, J_t, L);
        T = rotation_matrix_3d(x1, x2);

        dofs = [6*(n1-1)+(1:6), 6*(n2-1)+(1:6)];

        u_loc_p = T * U_perp(dofs);
        u_loc_q = T * U_quar(dofs);

        f_loc_p = ke_local * u_loc_p;
        f_loc_q = ke_local * u_loc_q;

        % Axial force at node j (DOF 7 in local)
        forces_perp(e) = f_loc_p(7);
        forces_quar(e) = f_loc_q(7);
    end

    %% ===== PACK RESULTS =====
    fea_results = struct();
    fea_results.nodes = nodes;
    fea_results.elements = elements;
    fea_results.U_perp = U_perp;
    fea_results.U_quar = U_quar;
    fea_results.forces_perp = forces_perp;
    fea_results.forces_quar = forces_quar;
    fea_results.max_disp = max(max_u_perp, max_u_quar);
    fea_results.n_nodes = n_nodes;
    fea_results.n_elem = n_elem;

    brace_idx = find(elements(:,5) == 2);
    fea_results.brace_forces_perp = forces_perp(brace_idx);
    fea_results.brace_forces_quar = forces_quar(brace_idx);
    fea_results.brace_idx = brace_idx;
    fea_results.eq_error_perp = eq_error_perp;
    fea_results.eq_error_quar = max(eq_error_qx, eq_error_qy);
    fea_results.drift_ratio = drift_ratio;
    fea_results.delta_cantilever = delta_cant;

    %% ===== PLOTTING =====
    fig = figure('Name','3D FEA — Chevron Bracing System',...
                 'Position',[50 30 1500 900],'Color','k');

    scale_perp = 0.02 * W / max(max(abs(U_perp)), 1e-10);  % auto-scale
    scale_quar = 0.02 * W / max(max(abs(U_quar)), 1e-10);

    % --- Left: Deformed shape (perpendicular wind) ---
    ax1 = subplot(1,3,1); hold on;
    set(ax1,'Color','k','XColor','w','YColor','w','ZColor','w');
    plot_deformed(nodes, elements, U_perp, forces_perp, scale_perp);
    title(sprintf('Perpendicular Wind (0°)\nMax disp: %.1f in', max_u_perp),...
          'Color','w','FontSize',12);
    view([-37 25]); axis equal; grid on;
    set(ax1,'GridColor',[0.3 0.3 0.3]);
    xlabel('X (ft)','Color','w'); ylabel('Y (ft)','Color','w'); zlabel('Z (ft)','Color','w');

    % --- Center: Deformed shape (quartering wind) ---
    ax2 = subplot(1,3,2); hold on;
    set(ax2,'Color','k','XColor','w','YColor','w','ZColor','w');
    plot_deformed(nodes, elements, U_quar, forces_quar, scale_quar);
    title(sprintf('Quartering Wind (45°)\nMax disp: %.1f in', max_u_quar),...
          'Color','w','FontSize',12);
    view([-37 25]); axis equal; grid on;
    set(ax2,'GridColor',[0.3 0.3 0.3]);
    xlabel('X (ft)','Color','w'); ylabel('Y (ft)','Color','w'); zlabel('Z (ft)','Color','w');

    % --- Right: Force comparison in brace elements ---
    ax3 = subplot(1,3,3); hold on;
    set(ax3,'Color','k','XColor','w','YColor','w');

    n_braces = length(brace_idx);
    f_p = fea_results.brace_forces_perp;
    f_q = fea_results.brace_forces_quar;

    % Sort by quartering force magnitude
    [~, sort_idx] = sort(abs(f_q), 'descend');
    f_p_sorted = f_p(sort_idx);
    f_q_sorted = f_q(sort_idx);

    barh(1:n_braces, [f_p_sorted, f_q_sorted], 'grouped');
    colormap(ax3, [0.2 0.6 1; 1 0.3 0.3]);

    legend({'0° (Perp.)','45° (Quart.)'},...
           'TextColor','w','Color',[0.15 0.15 0.15]);
    xlabel('Axial Force (kips)','Color','w');
    ylabel('Brace Element #','Color','w');
    title('Chevron Brace Axial Forces','Color','w','FontSize',12);
    grid on; set(ax3,'GridColor',[0.3 0.3 0.3]);

    % Highlight critical braces
    ratio = abs(f_q_sorted) ./ max(abs(f_p_sorted), 0.001);
    critical = find(ratio > 1.3);
    for c = critical'
        plot(f_q_sorted(c), c, 'o', 'Color', [1 1 0], 'MarkerSize', 8, 'LineWidth', 2);
    end

    sgtitle('Module 3: 3D Frame FEA — Citicorp Center Lateral System',...
            'Color','w','FontSize',16,'FontWeight','bold');
end

%% ===== HELPER FUNCTIONS =====

function ke = beam3d_stiffness(E, A, Iz, Iy, G, J, L)
%BEAM3D_STIFFNESS  12x12 stiffness matrix for 3D Euler-Bernoulli beam
%  DOF order per node: [ux, uy, uz, rx, ry, rz]

    L2 = L^2;
    L3 = L^3;

    ke = zeros(12, 12);

    % Axial (DOFs 1, 7)
    ea = E*A/L;
    ke(1,1) = ea;   ke(1,7) = -ea;
    ke(7,1) = -ea;  ke(7,7) = ea;

    % Torsion (DOFs 4, 10)
    gj = G*J/L;
    ke(4,4) = gj;   ke(4,10) = -gj;
    ke(10,4) = -gj; ke(10,10) = gj;

    % Bending in x-y plane (about z-axis) — DOFs 2,6,8,12
    c1 = 12*E*Iz/L3;  c2 = 6*E*Iz/L2;  c3 = 4*E*Iz/L;  c4 = 2*E*Iz/L;
    ke(2,2) = c1;    ke(2,6) = c2;    ke(2,8) = -c1;   ke(2,12) = c2;
    ke(6,2) = c2;    ke(6,6) = c3;    ke(6,8) = -c2;   ke(6,12) = c4;
    ke(8,2) = -c1;   ke(8,6) = -c2;   ke(8,8) = c1;    ke(8,12) = -c2;
    ke(12,2) = c2;   ke(12,6) = c4;   ke(12,8) = -c2;  ke(12,12) = c3;

    % Bending in x-z plane (about y-axis) — DOFs 3,5,9,11
    d1 = 12*E*Iy/L3;  d2 = 6*E*Iy/L2;  d3 = 4*E*Iy/L;  d4 = 2*E*Iy/L;
    ke(3,3) = d1;    ke(3,5) = -d2;   ke(3,9) = -d1;   ke(3,11) = -d2;
    ke(5,3) = -d2;   ke(5,5) = d3;    ke(5,9) = d2;    ke(5,11) = d4;
    ke(9,3) = -d1;   ke(9,5) = d2;    ke(9,9) = d1;    ke(9,11) = d2;
    ke(11,3) = -d2;  ke(11,5) = d4;   ke(11,9) = d2;   ke(11,11) = d3;
end

function T = rotation_matrix_3d(x1, x2)
%ROTATION_MATRIX_3D  12x12 transformation for 3D beam element

    dx = x2 - x1;
    L = norm(dx);
    lx = dx / L;

    % Local y-axis: use global Z as reference (or global X if nearly vertical)
    if abs(lx(3)) > 0.95
        ref = [1, 0, 0];
    else
        ref = [0, 0, 1];
    end

    lz = cross(lx, ref);
    lz = lz / norm(lz);
    ly = cross(lz, lx);

    R = [lx; ly; lz];

    T = zeros(12, 12);
    T(1:3, 1:3) = R;
    T(4:6, 4:6) = R;
    T(7:9, 7:9) = R;
    T(10:12, 10:12) = R;
end

function plot_deformed(nodes, elements, U, forces, scale)
%PLOT_DEFORMED  Plot deformed shape with color-coded forces

    n_elem = size(elements, 1);
    max_f = max(abs(forces));
    if max_f == 0, max_f = 1; end

    % Undeformed (ghost)
    for e = 1:n_elem
        n1 = elements(e, 1); n2 = elements(e, 2);
        plot3([nodes(n1,1) nodes(n2,1)], [nodes(n1,2) nodes(n2,2)],...
              [nodes(n1,3) nodes(n2,3)], ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.5);
    end

    % Deformed (colored by force)
    for e = 1:n_elem
        n1 = elements(e, 1); n2 = elements(e, 2);
        etype = elements(e, 5);

        x1d = nodes(n1,:) + scale * U(6*(n1-1)+(1:3))';
        x2d = nodes(n2,:) + scale * U(6*(n2-1)+(1:3))';

        f_norm = forces(e) / max_f;
        if f_norm > 0
            col = [min(1, 0.3 + 0.7*abs(f_norm)), 0.15, 0.15];
        else
            col = [0.15, 0.15, min(1, 0.3 + 0.7*abs(f_norm))];
        end

        lw = 1 + etype * 0.5;  % columns thicker
        plot3([x1d(1) x2d(1)], [x1d(2) x2d(2)], [x1d(3) x2d(3)],...
              '-', 'Color', col, 'LineWidth', lw);
    end
end
