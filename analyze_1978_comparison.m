function [fig, results_1978] = analyze_1978_comparison(params, wind_results, fea_results)
%ANALYZE_1978_COMPARISON  1978-era vs modern wind analysis — the Citicorp crisis
%  Recreates the simplifying assumptions used during the original crisis:
%    Scenario A: Original design — face wind only (bolts sized for this)
%    Scenario B: LeMessurier 1978 — quartering via static superposition
%    Scenario C: Modern ASCE 7-22 — proper quartering wind analysis
%    Scenario D: NIST 2021 DAD — face winds govern (Duthinh et al.)
%
%  Reference: Duthinh (NIST, 2021): "The simplifying assumption was made of
%  decomposing the effects of a corner wind into the sum of the effects on
%  two adjacent faces, as if these wind loads were static and perfectly
%  correlated, which they are not."

    H  = params.height;         % 915 ft
    W  = params.width;          % 157 ft
    Hs = params.stilt_height;   % 114 ft
    z  = wind_results.z;        % height array from Module 2

    %% ===== 1972 ANSI A58.1 WIND PROFILE =====
    % Original NYC Building Code used ANSI A58.1-1972:
    %   - Wind speed: 100 mph (fastest-mile at 30 ft in open terrain)
    %   - Exposure B (urban): alpha=4.5, zg=1200 ft
    %   - No directionality factor Kd (introduced later)
    %   - Simpler gust factor G (~1.0 to 1.1 for tall buildings)

    V_1978  = 100;              % mph — fastest-mile
    zg      = 1200;             % ft  — gradient height (Exposure B)
    alpha72 = 4.5;              % power-law exponent (ANSI A58.1-1972)
    alpha22 = 7.0;              % power-law exponent (ASCE 7-22)

    % 1972 velocity pressure exposure coefficient
    Kz_1972 = 2.58 * (max(z, 15) / zg).^(2 / alpha72);
    % Modern ASCE 7-22 coefficient (from Module 2)
    Kz_2022 = 2.01 * (max(z, 15) / zg).^(2 / alpha22);

    % Velocity pressure profiles
    q_1972 = 0.00256 * Kz_1972 * V_1978^2;                       % 1972 (no Kd)
    q_2022 = 0.00256 * Kz_2022 * params.Kd * V_1978^2;          % ASCE 7-22 w/ Kd

    % Gust factors
    G_1972  = 1.1;              % approximate for tall buildings in 1970s code
    Gf_2022 = wind_results.Gf;  % flexible building Gf from Module 2

    % Pressure coefficients (same both eras)
    Cp_ww = 0.8;   Cp_lw = -0.5;

    %% ===== SCENARIO A: ORIGINAL DESIGN (FACE WIND ONLY) =====
    % NYC Building Code at time of construction checked perpendicular wind.
    % Connections were designed (as welds) for these forces.
    % Later substituted with bolts — no recheck performed.

    p_face_72 = q_1972 * G_1972 * (Cp_ww + abs(Cp_lw));   % psf, net
    V_face_72 = trapz(z, p_face_72 * W) / 1000;            % kips
    M_face_72 = trapz(z, p_face_72 * W .* z) / 1000;       % kip-ft

    fprintf('    --- Scenario A: Original Design (1972 ANSI, face wind) ---\n');
    fprintf('      Base shear:       %.0f kips\n', V_face_72);
    fprintf('      Overturning M:    %.0f kip-ft\n', M_face_72);

    %% ===== SCENARIO B: LEMESSURIER 1978 QUARTERING (STATIC SUPERPOSITION) =====
    % LeMessurier's graduate student found that quartering winds were never
    % checked. LeMessurier recalculated using static superposition:
    %   - Full face pressure applied to BOTH adjacent faces simultaneously
    %   - Treated as perfectly correlated (they are NOT in reality)
    %   - Resultant: sqrt(V1^2 + V2^2) = sqrt(2) * V_face for square building
    %   - Most-loaded braces: ~1.40x face brace forces (load path through diaphragm)

    V_quar_72     = sqrt(2) * V_face_72;                    % kips (static superposition)
    M_quar_72     = sqrt(2) * M_face_72;                    % kip-ft
    amp_static_sup = sqrt(2);                                % 1.414 (theoretical)
    amp_brace_1978 = params.quartering_force_factor;         % 1.40 (LeMessurier's finding)

    fprintf('    --- Scenario B: LeMessurier 1978 (static superposition) ---\n');
    fprintf('      Quartering shear: %.0f kips (x%.3f of face)\n', V_quar_72, amp_static_sup);
    fprintf('      Brace amp factor: %.2f\n', amp_brace_1978);

    %% ===== SCENARIO C: MODERN ASCE 7-22 =====
    % Proper quartering analysis with reduced Cp for oblique incidence,
    % directionality factor Kd=0.85, and flexible building gust factor Gf.

    V_face_22 = wind_results.V_base_perp / 1000;            % kips
    V_quar_22 = wind_results.V_base_quar / 1000;            % kips
    M_face_22 = wind_results.M_perp / 1000;                 % kip-ft
    M_quar_22 = wind_results.M_quar / 1000;                 % kip-ft
    amp_asce7 = V_quar_22 / max(V_face_22, 1);

    fprintf('    --- Scenario C: ASCE 7-22 (modern, V=%d mph) ---\n', V_1978);
    fprintf('      Face shear:       %.0f kips\n', V_face_22);
    fprintf('      Quartering shear: %.0f kips (x%.2f of face)\n', V_quar_22, amp_asce7);

    %% ===== SCENARIO D: NIST 2021 DAD REASSESSMENT =====
    % Duthinh et al. used Database-Assisted Design (wind tunnel data) and
    % found that FACE winds govern the overturning moment — NOT quartering.
    % The 1978 static superposition was overly conservative.

    amp_nist = 0.85;  % approximate: face governs, quartering is less critical
    V_quar_nist = amp_nist * V_face_22;

    fprintf('    --- Scenario D: NIST 2021 (DAD, face governs) ---\n');
    fprintf('      Face governs: quartering ratio ~ %.2f (less than face)\n', amp_nist);

    %% ===== BRACE FORCE & CONNECTION DEMAND =====
    % Scale FEA brace forces to 1978 loading level
    max_brace_perp = max(abs(fea_results.brace_forces_perp));
    max_brace_quar = max(abs(fea_results.brace_forces_quar));

    % Scale factor: 1972 face loads vs current FEA loads (ASCE 7-22)
    scale_72 = V_face_72 / max(V_face_22, 1);

    brace_face_72  = max_brace_perp * scale_72;             % 1972 face demand
    brace_quar_72  = brace_face_72  * amp_brace_1978;       % 1978 quartering demand
    brace_quar_mod = max_brace_quar;                        % ASCE 7-22 quartering (from FEA)
    brace_nist     = max_brace_perp * amp_nist * scale_72;  % NIST (face governs)

    % Bolt capacity (A325-N, from params)
    Ab      = params.bolt_Ab;
    n_bolts = params.n_bolts_per_splice;
    Fnv     = params.bolt_Fnv;          % 54 ksi
    Fnt     = params.bolt_Fnt;          % 90 ksi
    phi     = 0.75;

    % Weld capacity (E70XX fillet, original design)
    Fw = 0.60 * params.weld_FEXX;       % 42 ksi (weld metal shear)
    te = params.weld_size / sqrt(2);     % effective throat (0.354 in)
    weld_cap_per_in = phi * Fw * te;     % kips/in
    weld_length = params.weld_length * 4; % total weld (both sides, both flanges)
    weld_capacity = weld_cap_per_in * weld_length;  % kips total

    % Bolt capacity in pure shear (simplified — worst case)
    bolt_capacity = phi * Fnv * Ab * n_bolts;        % kips

    % Demand/Capacity ratios
    DC_weld_face  = brace_face_72  / max(weld_capacity, 1);
    DC_bolt_face  = brace_face_72  / max(bolt_capacity, 1);
    DC_bolt_quar  = brace_quar_72  / max(bolt_capacity, 1);
    DC_bolt_mod   = brace_quar_mod / max(bolt_capacity, 1);

    fprintf('\n    --- Connection Demand/Capacity ---\n');
    fprintf('      Weld capacity:     %.0f kips\n', weld_capacity);
    fprintf('      Bolt capacity:     %.0f kips\n', bolt_capacity);
    fprintf('      Weld D/C (face):   %.3f %s\n', DC_weld_face, status_str(DC_weld_face));
    fprintf('      Bolt D/C (face):   %.3f %s\n', DC_bolt_face, status_str(DC_bolt_face));
    fprintf('      Bolt D/C (quar.):  %.3f %s\n', DC_bolt_quar, status_str(DC_bolt_quar));

    %% ===== PACK RESULTS =====
    results_1978 = struct();
    results_1978.V_face_72     = V_face_72;
    results_1978.V_quar_72     = V_quar_72;
    results_1978.V_face_22     = V_face_22;
    results_1978.V_quar_22     = V_quar_22;
    results_1978.amp_static    = amp_static_sup;
    results_1978.amp_brace     = amp_brace_1978;
    results_1978.amp_asce7     = amp_asce7;
    results_1978.amp_nist      = amp_nist;
    results_1978.DC_weld_face  = DC_weld_face;
    results_1978.DC_bolt_face  = DC_bolt_face;
    results_1978.DC_bolt_quar  = DC_bolt_quar;
    results_1978.DC_bolt_mod   = DC_bolt_mod;
    results_1978.Kz_1972       = Kz_1972;
    results_1978.Kz_2022       = Kz_2022;
    results_1978.q_1972        = q_1972;
    results_1978.q_2022        = q_2022;
    results_1978.weld_capacity = weld_capacity;
    results_1978.bolt_capacity = bolt_capacity;

    %% ===== VISUALIZATION (6 panels) =====
    fig = figure('Name','1978 vs Modern — Citicorp Wind Analysis Comparison',...
                 'Position',[30 30 1700 950],'Color','k');

    % --- Panel 1: Kz Profile Comparison ---
    ax1 = subplot(2,3,1); hold on;
    set(ax1,'Color','k','XColor','w','YColor','w');

    plot(Kz_1972, z, '-', 'Color', [1 0.5 0.2], 'LineWidth', 2);
    plot(Kz_2022, z, '-', 'Color', [0.3 0.7 1], 'LineWidth', 2);
    yline(Hs, ':', 'Color', [0.9 0.6 0.1], 'LineWidth', 1,...
          'Label', 'Stilt top', 'LabelColor', [0.9 0.6 0.1]);

    xlabel('K_z','Color','w'); ylabel('Height (ft)','Color','w');
    title('Exposure Coefficient Profile','Color','w','FontSize',11);
    legend({sprintf('1972 ANSI (\\alpha=%.1f)', alpha72),...
            sprintf('ASCE 7-22 (\\alpha=%.1f)', alpha22)},...
           'TextColor','w','Color',[0.15 0.15 0.15],'Location','southeast');
    grid on; set(ax1,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 2: Velocity Pressure Comparison ---
    ax2 = subplot(2,3,2); hold on;
    set(ax2,'Color','k','XColor','w','YColor','w');

    plot(q_1972, z, '-', 'Color', [1 0.5 0.2], 'LineWidth', 2);
    plot(q_2022, z, '-', 'Color', [0.3 0.7 1], 'LineWidth', 2);
    yline(Hs, ':', 'Color', [0.9 0.6 0.1], 'LineWidth', 1);

    xlabel('q_z (psf)','Color','w'); ylabel('Height (ft)','Color','w');
    title('Velocity Pressure (100 mph)','Color','w','FontSize',11);
    legend({'1972 ANSI (no K_d)', 'ASCE 7-22 (K_d=0.85)'},...
           'TextColor','w','Color',[0.15 0.15 0.15],'Location','southeast');
    grid on; set(ax2,'GridColor',[0.3 0.3 0.3]);
    text(max(q_1972)*0.5, H*0.25,...
         sprintf('1972 is %.0f%% higher\nat rooftop (no K_d)',...
         (q_1972(end)/q_2022(end) - 1)*100),...
         'Color', [1 0.7 0.3], 'FontSize', 9);

    % --- Panel 3: The Quartering Assumption (schematic) ---
    ax3 = subplot(2,3,3);
    set(ax3,'Color','k','XColor','k','YColor','k');
    axis off; hold on;

    % Draw building plan view (square)
    bx = [0.15 0.55 0.55 0.15 0.15];
    by = [0.55 0.55 0.95 0.95 0.55];
    plot(ax3, bx, by, '-w', 'LineWidth', 2);
    text(0.35, 0.75, 'PLAN', 'Color', [0.5 0.5 0.5], 'FontSize', 9,...
         'HorizontalAlignment','center','Units','data');

    % 1978 assumption: full arrows on two faces
    % South face arrows (full)
    for yy = [0.60 0.70 0.80 0.90]
        annotation('arrow', [0.72 0.72], [yy-0.02 yy-0.02]+[-0.05 0],...
                   'Color', [1 0.4 0.3], 'LineWidth', 1.5);
    end
    % West face arrows (full)
    for xx = [0.20 0.30 0.40 0.50]
        annotation('arrow', [xx-0.02 xx-0.02]+[-0.05 0], [0.92 0.92],...
                   'Color', [1 0.4 0.3], 'LineWidth', 1.5);
    end

    text(0.35, 0.48, '1978: Full p on BOTH faces', 'Color', [1 0.4 0.3],...
         'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center',...
         'Units', 'data');
    text(0.35, 0.40, '"static + perfectly correlated"', 'Color', [1 0.6 0.4],...
         'FontSize', 9, 'HorizontalAlignment', 'center', 'FontAngle', 'italic',...
         'Units', 'data');

    % Modern: reduced arrows
    text(0.35, 0.28, 'ASCE 7-22: 0.75 \times cos^2(45\circ) p', 'Color', [0.3 0.7 1],...
         'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center',...
         'Units', 'data');
    text(0.35, 0.20, 'Properly reduced for oblique incidence', 'Color', [0.5 0.7 1],...
         'FontSize', 9, 'HorizontalAlignment', 'center',...
         'Units', 'data');

    text(0.35, 0.08, 'NIST 2021: Face winds govern', 'Color', [0.3 1 0.3],...
         'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center',...
         'Units', 'data');
    text(0.35, 0.01, '(quartering < face, per wind tunnel data)', 'Color', [0.5 0.8 0.5],...
         'FontSize', 9, 'HorizontalAlignment', 'center',...
         'Units', 'data');

    title('Quartering Wind Assumption','Color','w','FontSize',11);

    % --- Panel 4: Base Shear Comparison (bar chart) ---
    ax4 = subplot(2,3,4); hold on;
    set(ax4,'Color','k','XColor','w','YColor','w');

    shear_data = [V_face_72, V_quar_72, V_face_22, V_quar_22, V_quar_nist];
    colors = [0.6 0.6 0.6;   % A: 1972 face (gray — original design)
              1.0 0.3 0.2;   % B: 1978 quartering (red — crisis)
              0.3 0.7 1.0;   % C: ASCE face (blue)
              0.5 0.5 1.0;   % D: ASCE quartering (light blue)
              0.3 0.9 0.3];  % E: NIST (green)
    labels = {'1972 Face', '1978 Quar.', 'ASCE Face', 'ASCE Quar.', 'NIST 2021'};

    b = bar(shear_data, 0.7, 'FaceColor', 'flat');
    b.CData = colors;
    set(ax4, 'XTick', 1:5, 'XTickLabel', labels, 'XTickLabelRotation', 30);
    ylabel('Base Shear (kips)','Color','w');
    title('Base Shear Comparison','Color','w','FontSize',11);

    % Annotate values and amplification factors
    for i = 1:5
        text(i, shear_data(i) + max(shear_data)*0.03,...
             sprintf('%.0f', shear_data(i)),...
             'Color','w','FontSize',9,'FontWeight','bold','HorizontalAlignment','center');
    end
    % Amplification arrow from face to quartering (1978)
    text(1.5, max(shear_data)*0.85, sprintf('\\times%.2f', amp_static_sup),...
         'Color', [1 0.4 0.3], 'FontSize', 14, 'FontWeight', 'bold',...
         'HorizontalAlignment', 'center');

    grid on; set(ax4,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 5: Connection Demand/Capacity ---
    ax5 = subplot(2,3,5); hold on;
    set(ax5,'Color','k','XColor','w','YColor','w');

    % Bar chart: demand vs capacity for each scenario
    scenarios = {'Weld\n(face)', 'Bolt\n(face)', 'Bolt\n(1978 quar.)', 'Bolt\n(ASCE quar.)'};
    demands   = [brace_face_72, brace_face_72, brace_quar_72, brace_quar_mod];
    caps      = [weld_capacity, bolt_capacity, bolt_capacity, bolt_capacity];
    dc_vals   = [DC_weld_face, DC_bolt_face, DC_bolt_quar, DC_bolt_mod];

    x = 1:4;
    bw = 0.35;
    % Capacity bars (green)
    bar(x - bw/2, caps, bw, 'FaceColor', [0.2 0.7 0.2], 'EdgeColor', 'none');
    % Demand bars (colored by D/C)
    for i = 1:4
        if dc_vals(i) > 1.0
            dcol = [1 0.3 0.2];  % red — exceeds
        elseif dc_vals(i) > 0.9
            dcol = [1 0.7 0.2];  % amber — marginal
        else
            dcol = [0.3 0.7 1];  % blue — OK
        end
        bar(x(i) + bw/2, demands(i), bw, 'FaceColor', dcol, 'EdgeColor', 'none');
    end

    set(ax5, 'XTick', 1:4, 'XTickLabel',...
        {'Weld (face)', 'Bolt (face)', 'Bolt (1978 Q)', 'Bolt (ASCE Q)'},...
        'XTickLabelRotation', 25);
    ylabel('Force (kips)','Color','w');
    title('Connection Demand vs Capacity','Color','w','FontSize',11);
    legend({'Capacity','Demand'},'TextColor','w','Color',[0.15 0.15 0.15],...
           'Location','northwest');

    % D/C labels
    for i = 1:4
        if dc_vals(i) > 1.0
            st = sprintf('D/C=%.2f\nOVERSTRESS', dc_vals(i));
            sc = [1 0.3 0.2];
        else
            st = sprintf('D/C=%.2f\nOK', dc_vals(i));
            sc = [0.3 1 0.3];
        end
        text(i, max(demands(i), caps(i)) + max(demands)*0.06, st,...
             'Color', sc, 'FontSize', 8, 'FontWeight', 'bold',...
             'HorizontalAlignment', 'center');
    end
    grid on; set(ax5,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 6: Key Findings ---
    ax6 = subplot(2,3,6);
    set(ax6,'Color','k','XColor','k','YColor','k');
    axis off;

    lines = {
        '1978 CITICORP CRISIS — KEY FINDINGS',            [1 1 1],       13, 'bold';
        '',                                                 [1 1 1],       8,  'normal';
        'ORIGINAL DESIGN (1972):',                          [0.6 0.6 0.6], 10, 'bold';
        '  Welds specified, face wind only checked',        [0.6 0.6 0.6], 9,  'normal';
        '  Connections adequate (D/C < 1.0)',                [0.3 0.8 0.3], 9,  'normal';
        '',                                                 [1 1 1],       6,  'normal';
        'CONSTRUCTION CHANGE:',                             [1 0.7 0.3],   10, 'bold';
        '  Welds changed to bolts (cost savings)',          [1 0.7 0.3],   9,  'normal';
        '  No design recheck performed',                    [1 0.7 0.3],   9,  'normal';
        '',                                                 [1 1 1],       6,  'normal';
        'LEMESSURIER 1978 DISCOVERY:',                      [1 0.3 0.2],   10, 'bold';
        sprintf('  Quartering amp = x%.2f (static superpn)', amp_brace_1978), ...
                                                            [1 0.3 0.2],   9,  'normal';
        sprintf('  Bolt D/C = %.2f — OVERSTRESSED', DC_bolt_quar), ...
                                                            [1 0.3 0.2],   9,  'bold';
        '  Emergency repairs (welded cover plates)',         [1 0.3 0.2],   9,  'normal';
        '',                                                 [1 1 1],       6,  'normal';
        'NIST 2021 REASSESSMENT:',                          [0.3 1 0.3],   10, 'bold';
        '  Face winds govern (not quartering)',              [0.3 1 0.3],   9,  'normal';
        '  Static superposition was conservative',           [0.3 1 0.3],   9,  'normal';
        '  But bolts were still under-designed',             [0.8 0.8 0.2], 9,  'bold';
    };

    y = 0.95;
    for i = 1:size(lines, 1)
        if ~isempty(lines{i,1})
            text(0.03, y, lines{i,1}, 'Color', lines{i,2},...
                 'FontSize', lines{i,3}, 'FontWeight', lines{i,4},...
                 'Units', 'normalized', 'FontName', 'FixedWidth');
        end
        y = y - 0.052;
    end

    sgtitle('Module 10: 1978-Era vs Modern Wind Analysis — The Citicorp Crisis',...
            'Color','w','FontSize',16,'FontWeight','bold');
end

%% ===== HELPER =====
function s = status_str(dc)
    if dc > 1.0
        s = '<-- OVERSTRESSED';
    elseif dc > 0.9
        s = '<-- marginal';
    else
        s = '(OK)';
    end
end
