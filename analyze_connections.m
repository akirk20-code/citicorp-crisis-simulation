function fig = analyze_connections(params, fea_results)
%ANALYZE_CONNECTIONS  AISC 360-22 connection capacity analysis
%  Compares demand (from FEA) vs capacity for bolted and welded connections.
%  Uses FEA-derived quartering forces (no hardcoded amplification factor).
%  Implements bolt combined tension-shear interaction per AISC Eq. J3-3a.

    %% ===== AISC 360-22 CONNECTION CAPACITIES =====

    % --- A325-N Bolt Capacity (per bolt) ---
    Ab   = params.bolt_Ab;       % 0.785 in^2 (1" diameter)
    Fnv  = params.bolt_Fnv;      % 54 ksi — nominal shear
    Fnt  = params.bolt_Fnt;      % 90 ksi — nominal tensile
    phi  = 0.75;                 % resistance factor (LRFD)

    Rn_shear_bolt = phi * Fnv * Ab;    % kips per bolt
    Rn_tension_bolt = phi * Fnt * Ab;  % kips per bolt

    fprintf('    A325-N Bolt (d=%.2f in):\n', params.bolt_diam);
    fprintf('      phi*Rn_v per bolt  = %.1f kips\n', Rn_shear_bolt);
    fprintf('      phi*Rn_t per bolt  = %.1f kips\n', Rn_tension_bolt);

    % --- E70XX Fillet Weld Capacity ---
    FEXX = params.weld_FEXX;     % 70 ksi
    w    = params.weld_size;     % 0.5 in leg
    Lw   = params.weld_length;   % 12 in per side
    te   = w * 0.707;            % effective throat
    phi_w = 0.75;
    Rn_weld_per_inch = phi_w * 0.60 * FEXX * te;  % kips/in
    Rn_weld_total = 2 * Rn_weld_per_inch * Lw;

    % Full penetration groove weld (CJP) — capacity = base metal
    Rn_CJP = phi * params.Fy * params.A_brace;  % kips

    fprintf('    CJP Groove Weld:  phi*Rn = %.1f kips (base metal governs)\n', Rn_CJP);

    % --- Repaired Connection (2" welded cover plates) ---
    Lw_repair = 18;  % inches per side
    Rn_repair = Rn_CJP + 4 * Rn_weld_per_inch * Lw_repair;
    fprintf('    Repaired Connection: %.1f kips\n', Rn_repair);

    %% ===== DEMAND FROM FEA =====
    % Perpendicular brace forces (from 3D FEA)
    F_brace_perp = abs(fea_results.brace_forces_perp);  % kips

    % Quartering brace forces — DERIVED directly from FEA (no hardcoded factor)
    F_brace_quar = abs(fea_results.brace_forces_quar);  % kips

    % DERIVED amplification factor (compare with historical 1.40)
    amp_derived = max(F_brace_quar) / max(max(F_brace_perp), 0.001);
    amp_historical = params.quartering_force_factor;  % 1.40 (LeMessurier)

    fprintf('    Quartering brace forces (from FEA, no hardcoded factor):\n');
    fprintf('      Max perp. force:       %.1f kips\n', max(F_brace_perp));
    fprintf('      Max quar. force (FEA): %.1f kips\n', max(F_brace_quar));
    fprintf('      DERIVED amplification: %.2f (quartering / perpendicular)\n', amp_derived);
    fprintf('      VALIDATION: historical = %.2f (LeMessurier 1978)\n', amp_historical);

    % Brace geometry — diagonal angle affects force decomposition at gusset
    theta = atand(params.tier_stories * params.story_height / (params.width/2));
    fprintf('    Brace angle from horizontal: %.1f deg\n', theta);

    %% ===== AUTO-SIZE BOLT GROUP FOR PERPENDICULAR =====
    % The real story: bolt groups were designed for perpendicular wind only
    % (1970 NYC code). They PASSED the perpendicular check (D/C ~ 0.90).
    % The quartering check was never done until 1978.

    F_max_perp = max(F_brace_perp);
    F_max_quar = max(F_brace_quar);
    target_DC = 0.90;

    n_bolts = max(4, ceil(F_max_perp / (target_DC * Rn_shear_bolt)));
    Rn_shear_group = n_bolts * Rn_shear_bolt;
    Rn_tension_group = n_bolts * Rn_tension_bolt;

    fprintf('    Auto-sized bolt group: %d bolts (D/C=%.2f under perp.)\n',...
            n_bolts, target_DC);
    fprintf('    Max brace force (perp.):   %.1f kips\n', F_max_perp);
    fprintf('    Max brace force (quart.):  %.1f kips (FEA-derived, amp=%.2f)\n',...
            F_max_quar, amp_derived);

    %% ===== D/C RATIOS WITH AISC J3-3a INTERACTION =====
    n_braces = length(F_brace_perp);

    % At the gusset plate, the diagonal brace force resolves into:
    %   Shear component (horizontal): F * cos(theta)
    %   Tension component (vertical): F * sin(theta)
    % The bolts see COMBINED tension + shear at this connection.

    % Pre-allocate
    DC_bolt_perp  = zeros(n_braces, 1);  % J3-3a combined interaction
    DC_bolt_quar  = zeros(n_braces, 1);
    DC_weld_perp  = zeros(n_braces, 1);
    DC_weld_quar  = zeros(n_braces, 1);
    DC_repair_perp = zeros(n_braces, 1);
    DC_repair_quar = zeros(n_braces, 1);

    % Store stress demands for plotting
    frv_perp = zeros(n_braces, 1);  % shear stress on bolt
    frt_perp = zeros(n_braces, 1);  % tension stress on bolt
    frv_quar = zeros(n_braces, 1);
    frt_quar = zeros(n_braces, 1);

    for b = 1:n_braces
        % --- PERPENDICULAR: bolt interaction (AISC J3-3a) ---
        V_p = F_brace_perp(b) * cosd(theta);   % shear demand at gusset
        T_p = F_brace_perp(b) * sind(theta);   % tension demand at gusset
        fv_p = V_p / (n_bolts * Ab);           % required shear stress
        ft_p = T_p / (n_bolts * Ab);           % required tension stress

        % Eq. J3-3a: reduced available tension under combined loading
        Fnt_prime_p = 1.3*Fnt - (Fnt/(phi*Fnv))*fv_p;
        Fnt_prime_p = max(min(Fnt_prime_p, Fnt), 0);

        % D/C = max of shear ratio and tension interaction ratio
        dc_shear_p = fv_p / (phi * Fnv);
        dc_tension_p = ft_p / max(phi * Fnt_prime_p, 0.001);
        DC_bolt_perp(b) = max(dc_shear_p, dc_tension_p);

        frv_perp(b) = fv_p;
        frt_perp(b) = ft_p;

        % --- QUARTERING: bolt interaction (amplified forces) ---
        V_q = F_brace_quar(b) * cosd(theta);
        T_q = F_brace_quar(b) * sind(theta);
        fv_q = V_q / (n_bolts * Ab);
        ft_q = T_q / (n_bolts * Ab);

        Fnt_prime_q = 1.3*Fnt - (Fnt/(phi*Fnv))*fv_q;
        Fnt_prime_q = max(min(Fnt_prime_q, Fnt), 0);

        dc_shear_q = fv_q / (phi * Fnv);
        dc_tension_q = ft_q / max(phi * Fnt_prime_q, 0.001);
        DC_bolt_quar(b) = max(dc_shear_q, dc_tension_q);

        frv_quar(b) = fv_q;
        frt_quar(b) = ft_q;

        % --- CJP Groove Weld D/C ---
        DC_weld_perp(b) = F_brace_perp(b) / Rn_CJP;
        DC_weld_quar(b) = F_brace_quar(b) / Rn_CJP;

        % --- Repaired Connection D/C ---
        DC_repair_perp(b) = F_brace_perp(b) / Rn_repair;
        DC_repair_quar(b) = F_brace_quar(b) / Rn_repair;
    end

    % Summary statistics
    fprintf('\n    D/C Ratio Summary (AISC J3-3a interaction, max across braces):\n');
    fprintf('    %-30s  Perp.    Quart.(FEA)\n', 'Connection Type');
    fprintf('    %s\n', repmat('-', 1, 60));
    fprintf('    %-30s  %.3f    %.3f\n', 'A325 Bolted (J3-3a)', max(DC_bolt_perp), max(DC_bolt_quar));
    fprintf('    %-30s  %.3f    %.3f\n', 'CJP Groove Weld', max(DC_weld_perp), max(DC_weld_quar));
    fprintf('    %-30s  %.3f    %.3f\n', 'Repaired (cover plates)', max(DC_repair_perp), max(DC_repair_quar));

    n_fail = sum(DC_bolt_quar >= 1.0);
    fprintf('    Braces exceeding bolt capacity under quartering: %d of %d\n', n_fail, n_braces);

    %% ===== PLOTTING =====
    fig = figure('Name','Connection Analysis (AISC 360-22)',...
                 'Position',[50 30 1500 900],'Color','k');

    % --- Panel 1: D/C bar chart ---
    ax1 = subplot(2,3,1); hold on;
    set(ax1,'Color','k','XColor','w','YColor','w');

    categories = {'Bolted\newline(perp.)','Bolted\newline(quart.)',...
                  'Welded\newline(perp.)','Welded\newline(quart.)',...
                  'Repaired\newline(perp.)','Repaired\newline(quart.)'};
    DC_vals = [max(DC_bolt_perp), max(DC_bolt_quar),...
               max(DC_weld_perp), max(DC_weld_quar),...
               max(DC_repair_perp), max(DC_repair_quar)];
    bar_colors = [0.2 0.6 1;  1 0.3 0.3;
                  0.2 0.6 1;  1 0.3 0.3;
                  0.2 0.6 1;  1 0.3 0.3];

    b1 = bar(1:6, DC_vals, 0.7);
    b1.FaceColor = 'flat';
    b1.CData = bar_colors;

    yline(1.0, '--', 'Color', [1 1 0], 'LineWidth', 2,...
          'Label', 'D/C = 1.0 (FAILURE)', 'LabelColor', [1 1 0],...
          'FontSize', 10, 'LabelHorizontalAlignment', 'left');

    for i = 1:6
        if DC_vals(i) >= 1.0
            text(i, DC_vals(i)+0.05, 'FAILS', 'Color', 'r',...
                 'FontSize', 10, 'FontWeight', 'bold',...
                 'HorizontalAlignment', 'center');
        else
            text(i, DC_vals(i)+0.05, sprintf('%.2f', DC_vals(i)),...
                 'Color', 'w', 'FontSize', 9,...
                 'HorizontalAlignment', 'center');
        end
    end

    set(ax1, 'XTick', 1:6, 'XTickLabel', categories, 'XTickLabelRotation', 0);
    ylabel('Demand / Capacity Ratio','Color','w');
    title('Max D/C Ratio by Connection Type','Color','w','FontSize',12);
    ylim([0, max(max(DC_vals)*1.3, 1.5)]);
    grid on; set(ax1,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 2: D/C distribution across braces (bolted, quartering) ---
    ax2 = subplot(2,3,2); hold on;
    set(ax2,'Color','k','XColor','w','YColor','w');

    DC_sorted = sort(DC_bolt_quar, 'descend');
    barh(1:n_braces, DC_sorted, 'FaceColor', [1 0.3 0.3], 'EdgeColor', 'none');
    xline(1.0, '--', 'Color', [1 1 0], 'LineWidth', 2);
    text(max(DC_sorted)*0.5, n_braces*0.1,...
         sprintf('%d of %d braces\nexceed capacity', n_fail, n_braces),...
         'Color', [1 1 0], 'FontSize', 11, 'FontWeight', 'bold');

    xlabel('D/C Ratio (Bolted, Quartering)','Color','w');
    ylabel('Brace # (sorted)','Color','w');
    title('Bolted D/C — Quartering Wind','Color','w','FontSize',12);
    grid on; set(ax2,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 3: Bolt interaction diagram (AISC J3-3a) ---
    ax3 = subplot(2,3,3); hold on;
    set(ax3,'Color','k','XColor','w','YColor','w');

    % AISC J3-3a interaction envelope
    fv_range = linspace(0, Fnv, 100);
    ft_limit = zeros(size(fv_range));
    for i = 1:length(fv_range)
        ft_limit(i) = 1.3*Fnt - (Fnt/(phi*Fnv))*fv_range(i);
        ft_limit(i) = min(ft_limit(i), Fnt);
        ft_limit(i) = max(ft_limit(i), 0);
    end

    fill([fv_range fliplr(fv_range)], [ft_limit zeros(size(ft_limit))],...
         [0.2 0.4 0.2], 'FaceAlpha', 0.3, 'EdgeColor', [0.3 0.8 0.3], 'LineWidth', 2);

    % Demand points
    plot(frv_perp, frt_perp, 'o', 'Color', [0.2 0.6 1], 'MarkerSize', 6,...
         'MarkerFaceColor', [0.2 0.6 1]);
    plot(frv_quar, frt_quar, 's', 'Color', [1 0.3 0.3], 'MarkerSize', 8,...
         'MarkerFaceColor', [1 0.3 0.3]);

    % Mark points outside the envelope
    has_exceed = false;
    for b = 1:n_braces
        fv_q = frv_quar(b);
        ft_q = frt_quar(b);
        Fnt_prime_check = min(1.3*Fnt - (Fnt/(phi*Fnv))*fv_q, Fnt);
        if ft_q > phi * Fnt_prime_check || fv_q > phi * Fnv
            plot(fv_q, ft_q, 'x', 'Color', [1 1 0], 'MarkerSize', 14, 'LineWidth', 3);
            has_exceed = true;
        end
    end

    xlabel('Required Shear Stress f_{rv} (ksi)','Color','w');
    ylabel('Required Tensile Stress f_{rt} (ksi)','Color','w');
    title('AISC J3-3a Bolt Interaction','Color','w','FontSize',12);
    leg_entries = {'Capacity envelope','Perp. wind','Quart. wind (FEA)'};
    if has_exceed, leg_entries{end+1} = 'Exceeds capacity'; end
    legend(leg_entries, 'TextColor','w','Color',[0.15 0.15 0.15],'Location','northeast');
    grid on; set(ax3,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 4: Force amplification visualization ---
    ax4 = subplot(2,3,4); hold on;
    set(ax4,'Color','k','XColor','w','YColor','w');

    % Normalize forces and capacities to perpendicular demand
    F_norm_perp = 1.0;
    F_norm_quar = F_max_quar / F_max_perp;
    C_weld_norm = Rn_CJP / F_max_perp;
    C_bolt_norm = Rn_shear_group / F_max_perp;
    C_repair_norm = Rn_repair / F_max_perp;

    x_pos = [1, 2, 4, 5, 6];
    heights = [F_norm_perp, F_norm_quar, C_weld_norm, C_bolt_norm, C_repair_norm];
    bar_c = [0.2 0.6 1; 1 0.3 0.3; 0.2 0.8 0.2; 1 0.7 0.1; 0.5 0.8 1];

    for i = 1:5
        bar(x_pos(i), heights(i), 0.6, 'FaceColor', bar_c(i,:));
    end

    yline(1.0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

    set(ax4, 'XTick', x_pos, 'XTickLabel', ...
        {'Perp.',sprintf('Quart.\n(%.0f%%)',amp_derived*100-100),'Weld','Bolt','Repair'});
    ylabel('Normalized to Perp. Demand','Color','w');
    title('Demand vs Capacity','Color','w','FontSize',12);
    grid on; set(ax4,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 5-6: Crisis narrative ---
    ax5 = subplot(2,3,[5 6]);
    set(ax5,'Color','k','XColor','k','YColor','k','XTick',[],'YTick',[]);
    box off; axis off;

    narrative = {
        '\bf{THE CITICORP CRISIS — Connection Failure Mechanism}'
        ''
        '\bf{Step 1:} Original design specified \bf{CJP welded} connections'
        sprintf('  Capacity: %.0f kips — adequate for ALL wind angles', Rn_CJP)
        ''
        '\bf{Step 2:} Bethlehem Steel substituted \bf{A325 bolted} splices'
        sprintf('  Bolt group: %d bolts, capacity = %.0f kips (%.0f%% of weld)', ...
                n_bolts, Rn_shear_group, 100*Rn_shear_group/Rn_CJP)
        ''
        '\bf{Step 3:} 1970 NYC code only checked \bf{perpendicular} wind'
        sprintf('  Perp. D/C = %.2f (AISC J3-3a)  \\rightarrow  PASSES', max(DC_bolt_perp))
        ''
        '\bf{Step 4:} LeMessurier (1978) checked \bf{quartering winds}'
        sprintf('  Quartering amplifies critical brace forces by %.0f%% (FEA-derived)', (amp_derived-1)*100)
        sprintf('  Quart. D/C = %.2f (J3-3a)  \\rightarrow  \\bf{FAILS}', max(DC_bolt_quar))
        ''
        '\bf{Step 5:} Emergency repair — 2" welded cover plates'
        sprintf('  Repaired D/C = %.2f  \\rightarrow  SAFE', max(DC_repair_quar))
        ''
        '\bf{AISC 360-22:} Eq. J3-3a combined tension-shear interaction'
        sprintf('  A325-N: F_{nv}=%d ksi, F_{nt}=%d ksi, \\theta=%.0f\\circ', Fnv, Fnt, theta)
    };

    text(0.02, 0.97, narrative, 'Color', 'w', 'FontSize', 10,...
         'VerticalAlignment', 'top', 'Interpreter', 'tex',...
         'FontName', 'FixedWidth');

    sgtitle('Module 4: AISC 360-22 Connection Analysis — Citicorp Center',...
            'Color','w','FontSize',16,'FontWeight','bold');
end
