function [fig, mc_results] = monte_carlo_reliability(params, fea_results, wind_results)
%MONTE_CARLO_RELIABILITY  16-scenario reliability analysis of Citicorp Center
%  Full factorial design: 2 methodologies x 2 connections x 2 wind dirs x 2 TMD states
%
%  METHODOLOGY (no circular reasoning):
%    1. Gumbel calibrated from ASCE 7 wind map data (V_50, V_100)
%    2. Failure thresholds DERIVED from FEA brace forces + connection capacity
%    3. Historical values (LeMessurier T=16/55 yr) used only for VALIDATION
%
%  FACTORIAL DIMENSIONS:
%    Methodology:  1970s (ANSI A58.1, quartering=1.40x) vs Modern (ASCE 7-22, FEA)
%    Connection:   Welded (CJP) vs Bolted (A325-N)
%    Direction:    Perpendicular (0 deg) vs Quartering (45 deg)
%    TMD:          ON (8% damping) vs OFF (1% damping)
%
%  Gumbel CDF: F(v) = exp(-exp(-(v - mu)/beta))

    %% ===== GUMBEL DISTRIBUTION CALIBRATION =====
    V_50  = 100;   % mph (50-yr MRI, NYC Building Code)
    V_100 = 110;   % mph (100-yr MRI, estimated from ASCE 7)

    y_50  = -log(-log(1 - 1/50));
    y_100 = -log(-log(1 - 1/100));

    beta_gumbel = (V_100 - V_50) / (y_100 - y_50);
    mu_gumbel   = V_50 - beta_gumbel * y_50;

    fprintf('    Gumbel calibration (from ASCE 7 wind map data):\n');
    fprintf('      mu = %.2f mph, beta = %.2f mph\n', mu_gumbel, beta_gumbel);
    fprintf('      Calibration points: V_50=%.0f, V_100=%.0f (independent of failure data)\n', V_50, V_100);

    V_10  = mu_gumbel - beta_gumbel * log(-log(1 - 1/10));
    V_500 = mu_gumbel - beta_gumbel * log(-log(1 - 1/500));
    fprintf('      Implied V_10=%.0f, V_50=%.0f, V_500=%.0f mph\n', V_10, V_50, V_500);

    %% ===== 1972 WIND PROFILE SCALING =====
    % Compute 1970s ANSI A58.1 base shear to scale FEA forces
    % (Same formulas as analyze_1978_comparison.m lines 26-38)
    z_wind   = wind_results.z;
    zg       = 1200;       % ft — gradient height (Exposure B, 1972)
    alpha72  = 4.5;        % power-law exponent (ANSI A58.1-1972)
    G_1972   = 1.1;        % gust factor (approx for tall buildings)
    Cp_ww    = 0.8;
    Cp_lw    = 0.5;

    Kz_1972  = 2.01 * (max(z_wind, 15) / zg).^(2 / alpha72);
    q_1972   = 0.00256 * Kz_1972 * params.V_design^2;  % no Kd
    p_face_72 = q_1972 * G_1972 * (Cp_ww + Cp_lw);     % net face pressure
    V_face_72 = trapz(z_wind, p_face_72 * params.width) / 1000;  % kips
    V_face_22 = wind_results.V_base_perp / 1000;        % kips (ASCE 7-22)
    scale_72  = V_face_72 / max(V_face_22, 1);

    fprintf('\n    1972 vs Modern wind scaling:\n');
    fprintf('      1972 ANSI face shear: %.0f kips\n', V_face_72);
    fprintf('      ASCE 7-22 face shear: %.0f kips\n', V_face_22);
    fprintf('      Scale factor (1972/modern): %.3f\n', scale_72);

    %% ===== DEMAND AND CAPACITY =====
    V_design = params.V_design;  % 100 mph
    phi      = 0.75;

    % FEA brace forces at V_design (ASCE 7-22 loads)
    F_perp_modern = max(abs(fea_results.brace_forces_perp));
    F_quar_modern = max(abs(fea_results.brace_forces_quar));

    % 1970s demands (scaled from FEA + static superposition for quartering)
    F_perp_1970s = F_perp_modern * scale_72;
    F_quar_1970s = F_perp_1970s * params.quartering_force_factor;  % x1.40

    % TMD dynamic amplification factor (background/resonant separation)
    % TMD damping only reduces the resonant component; background is
    % independent of damping. Use gust factor decomposition.
    if isfield(wind_results, 'B_sq') && isfield(wind_results, 'R_sq')
        B_sq_w = wind_results.B_sq;
        R_sq_w = wind_results.R_sq;
        gQ_w   = wind_results.gQ;
        gR_w   = wind_results.gR;
        Iz_w   = wind_results.Iz_bar;

        % Gf numerator with TMD on (as computed in wind module):
        Gf_num_tmd = 1 + 1.7 * Iz_w * sqrt(gQ_w^2 * B_sq_w + gR_w^2 * R_sq_w);

        % R² scales inversely with damping: R²_bare = R²_tmd * (beta_tmd / beta_bare)
        R_sq_bare = R_sq_w * (params.damping_tmd / params.damping_bare);
        Gf_num_bare = 1 + 1.7 * Iz_w * sqrt(gQ_w^2 * B_sq_w + gR_w^2 * R_sq_bare);

        damp_amp = Gf_num_bare / Gf_num_tmd;
    else
        % Fallback if wind_results doesn't have gust factor components
        damp_amp = sqrt(params.damping_tmd / params.damping_bare);
    end

    % Connection capacities
    weld_capacity = phi * params.Fy * params.A_brace;  % CJP groove weld (kips)

    % Bolt group: historical splice size (6 rows x 8 bolts)
    Ab  = params.bolt_Ab;
    Fnv = params.bolt_Fnv;
    Rn_per_bolt = phi * Fnv * Ab;
    n_bolts     = params.n_bolts_per_splice;  % 48 (from historical records)
    bolt_capacity = n_bolts * Rn_per_bolt;

    fprintf('\n    Connection capacities:\n');
    fprintf('      Weld (CJP): %.1f kips\n', weld_capacity);
    fprintf('      Bolt (%d bolts): %.1f kips (D/C=%.2f under modern perp.)\n',...
            n_bolts, bolt_capacity, F_perp_modern / max(bolt_capacity, 0.01));
    fprintf('      FEA demands: perp=%.1f, quar=%.1f kips (modern)\n', F_perp_modern, F_quar_modern);
    fprintf('      1970s demands: perp=%.1f, quar=%.1f kips (scaled+amp)\n', F_perp_1970s, F_quar_1970s);

    %% ===== 16-SCENARIO FACTORIAL =====
    % Dimension order: Methodology(2) x Connection(2) x Direction(2) x TMD(2)
    N_scenarios = 16;
    V_fail      = zeros(1, N_scenarios);
    labels      = cell(1, N_scenarios);
    short_labels = cell(1, N_scenarios);
    demands     = zeros(1, N_scenarios);
    capacities  = zeros(1, N_scenarios);

    method_names = {'1970s', 'Modern'};
    conn_names   = {'Welded', 'Bolted'};
    dir_names    = {'Perp', 'Quar'};
    tmd_names    = {'TMD ON', 'TMD OFF'};

    % Base demands: [method][direction]
    %   1970s perp, 1970s quar, modern perp, modern quar
    base_demand = [F_perp_1970s, F_quar_1970s; ...
                   F_perp_modern, F_quar_modern];

    s = 0;
    for m = 1:2        % 1=1970s, 2=Modern
        for c = 1:2    % 1=Welded, 2=Bolted
            for d = 1:2    % 1=Perp, 2=Quartering
                for t = 1:2   % 1=TMD ON, 2=TMD OFF
                    s = s + 1;

                    % Select demand
                    F_base = base_demand(m, d);

                    % TMD OFF amplification
                    if t == 2
                        F_demand = F_base * damp_amp;
                    else
                        F_demand = F_base;
                    end

                    % Select capacity
                    if c == 1
                        cap = weld_capacity;
                    else
                        cap = bolt_capacity;
                    end

                    % Failure wind speed
                    V_fail(s) = V_design * sqrt(cap / max(F_demand, 0.01));
                    V_fail(s) = min(V_fail(s), 300);

                    demands(s)    = F_demand;
                    capacities(s) = cap;
                    labels{s}     = sprintf('%s, %s, %s, %s',...
                        method_names{m}, conn_names{c}, dir_names{d}, tmd_names{t});
                    short_labels{s} = sprintf('%s-%s-%s-%s',...
                        method_names{m}(1:3), conn_names{c}(1:3), ...
                        dir_names{d}(1:3), tmd_names{t}(5:end));
                end
            end
        end
    end

    %% ===== MONTE CARLO SIMULATION =====
    N_samples = 100000;
    rng(42);
    U = rand(N_samples, 1);
    V_annual = mu_gumbel - beta_gumbel * log(-log(U));

    % Monte Carlo failure probabilities
    P_fail_mc = zeros(1, N_scenarios);
    for si = 1:N_scenarios
        P_fail_mc(si) = sum(V_annual >= V_fail(si)) / N_samples;
    end

    % Analytical (exact) Gumbel probabilities
    P_fail = 1 - exp(-exp(-(V_fail - mu_gumbel) / beta_gumbel));
    T_return = 1 ./ max(P_fail, 1e-20);

    % Reliability index
    beta_rel = zeros(1, N_scenarios);
    for si = 1:N_scenarios
        if P_fail(si) > 0 && P_fail(si) < 1
            beta_rel(si) = -sqrt(2) * erfinv(2*P_fail(si) - 1);
        else
            beta_rel(si) = Inf;
        end
    end

    %% ===== CONSOLE OUTPUT =====
    fprintf('\n    ====== 1970s METHODOLOGY (ANSI A58.1-1972, scale=%.3f, quartering=x%.2f) ======\n',...
            scale_72, params.quartering_force_factor);
    fprintf('    %-42s  V_fail  Pf          T_return     beta\n', 'Scenario');
    fprintf('    %s\n', repmat('-', 1, 90));
    for si = 1:8
        fprintf('    %-42s  %5.0f   %.4e   %9.0f yr  %.2f\n',...
                labels{si}, V_fail(si), P_fail(si), min(T_return(si), 999999), beta_rel(si));
    end

    fprintf('\n    ====== MODERN METHODOLOGY (ASCE 7-22, Kd=0.85, FEA quartering) ======\n');
    fprintf('    %-42s  V_fail  Pf          T_return     beta\n', 'Scenario');
    fprintf('    %s\n', repmat('-', 1, 90));
    for si = 9:16
        fprintf('    %-42s  %5.0f   %.4e   %9.0f yr  %.2f\n',...
                labels{si}, V_fail(si), P_fail(si), min(T_return(si), 999999), beta_rel(si));
    end

    % Validation against LeMessurier
    fprintf('\n    VALIDATION vs LeMessurier:\n');
    fprintf('      Scenario  8 (1970s, Bolt, Quar, TMD OFF): T=%.0f yr  vs  historical T=16 yr\n', T_return(8));
    fprintf('      Scenario  7 (1970s, Bolt, Quar, TMD ON):  T=%.0f yr  vs  historical T=55 yr\n', T_return(7));
    fprintf('      Scenario 16 (Modern, Bolt, Quar, TMD OFF): T=%.0f yr  (modern reassessment)\n', T_return(16));
    fprintf('      (Differences reflect independent calibration - NOT circular.)\n');

    %% ===== PLOTTING (8 panels, 2x4) =====
    fig = figure('Name','Monte Carlo Reliability — 16 Scenarios (2x2x2x2)',...
                 'Position',[30 20 1800 1000],'Color','k');

    % Reshape results into 4x4 heatmap matrix
    % Rows: (1970s-Weld, 1970s-Bolt, Modern-Weld, Modern-Bolt)
    % Cols: (Perp+TMD, Perp-TMD, Quar+TMD, Quar-TMD)
    T_heatmap    = zeros(4, 4);
    beta_heatmap = zeros(4, 4);
    Vf_heatmap   = zeros(4, 4);
    si = 0;
    for m = 1:2
        for c = 1:2
            row = 2*(m-1) + c;
            for d = 1:2
                for t = 1:2
                    si = si + 1;
                    col = 2*(d-1) + t;
                    T_heatmap(row, col)    = min(T_return(si), 99999);
                    beta_heatmap(row, col) = min(beta_rel(si), 10);
                    Vf_heatmap(row, col)   = V_fail(si);
                end
            end
        end
    end

    row_labels = {'1970s Weld', '1970s Bolt', 'Modern Weld', 'Modern Bolt'};
    col_labels = {'Perp+TMD', 'Perp-TMD', 'Quar+TMD', 'Quar-TMD'};

    % --- Panel 1: Return Period Heatmap ---
    ax1 = subplot(2,4,1);
    T_log = log10(max(T_heatmap, 1));  % log scale for color
    imagesc(T_log); hold on;
    colormap(ax1, risk_colormap());
    cb1 = colorbar('Color','w');
    cb1.Label.String = 'log_{10}(T)';
    cb1.Label.Color = 'w';
    clim([0.5 5]);
    set(ax1,'Color','k','XColor','w','YColor','w');
    set(ax1,'XTick',1:4,'XTickLabel',col_labels,'XTickLabelRotation',30);
    set(ax1,'YTick',1:4,'YTickLabel',row_labels);
    for r = 1:4
        for ci = 1:4
            T_val = T_heatmap(r, ci);
            if T_val >= 10000
                str = '>10k';
            else
                str = sprintf('%.0f', T_val);
            end
            % Use dark text on light cells, white on dark
            if T_log(r,ci) > 3
                tcol = 'k';
            else
                tcol = 'w';
            end
            text(ci, r, str, 'Color', tcol, 'FontSize', 9, 'FontWeight', 'bold',...
                 'HorizontalAlignment', 'center');
        end
    end
    title('Return Period (years)','Color','w','FontSize',11);

    % --- Panel 2: Reliability Index Heatmap ---
    ax2 = subplot(2,4,2);
    imagesc(beta_heatmap); hold on;
    colormap(ax2, beta_colormap());
    cb2 = colorbar('Color','w');
    cb2.Label.String = '\beta';
    cb2.Label.Color = 'w';
    clim([0 6]);
    set(ax2,'Color','k','XColor','w','YColor','w');
    set(ax2,'XTick',1:4,'XTickLabel',col_labels,'XTickLabelRotation',30);
    set(ax2,'YTick',1:4,'YTickLabel',row_labels);
    for r = 1:4
        for ci = 1:4
            bval = beta_heatmap(r, ci);
            if bval >= 10
                str = '>10';
            else
                str = sprintf('%.1f', bval);
            end
            if bval > 4
                tcol = 'k';
            else
                tcol = 'w';
            end
            text(ci, r, str, 'Color', tcol, 'FontSize', 9, 'FontWeight', 'bold',...
                 'HorizontalAlignment', 'center');
        end
    end
    % Add beta=3.0 reference contour note
    title('Reliability Index \beta','Color','w','FontSize',11);

    % --- Panel 3: Gumbel PDF with key thresholds ---
    ax3 = subplot(2,4,3); hold on;
    set(ax3,'Color','k','XColor','w','YColor','w');

    v_lo = max(0, mu_gumbel - 3*beta_gumbel);
    v_hi = mu_gumbel + 8*beta_gumbel;
    edges = linspace(v_lo, v_hi, 80);
    histogram(V_annual, edges, 'FaceColor', [0.3 0.3 0.6], 'EdgeColor', 'none',...
              'FaceAlpha', 0.7, 'Normalization', 'pdf');

    v_plot = linspace(v_lo, v_hi, 500);
    z_g = (v_plot - mu_gumbel) / beta_gumbel;
    pdf_g = (1/beta_gumbel) .* exp(-(z_g + exp(-z_g)));
    plot(v_plot, pdf_g, '-w', 'LineWidth', 2);

    % Key scenarios only (4 highlighted)
    key_idx    = [8, 7, 16, 9];  % worst 1970s, 1970s+TMD, worst modern, best
    key_colors = {[1 0.2 0.2], [1 0.7 0.1], [0.8 0.4 1], [0.2 0.8 0.2]};
    key_names  = {'70s-B-Q-off', '70s-B-Q-on', 'Mod-B-Q-off', 'Mod-W-P-on'};
    for ki = 1:length(key_idx)
        si = key_idx(ki);
        if V_fail(si) < v_hi
            xline(V_fail(si), '-', 'Color', key_colors{ki}, 'LineWidth', 2,...
                  'Label', sprintf('%.0f', V_fail(si)),...
                  'LabelColor', key_colors{ki}, 'FontSize', 7,...
                  'LabelHorizontalAlignment', 'left',...
                  'LabelVerticalAlignment', 'top');
        end
    end

    % Historical markers
    xline(params.V_quartering_fail_noTMD, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1,...
          'Label', 'Hist.70', 'LabelColor', [0.5 0.5 0.5], 'FontSize', 7);
    xline(params.V_quartering_fail_TMD, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1,...
          'Label', 'Hist.85', 'LabelColor', [0.5 0.5 0.5], 'FontSize', 7);

    xlabel('Annual Max Wind (mph)','Color','w');
    ylabel('PDF','Color','w');
    title('Gumbel PDF + Key V_{fail}','Color','w','FontSize',11);
    grid on; set(ax3,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 4: CDF Exceedance ---
    ax4 = subplot(2,4,4); hold on;
    set(ax4,'Color','k','XColor','w','YColor','w');

    cdf_g = exp(-exp(-z_g));
    plot(v_plot, 1 - cdf_g, '-w', 'LineWidth', 2);

    for ki = 1:length(key_idx)
        si = key_idx(ki);
        if P_fail(si) > 1e-10
            plot(V_fail(si), P_fail(si), 'o', 'Color', key_colors{ki},...
                 'MarkerSize', 10, 'MarkerFaceColor', key_colors{ki}, 'LineWidth', 2);
            text(V_fail(si)+2, P_fail(si)*1.5, key_names{ki},...
                 'Color', key_colors{ki}, 'FontSize', 8);
        end
    end

    set(ax4, 'YScale', 'log');
    xlabel('Wind Speed (mph)','Color','w');
    ylabel('P(V_{ann} \geq v)','Color','w');
    title('Exceedance Probability','Color','w','FontSize',11);
    ylim([1e-5 1]);
    grid on; set(ax4,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 5: 1970s vs Modern Grouped Bar Chart ---
    ax5 = subplot(2,4,5); hold on;
    set(ax5,'Color','k','XColor','w','YColor','w');

    % Show bolted scenarios (the interesting ones) side by side
    % 1970s bolted: scenarios 5-8; Modern bolted: scenarios 13-16
    T_1970s_bolt = min(T_return(5:8), 5000);
    T_modern_bolt = min(T_return(13:16), 5000);
    x = 1:4;
    bw = 0.35;
    bar(x - bw/2, T_1970s_bolt, bw, 'FaceColor', [1 0.4 0.3], 'EdgeColor', 'none');
    bar(x + bw/2, T_modern_bolt, bw, 'FaceColor', [0.3 0.7 1], 'EdgeColor', 'none');

    % Value labels
    for i = 1:4
        text(x(i)-bw/2, T_1970s_bolt(i)+max([T_1970s_bolt T_modern_bolt])*0.03,...
             sprintf('%.0f', T_return(4+i)),...
             'Color',[1 0.6 0.4],'FontSize',8,'FontWeight','bold','HorizontalAlignment','center');
        text(x(i)+bw/2, T_modern_bolt(i)+max([T_1970s_bolt T_modern_bolt])*0.03,...
             sprintf('%.0f', T_return(12+i)),...
             'Color',[0.5 0.8 1],'FontSize',8,'FontWeight','bold','HorizontalAlignment','center');
    end

    yline(50, ':', 'Color', [0.3 0.8 0.3], 'LineWidth', 1,...
          'Label', '50-yr target', 'LabelColor', [0.3 0.8 0.3]);
    yline(16, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1,...
          'Label', 'Hist. T=16', 'LabelColor', [0.5 0.5 0.5]);

    set(ax5,'XTick',1:4,'XTickLabel',{'P+TMD','P-TMD','Q+TMD','Q-TMD'});
    ylabel('Return Period (yr)','Color','w');
    title('Bolted: 1970s vs Modern','Color','w','FontSize',11);
    legend({'1970s','Modern'},'TextColor','w','Color',[0.15 0.15 0.15],'Location','northwest');
    grid on; set(ax5,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 6: MC Convergence (worst case = Scenario 8) ---
    ax6 = subplot(2,4,6); hold on;
    set(ax6,'Color','k','XColor','w','YColor','w');

    worst_s = 8;  % 1970s-Bolted-Quar-TMD OFF
    if P_fail(worst_s) > 0
        n_pts = 500;
        n_vec = round(logspace(2, log10(N_samples), n_pts));
        Pf_running = zeros(size(n_vec));
        failures_w = V_annual >= V_fail(worst_s);
        for i = 1:length(n_vec)
            Pf_running(i) = sum(failures_w(1:n_vec(i))) / n_vec(i);
        end

        plot(n_vec, Pf_running, '-', 'Color', [1 0.3 0.3], 'LineWidth', 1.5);
        yline(P_fail(worst_s), '--w', 'LineWidth', 1.5);

        ci_u = P_fail(worst_s) + 1.96*sqrt(P_fail(worst_s)*(1-P_fail(worst_s))./n_vec);
        ci_l = P_fail(worst_s) - 1.96*sqrt(P_fail(worst_s)*(1-P_fail(worst_s))./n_vec);
        fill([n_vec fliplr(n_vec)], [ci_u fliplr(ci_l)],...
             [1 0.3 0.3], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

        set(ax6, 'XScale', 'log');
        xlabel('Samples','Color','w');
        ylabel('P_f','Color','w');
        title(sprintf('MC Convergence (Sc.%d)', worst_s),'Color','w','FontSize',11);
        legend({'MC','Exact','95% CI'},'TextColor','w','Color',[0.15 0.15 0.15],'Location','northeast');
    end
    grid on; set(ax6,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 7: Validation (Derived vs Historical) ---
    ax7 = subplot(2,4,7); hold on;
    set(ax7,'Color','k','XColor','w','YColor','w');

    % Compare derived vs LeMessurier for the two critical scenarios
    derived_T = [T_return(8), T_return(7)];
    histor_T  = [16, 55];
    x = 1:2;
    bw = 0.3;
    bar(x - bw/2, min(derived_T, 200), bw, 'FaceColor', [0.3 0.7 1], 'EdgeColor', 'none');
    bar(x + bw/2, histor_T, bw, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'none');

    for i = 1:2
        text(x(i)-bw/2, min(derived_T(i),200)+5, sprintf('%.0f', derived_T(i)),...
             'Color',[0.5 0.8 1],'FontSize',10,'FontWeight','bold','HorizontalAlignment','center');
        text(x(i)+bw/2, histor_T(i)+5, sprintf('%.0f', histor_T(i)),...
             'Color',[1 0.6 0.4],'FontSize',10,'FontWeight','bold','HorizontalAlignment','center');
    end

    set(ax7,'XTick',1:2,'XTickLabel',{'Bolt+Q-TMD','Bolt+Q+TMD'});
    ylabel('Return Period (yr)','Color','w');
    title('Derived vs Historical','Color','w','FontSize',11);
    legend({'Derived (this sim)','LeMessurier 1978'},...
           'TextColor','w','Color',[0.15 0.15 0.15],'Location','northeast');
    grid on; set(ax7,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 8: Summary Text ---
    ax8 = subplot(2,4,8);
    set(ax8,'Color','k','XColor','k','YColor','k','XTick',[],'YTick',[]);
    box off; axis off;

    summary = {
        '\bf{16-Scenario Reliability Summary}'
        ''
        sprintf('Gumbel: \\mu=%.1f, \\beta_G=%.1f mph', mu_gumbel, beta_gumbel)
        sprintf('Calibration: V_{50}=%d, V_{100}=%d', V_50, V_100)
        ''
        sprintf('Weld cap: %.0f kips', weld_capacity)
        sprintf('Bolt cap: %.0f kips (%d bolts)', bolt_capacity, n_bolts)
        sprintf('1972 scale: %.3f', scale_72)
        sprintf('TMD amp: %.2f', damp_amp)
        ''
        '\bf{CRITICAL SCENARIOS:}'
        sprintf(' 8: 70s-B-Q-off  T=%5.0f yr', T_return(8))
        sprintf(' 7: 70s-B-Q-on   T=%5.0f yr', T_return(7))
        sprintf('16: Mod-B-Q-off  T=%5.0f yr', T_return(16))
        sprintf('15: Mod-B-Q-on   T=%5.0f yr', T_return(15))
        ''
        '\bf{VALIDATION:}'
        sprintf(' Sc 8 vs hist: %.0f vs 16 yr', T_return(8))
        sprintf(' Sc 7 vs hist: %.0f vs 55 yr', T_return(7))
    };
    text(0.02, 0.97, summary, 'Color', 'w', 'FontSize', 8,...
         'VerticalAlignment', 'top', 'FontName', 'FixedWidth',...
         'Interpreter', 'tex');

    sgtitle('Module 5: 16-Scenario Monte Carlo Reliability (1970s vs Modern)',...
            'Color','w','FontSize',15,'FontWeight','bold');

    %% ===== EXPORT JSON FOR SLIDES =====
    json_struct = struct();
    json_struct.N_scenarios = N_scenarios;
    json_struct.mu_gumbel = mu_gumbel;
    json_struct.beta_gumbel = beta_gumbel;
    json_struct.scale_72 = scale_72;
    json_struct.damp_amp = damp_amp;
    json_struct.weld_capacity = weld_capacity;
    json_struct.bolt_capacity = bolt_capacity;
    json_struct.n_bolts = n_bolts;
    json_struct.V_fail = V_fail;
    json_struct.P_fail = P_fail;
    json_struct.T_return = T_return;
    json_struct.beta_rel = beta_rel;
    json_struct.labels = labels;
    json_struct.short_labels = short_labels;

    json_path = fullfile(pwd, 'figures', 'mc_results_16.json');
    try
        json_str = jsonencode(json_struct);
        fid = fopen(json_path, 'w');
        if fid > 0
            fprintf(fid, '%s', json_str);
            fclose(fid);
            fprintf('    JSON results exported: %s\n', json_path);
        end
    catch
        fprintf('    (JSON export skipped — jsonencode not available)\n');
    end

    %% ===== PACK OUTPUT STRUCT =====
    mc_results = struct();
    mc_results.V_fail     = V_fail;
    mc_results.P_fail     = P_fail;
    mc_results.T_return   = T_return;
    mc_results.beta_rel   = beta_rel;
    mc_results.labels     = labels;
    mc_results.short_labels = short_labels;
    mc_results.N_scenarios  = N_scenarios;
    mc_results.mu_gumbel    = mu_gumbel;
    mc_results.beta_gumbel  = beta_gumbel;
    mc_results.scale_72     = scale_72;
    mc_results.damp_amp     = damp_amp;
    mc_results.weld_capacity  = weld_capacity;
    mc_results.bolt_capacity  = bolt_capacity;
    mc_results.n_bolts        = n_bolts;
end

%% ===== LOCAL HELPER FUNCTIONS =====

function cmap = risk_colormap()
%RISK_COLORMAP  Red (low T) -> Yellow -> Green (high T)
    n = 256;
    r = [linspace(0.8, 1, n/2), linspace(1, 0.2, n/2)];
    g = [linspace(0.1, 1, n/2), linspace(1, 0.8, n/2)];
    b = [linspace(0.1, 0.2, n/2), linspace(0.2, 0.2, n/2)];
    cmap = [r', g', b'];
end

function cmap = beta_colormap()
%BETA_COLORMAP  Red (low beta) -> Yellow (beta~3) -> Green (high beta)
    n = 256;
    r = [linspace(0.9, 1, n/2), linspace(1, 0.1, n/2)];
    g = [linspace(0.1, 1, n/2), linspace(1, 0.7, n/2)];
    b = [linspace(0.1, 0.1, n/2), linspace(0.1, 0.3, n/2)];
    cmap = [r', g', b'];
end
