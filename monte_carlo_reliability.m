function fig = monte_carlo_reliability(params)
%MONTE_CARLO_RELIABILITY  Reliability analysis of Citicorp Center
%  Uses Gumbel Type I extreme value distribution for annual max wind speed.
%  Compares failure probabilities across 4 scenarios.
%
%  Calibration: Gumbel parameters derived from ASCE 7 wind speeds
%    V_50  = 100 mph  (50-year MRI, ~1970 code)
%    V_100 = 110 mph  (100-year MRI, estimated)
%  Gumbel CDF: F(v) = exp(-exp(-(v - mu)/beta))
%    => V_N = mu - beta * ln(-ln(1 - 1/N))

    %% ===== GUMBEL DISTRIBUTION CALIBRATION =====
    % Calibrate Gumbel to match LeMessurier's 1978 analysis:
    %   - P(V >= 70 mph) = 1/16  (no TMD, quartering failure)
    %   - P(V >= 85 mph) = 1/55  (with TMD, quartering failure)
    % These are "effective wind speeds at building reference height" using
    % 1970s fastest-mile conventions. The Gumbel represents the annual
    % maximum of this effective wind speed in NYC.
    %
    % Gumbel CDF: F(v) = exp(-exp(-(v - mu)/beta))
    % P(V >= v) = 1 - F(v) = 1 - exp(-exp(-(v-mu)/beta))
    %
    % Two equations from LeMessurier's return periods:
    %   (V_fail3 - mu)/beta = -ln(-ln(1 - 1/T3))  where T3 = 16
    %   (V_fail4 - mu)/beta = -ln(-ln(1 - 1/T4))  where T4 = 55

    T_noTMD  = 16;   % year — historical return period (no TMD)
    T_TMD    = 55;   % year — historical return period (with TMD)
    V_f3     = params.V_quartering_fail_noTMD;   % 70 mph
    V_f4     = params.V_quartering_fail_TMD;     % 85 mph

    y_3 = -log(-log(1 - 1/T_noTMD));   % reduced variate for T=16
    y_4 = -log(-log(1 - 1/T_TMD));     % reduced variate for T=55

    beta_gumbel = (V_f4 - V_f3) / (y_4 - y_3);  % scale parameter
    mu_gumbel   = V_f3 - beta_gumbel * y_3;       % location parameter

    % Derived statistics for context
    V_50_derived = mu_gumbel - beta_gumbel * log(-log(1 - 1/50));
    V_100_derived = mu_gumbel - beta_gumbel * log(-log(1 - 1/100));

    fprintf('    Gumbel parameters: mu = %.2f mph, beta = %.2f mph\n',...
            mu_gumbel, beta_gumbel);
    fprintf('    Implied V_50 = %.1f mph, V_100 = %.1f mph\n',...
            V_50_derived, V_100_derived);

    %% ===== FAILURE THRESHOLDS =====
    % Critical wind speeds causing structural failure for each scenario

    % Scenario 1: Original design intent (welded connections, perp. wind only)
    %   Building was designed for V=100mph perpendicular. CJP groove welds
    %   develop full base metal strength — capacity ~2.5x the bolt group.
    %   Failure speed well above any realistic NYC wind.
    V_fail_1 = 120;

    % Scenario 2: As-built (bolted), but only perpendicular wind considered
    %   Bolt group was sized to D/C ~ 0.9 under perpendicular design wind.
    %   Failure occurs when V^2 exceeds capacity, i.e. V/V_design = 1/sqrt(0.9).
    %   With V_design ~ 100 mph: V_fail ~ 105 mph
    V_fail_2 = 105;

    % Scenario 3: As-built (bolted), quartering wind, NO TMD
    %   The actual crisis: bolted connections + 40% force increase from
    %   quartering = failure at ~70 mph. TMD is off (power failure).
    V_fail_3 = params.V_quartering_fail_noTMD;  % 70 mph

    % Scenario 4: As-built (bolted), quartering wind, WITH TMD
    %   TMD reduces dynamic response, raising failure threshold to ~85 mph
    V_fail_4 = params.V_quartering_fail_TMD;  % 85 mph

    V_fail = [V_fail_1, V_fail_2, V_fail_3, V_fail_4];
    scenario_names = {
        '1: Design intent (welded, perp.)'
        '2: As-built (bolted, perp. only)'
        '3: As-built (bolted, quart., NO TMD)'
        '4: As-built (bolted, quart., WITH TMD)'
    };
    scenario_colors = {[0.2 0.8 0.2], [0.2 0.6 1.0], [1 0.2 0.2], [1 0.7 0.1]};

    %% ===== MONTE CARLO SIMULATION =====
    N_samples = 100000;
    rng(42);  % reproducibility

    % Generate annual max wind speeds from Gumbel distribution
    % Inverse CDF: V = mu - beta * ln(-ln(U)), U ~ Uniform(0,1)
    U = rand(N_samples, 1);
    V_annual = mu_gumbel - beta_gumbel * log(-log(U));

    % Annual failure probability for each scenario
    P_fail = zeros(1, 4);
    for s = 1:4
        P_fail(s) = sum(V_annual >= V_fail(s)) / N_samples;
    end

    % Return periods
    T_return = 1 ./ P_fail;
    T_return(P_fail == 0) = Inf;

    % Reliability index: beta = -Phi_inv(Pf)
    % Using erfinv (base MATLAB) instead of norminv (Statistics Toolbox):
    %   norminv(p) = sqrt(2) * erfinv(2*p - 1)
    beta_rel = zeros(1,4);
    for s = 1:4
        if P_fail(s) > 0 && P_fail(s) < 1
            beta_rel(s) = -sqrt(2) * erfinv(2*P_fail(s) - 1);
        else
            beta_rel(s) = Inf;
        end
    end

    % Print results
    fprintf('    %-45s  Pf         T_return   beta\n', 'Scenario');
    fprintf('    %s\n', repmat('-', 1, 80));
    for s = 1:4
        fprintf('    %-45s  %.4e   %.1f yr   %.2f\n',...
                scenario_names{s}, P_fail(s), T_return(s), beta_rel(s));
    end

    % --- Analytical Gumbel probabilities (for validation) ---
    % P(V >= v) = 1 - F(v) = 1 - exp(-exp(-(v-mu)/beta))
    P_fail_exact = 1 - exp(-exp(-(V_fail - mu_gumbel)/beta_gumbel));
    fprintf('\n    Analytical (exact Gumbel) validation:\n');
    for s = 1:4
        fprintf('    Scenario %d: Pf_MC = %.4e, Pf_exact = %.4e\n',...
                s, P_fail(s), P_fail_exact(s));
    end

    %% ===== PLOTTING =====
    fig = figure('Name','Monte Carlo Reliability Analysis',...
                 'Position',[50 30 1500 950],'Color','k');

    % --- Panel 1: Wind speed PDF with failure thresholds ---
    ax1 = subplot(2,3,1); hold on;
    set(ax1,'Color','k','XColor','w','YColor','w');

    % Histogram of Monte Carlo samples
    v_lo = max(0, mu_gumbel - 3*beta_gumbel);
    v_hi = mu_gumbel + 8*beta_gumbel;
    edges = linspace(v_lo, v_hi, 80);
    histogram(V_annual, edges, 'FaceColor', [0.3 0.3 0.6], 'EdgeColor', 'none',...
              'FaceAlpha', 0.7, 'Normalization', 'pdf');

    % Overlay analytical Gumbel PDF
    v_plot = linspace(v_lo, v_hi, 500);
    z_gumbel = (v_plot - mu_gumbel) / beta_gumbel;
    pdf_gumbel = (1/beta_gumbel) .* exp(-(z_gumbel + exp(-z_gumbel)));
    plot(v_plot, pdf_gumbel, '-w', 'LineWidth', 2);

    % Failure threshold lines
    for s = 1:4
        xline(V_fail(s), '-', 'Color', scenario_colors{s}, 'LineWidth', 2,...
              'Label', sprintf('V_{fail}=%d', V_fail(s)),...
              'LabelColor', scenario_colors{s}, 'FontSize', 8,...
              'LabelHorizontalAlignment', 'left',...
              'LabelVerticalAlignment', 'top');
    end

    xlabel('Annual Max Wind Speed (mph)','Color','w');
    ylabel('PDF','Color','w');
    title('Annual Max Wind Speed Distribution','Color','w','FontSize',12);
    grid on; set(ax1,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 2: CDF with failure probabilities ---
    ax2 = subplot(2,3,2); hold on;
    set(ax2,'Color','k','XColor','w','YColor','w');

    cdf_gumbel = exp(-exp(-z_gumbel));
    plot(v_plot, 1 - cdf_gumbel, '-w', 'LineWidth', 2);  % exceedance probability

    for s = 1:4
        P_exceed = 1 - exp(-exp(-(V_fail(s) - mu_gumbel)/beta_gumbel));
        plot(V_fail(s), P_exceed, 'o', 'Color', scenario_colors{s},...
             'MarkerSize', 12, 'MarkerFaceColor', scenario_colors{s}, 'LineWidth', 2);
        text(V_fail(s)+2, P_exceed*1.3, sprintf('P_f=%.3f', P_exceed),...
             'Color', scenario_colors{s}, 'FontSize', 9);
    end

    set(ax2, 'YScale', 'log');
    xlabel('Wind Speed (mph)','Color','w');
    ylabel('P(V_{annual} \geq v)','Color','w');
    title('Exceedance Probability (CDF)','Color','w','FontSize',12);
    ylim([1e-4 1]);
    grid on; set(ax2,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 3: Return period bar chart ---
    ax3 = subplot(2,3,3); hold on;
    set(ax3,'Color','k','XColor','w','YColor','w');

    T_return_exact = 1 ./ P_fail_exact;
    b = bar(1:4, T_return_exact, 0.6);
    b.FaceColor = 'flat';
    for s = 1:4
        b.CData(s,:) = scenario_colors{s};
    end

    % Add value labels
    for s = 1:4
        text(s, T_return_exact(s)+5, sprintf('%.0f yr', T_return_exact(s)),...
             'Color','w','FontSize',11,'FontWeight','bold',...
             'HorizontalAlignment','center');
    end

    % Target reliability line
    yline(50, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5,...
          'Label', '50-yr target', 'LabelColor', [0.7 0.7 0.7]);

    set(ax3, 'XTick', 1:4, 'XTickLabel', {'Design','Bolted\newline(perp.)','Bolted\newline(quart.)','Bolted\newline+TMD'});
    set(ax3, 'XTickLabelRotation', 0);
    ylabel('Return Period (years)','Color','w');
    title('Failure Return Periods','Color','w','FontSize',12);
    grid on; set(ax3,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 4: Reliability index comparison ---
    ax4 = subplot(2,3,4); hold on;
    set(ax4,'Color','k','XColor','w','YColor','w');

    beta_exact = -sqrt(2) * erfinv(2*P_fail_exact - 1);
    b2 = bar(1:4, beta_exact, 0.6);
    b2.FaceColor = 'flat';
    for s = 1:4
        b2.CData(s,:) = scenario_colors{s};
    end

    for s = 1:4
        text(s, beta_exact(s)+0.1, sprintf('\\beta=%.2f', beta_exact(s)),...
             'Color','w','FontSize',11,'FontWeight','bold',...
             'HorizontalAlignment','center');
    end

    % ASCE target reliability lines
    yline(3.0, '--', 'Color', [0.5 1 0.5], 'LineWidth', 1.5,...
          'Label', '\beta_{target}=3.0 (ASCE)', 'LabelColor', [0.5 1 0.5]);
    yline(2.5, ':', 'Color', [1 1 0.5], 'LineWidth', 1,...
          'Label', '\beta=2.5 (marginal)', 'LabelColor', [1 1 0.5]);

    set(ax4, 'XTick', 1:4, 'XTickLabel', {'Design','Bolted\newline(perp.)','Bolted\newline(quart.)','Bolted\newline+TMD'});
    ylabel('Reliability Index \beta','Color','w');
    title('Reliability Index Comparison','Color','w','FontSize',12);
    grid on; set(ax4,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 5: Convergence of Monte Carlo estimate ---
    ax5 = subplot(2,3,5); hold on;
    set(ax5,'Color','k','XColor','w','YColor','w');

    % Show convergence for the critical scenario (Scenario 3)
    n_pts = 500;
    n_vec = round(logspace(2, log10(N_samples), n_pts));
    Pf_running = zeros(size(n_vec));
    failures_3 = V_annual >= V_fail_3;
    for i = 1:length(n_vec)
        Pf_running(i) = sum(failures_3(1:n_vec(i))) / n_vec(i);
    end

    plot(n_vec, Pf_running, '-', 'Color', scenario_colors{3}, 'LineWidth', 1.5);
    yline(P_fail_exact(3), '--w', 'LineWidth', 1.5);

    % 95% confidence interval
    ci_upper = P_fail_exact(3) + 1.96*sqrt(P_fail_exact(3)*(1-P_fail_exact(3))./n_vec);
    ci_lower = P_fail_exact(3) - 1.96*sqrt(P_fail_exact(3)*(1-P_fail_exact(3))./n_vec);
    fill([n_vec fliplr(n_vec)], [ci_upper fliplr(ci_lower)],...
         scenario_colors{3}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

    set(ax5, 'XScale', 'log');
    xlabel('Number of Samples','Color','w');
    ylabel('P_f estimate','Color','w');
    title('MC Convergence (Scenario 3)','Color','w','FontSize',12);
    legend({'MC estimate','Exact Gumbel','95% CI'},...
           'TextColor','w','Color',[0.15 0.15 0.15],'Location','northeast');
    grid on; set(ax5,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 6: Risk narrative / summary ---
    ax6 = subplot(2,3,6);
    set(ax6,'Color','k','XColor','k','YColor','k','XTick',[],'YTick',[]);
    box off; axis off;

    summary = {
        '\bf{Monte Carlo Reliability Summary}'
        ''
        sprintf('Samples: N = %d', N_samples)
        sprintf('Distribution: Gumbel Type I')
        sprintf('  \\mu = %.2f mph, \\beta = %.2f mph', mu_gumbel, beta_gumbel)
        ''
        '\bf{Annual Failure Probabilities:}'
        sprintf('  1. Design (welded):    P_f = %.2e  (T=%.0f yr)', P_fail_exact(1), T_return_exact(1))
        sprintf('  2. Bolted (perp.):     P_f = %.2e  (T=%.0f yr)', P_fail_exact(2), T_return_exact(2))
        sprintf('  3. Bolted+quart-TMD:   P_f = %.2e  (T=%.0f yr)', P_fail_exact(3), T_return_exact(3))
        sprintf('  4. Bolted+quart+TMD:   P_f = %.2e  (T=%.0f yr)', P_fail_exact(4), T_return_exact(4))
        ''
        '\bf{The Crisis:}'
        'Scenario 3 shows ~6% annual failure'
        'probability — a 1-in-16 year event.'
        'For a 30-year building lifetime:'
        sprintf('  P(fail in 30yr) = %.1f%%', (1-(1-P_fail_exact(3))^30)*100)
        ''
        '\bf{NIST 2021 Note:}'
        'Modern analysis suggests face winds'
        'may actually govern over quartering.'
    };
    text(0.05, 0.97, summary, 'Color', 'w', 'FontSize', 9.5,...
         'VerticalAlignment', 'top', 'FontName', 'FixedWidth',...
         'Interpreter', 'tex');

    sgtitle('Module 5: Monte Carlo Reliability Analysis — Citicorp Center',...
            'Color','w','FontSize',16,'FontWeight','bold');
end
