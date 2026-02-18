function [fig, wind_results] = analyze_wind_pressure(params)
%ANALYZE_WIND_PRESSURE  Full ASCE 7-22 wind pressure analysis
%  Implements Chapters 26-28 methodology for the Citicorp Center.
%  Compares perpendicular (0 deg) vs quartering (45 deg) wind loading.
%
%  Returns wind_results struct with height-varying pressures for FEA module.

    H  = params.height;        % 915 ft
    W  = params.width;         % 157 ft
    Hs = params.stilt_height;  % 114 ft

    %% ===== ASCE 7-22 VELOCITY PRESSURE PROFILE =====
    % qz = 0.00256 * Kz * Kzt * Kd * Ke * V^2  (psf, V in mph)
    % Ref: ASCE 7-22 Eq. 26.10-1

    z = linspace(15, H, 200)';  % height array (ft), min 15 ft per code

    % --- Kz: Velocity Pressure Exposure Coefficient (Table 26.10-1) ---
    % Exposure B (urban): alpha = 7.0, zg = 1200 ft
    alpha_exp = 7.0;
    zg        = 1200;  % ft — gradient height for Exposure B
    Kz = zeros(size(z));
    for i = 1:length(z)
        if z(i) <= 15
            Kz(i) = 2.01 * (15/zg)^(2/alpha_exp);
        else
            Kz(i) = 2.01 * (z(i)/zg)^(2/alpha_exp);
        end
    end

    % Velocity pressure at each height (psf)
    V = params.V_design;  % 100 mph (1970 NYC code)
    qz = 0.00256 * Kz * params.Kzt * params.Kd * params.Ke * V^2;

    % Also compute for ASCE 7-22 wind speed
    V_modern = params.V_asce7_catIII;  % 120 mph
    qz_modern = 0.00256 * Kz * params.Kzt * params.Kd * params.Ke * V_modern^2;

    %% ===== GUST EFFECT FACTOR (Flexible Building) =====
    % ASCE 7-22 Section 26.11.5
    % For buildings with T1 > 1 sec, use Gf (flexible)

    T1   = params.T1;          % 6.5 sec
    n1   = 1/T1;               % natural frequency (Hz)
    beta_damp = params.damping_bare;  % 1% damping
    z_bar = 0.6 * H;           % equivalent height (Eq. 26.11-3)
    if z_bar < 15, z_bar = 15; end

    % Turbulence intensity at z_bar (Eq. 26.11-7)
    c_exp = 0.30;  % Exposure B constant
    Iz_bar = c_exp * (33/z_bar)^(1/6);

    % Background response factor Q (Eq. 26.11-8)
    Lz_bar = 320 * (z_bar/33)^(1/3);  % integral length scale (ft)
    B_sq = 1 / (1 + 0.63 * ((W + H)/Lz_bar)^0.63);
    Q = sqrt(B_sq);

    % Resonant response factor R (Eq. 26.11-10 through 26.11-15)
    V_zbar = 88 * (z_bar/33)^(1/alpha_exp);  % mean hourly speed at z_bar (ft/s)
    % using V_design converted: 100 mph * 88/60 for ft/s, then profile
    V_zbar = (V * 88/60) * (z_bar/zg)^(1/alpha_exp);

    N1 = n1 * Lz_bar / V_zbar;  % reduced frequency

    Rn = 7.47 * N1 / (1 + 10.3 * N1)^(5/3);  % Eq. 26.11-11

    % Aerodynamic admittance functions
    eta_h = 4.6 * n1 * H / V_zbar;
    eta_B = 4.6 * n1 * W / V_zbar;
    Rh = (1/eta_h) - (1/(2*eta_h^2)) * (1 - exp(-2*eta_h));
    RB = (1/eta_B) - (1/(2*eta_B^2)) * (1 - exp(-2*eta_B));

    R_sq = (1/beta_damp) * Rn * Rh * RB;
    R = sqrt(R_sq);

    % Peak factor (Eq. 26.11-9)
    gR = sqrt(2 * log(3600 * n1)) + 0.577 / sqrt(2 * log(3600 * n1));
    gQ = 3.4;
    gv = 3.4;

    % Gust effect factor for flexible building (Eq. 26.11-6)
    Gf = 0.925 * (1 + 1.7 * Iz_bar * sqrt(gQ^2 * B_sq + gR^2 * R_sq)) / ...
         (1 + 1.7 * gv * Iz_bar);

    fprintf('    Gust effect factor Gf = %.3f (flexible building)\n', Gf);
    fprintf('    Background factor Q = %.3f, Resonant factor R = %.3f\n', Q, R);

    %% ===== EXTERNAL PRESSURE COEFFICIENTS =====
    % ASCE 7-22 Figure 27.3-1 for enclosed buildings
    % L/B = 1.0 (square plan), so Cp depends on surface orientation

    % Windward wall:  Cp = +0.8
    % Leeward wall:   Cp = -0.5 (for L/B = 1.0)
    % Side walls:     Cp = -0.7
    % Roof windward:  Cp = -0.7 to -0.18 (depending on theta)

    Cp_windward = 0.8;
    Cp_leeward  = -0.5;
    Cp_sidewall = -0.7;

    %% ===== DESIGN WIND PRESSURES =====
    % p = qz * Gf * Cp - qi * (GCpi)   (Eq. 27.3-1)
    % For MWFRS, internal pressure cancels for overall forces
    % Net pressure on windward + leeward = qz*Gf*Cp_w + qh*Gf*|Cp_l|

    qh = qz(end);  % velocity pressure at roof height

    % Windward face pressures (vary with height)
    p_windward = qz * Gf * Cp_windward;

    % Leeward face pressure (constant = qh based)
    p_leeward = qh * Gf * abs(Cp_leeward) * ones(size(z));

    % Total net pressure on MWFRS (windward + leeward)
    p_net_perp = p_windward + p_leeward;

    %% ===== QUARTERING WIND ANALYSIS =====
    % At 45 degrees, wind pressure acts on two faces simultaneously.
    % ASCE 7-22 Section 27.3.6 Case 3: apply 75% of full design pressure
    % on each face (no velocity decomposition — the 0.75 already accounts
    % for oblique incidence).
    %
    % The key insight: for the Citicorp chevron system, quartering winds
    % create asymmetric loading — half the braces on a face carry increased
    % load while the perpendicular face braces redistribute forces

    % ASCE 7-22 Sec. 27.3.6 Case 3: 75% of full design pressure on each face
    % No velocity decomposition — the 0.75 factor already accounts for
    % oblique incidence. Using cos²(45°) would double-reduce.
    p_face1_45 = 0.75 * qz * Gf * Cp_windward;
    p_face2_45 = 0.75 * qz * Gf * Cp_windward;

    % Net on each windward face for quartering
    p_leeward_45 = 0.75 * qh * Gf * abs(Cp_leeward) * ones(size(z));
    p_net_face1 = p_face1_45 + p_leeward_45;
    p_net_face2 = p_face2_45 + p_leeward_45;

    % Total base shear comparison
    % Perpendicular: acts on one face only (width W)
    V_base_perp = trapz(z, p_net_perp * W);  % lbs
    % Quartering: acts on two faces (each width W)
    V_base_quar = trapz(z, p_net_face1 * W) * sqrt(2);  % resultant

    fprintf('    Base shear (perpendicular): %.0f kips\n', V_base_perp/1000);
    fprintf('    Base shear (quartering):    %.0f kips\n', V_base_quar/1000);
    fprintf('    Ratio (quartering/perp):    %.2f\n', V_base_quar/V_base_perp);

    %% ===== OVERTURNING MOMENTS =====
    M_perp = trapz(z, p_net_perp * W .* z);  % ft-lbs
    M_quar = trapz(z, p_net_face1 * W .* z) * sqrt(2);

    %% ===== PACK RESULTS FOR FEA =====
    wind_results = struct();
    wind_results.z = z;
    wind_results.qz = qz;
    wind_results.Gf = Gf;
    wind_results.p_net_perp = p_net_perp;
    wind_results.p_net_face1_45 = p_net_face1;
    wind_results.p_net_face2_45 = p_net_face2;
    wind_results.V_base_perp = V_base_perp;
    wind_results.V_base_quar = V_base_quar;
    wind_results.M_perp = M_perp;
    wind_results.M_quar = M_quar;
    wind_results.Cp_windward = Cp_windward;
    wind_results.Cp_leeward = Cp_leeward;
    wind_results.Cp_sidewall = Cp_sidewall;
    wind_results.p_windward = p_windward;
    wind_results.p_leeward = p_leeward;
    wind_results.qz_modern = qz_modern;
    wind_results.B_sq = B_sq;
    wind_results.R_sq = R_sq;
    wind_results.gQ = gQ;
    wind_results.gR = gR;
    wind_results.gv = gv;
    wind_results.Iz_bar = Iz_bar;
    wind_results.beta_damp = beta_damp;

    %% ===== PLOTTING =====
    fig = figure('Name','Wind Pressure Analysis (ASCE 7-22)',...
                 'Position',[100 50 1400 900],'Color','k');

    % --- Panel 1: Velocity pressure profile ---
    ax1 = subplot(2,3,1); hold on;
    set(ax1,'Color','k','XColor','w','YColor','w');
    plot(qz, z, '-c', 'LineWidth', 2);
    plot(qz_modern, z, '--', 'Color', [1 0.5 0], 'LineWidth', 2);
    yline(Hs, ':', 'Color', [0.9 0.6 0.1], 'LineWidth', 1.5, 'Label', 'Stilt Top',...
          'LabelHorizontalAlignment', 'left', 'FontSize', 9, 'LabelColor', [0.9 0.6 0.1]);
    xlabel('q_z (psf)','Color','w');
    ylabel('Height z (ft)','Color','w');
    title('Velocity Pressure Profile','Color','w','FontSize',12);
    legend({'V=100 mph (1970 code)','V=120 mph (ASCE 7-22)'},...
           'TextColor','w','Color',[0.15 0.15 0.15],'Location','southeast');
    grid on; set(ax1,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 2: Kz profile ---
    ax2 = subplot(2,3,2); hold on;
    set(ax2,'Color','k','XColor','w','YColor','w');
    plot(Kz, z, '-g', 'LineWidth', 2);
    xlabel('K_z','Color','w');
    ylabel('Height z (ft)','Color','w');
    title('Exposure Coefficient K_z (Exp. B)','Color','w','FontSize',12);
    text(Kz(end)*0.6, H*0.3, sprintf('\\alpha = %.1f\nz_g = %d ft\nExposure B\n(Urban/Suburban)',...
         alpha_exp, zg), 'Color','g','FontSize',10);
    grid on; set(ax2,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 3: Net pressure on building faces ---
    ax3 = subplot(2,3,3); hold on;
    set(ax3,'Color','k','XColor','w','YColor','w');
    plot(p_windward, z, '-r', 'LineWidth', 2);
    plot(-p_leeward, z, '-b', 'LineWidth', 2);
    plot(p_net_perp, z, '-w', 'LineWidth', 2.5);
    xlabel('Pressure (psf)','Color','w');
    ylabel('Height z (ft)','Color','w');
    title('Face Pressures (0° Wind)','Color','w','FontSize',12);
    legend({'Windward (+0.8)','Leeward (-0.5)','Net MWFRS'},...
           'TextColor','w','Color',[0.15 0.15 0.15],'Location','southeast');
    grid on; set(ax3,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 4: Perpendicular vs quartering net pressure ---
    ax4 = subplot(2,3,4); hold on;
    set(ax4,'Color','k','XColor','w','YColor','w');
    area(p_net_perp, z, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.4, 'EdgeColor', 'c', 'LineWidth', 1.5);
    area(p_net_face1, z, 'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.4, 'EdgeColor', [1 0.4 0.4], 'LineWidth', 1.5);
    xlabel('Net Pressure per Face (psf)','Color','w');
    ylabel('Height z (ft)','Color','w');
    title('Perpendicular vs Quartering (per face)','Color','w','FontSize',12);
    legend({'0° (one face)','45° (each of two faces)'},...
           'TextColor','w','Color',[0.15 0.15 0.15],'Location','southeast');
    grid on; set(ax4,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 5: Polar pressure distribution (using polaraxes) ---
    ax5 = subplot(2,3,5);
    delete(ax5);  % remove Cartesian axes, replace with polar
    ax5 = polaraxes('Position', [0.37 0.07 0.26 0.38]);
    hold(ax5, 'on');

    theta_angles = linspace(0, 2*pi, 361);
    % Build pressure coefficient around the building perimeter
    Cp_polar = zeros(size(theta_angles));
    for i = 1:length(theta_angles)
        th_deg = mod(rad2deg(theta_angles(i)), 360);
        if th_deg >= 315 || th_deg < 45
            Cp_polar(i) = Cp_windward;       % South (windward)
        elseif th_deg >= 45 && th_deg < 135
            Cp_polar(i) = Cp_sidewall;        % East (side)
        elseif th_deg >= 135 && th_deg < 225
            Cp_polar(i) = Cp_leeward;         % North (leeward)
        else
            Cp_polar(i) = Cp_sidewall;        % West (side)
        end
    end
    polarplot(ax5, theta_angles, abs(Cp_polar), '-c', 'LineWidth', 2);

    % Quartering wind rotates the pattern 45 degrees
    Cp_polar_45 = zeros(size(theta_angles));
    for i = 1:length(theta_angles)
        th_deg = mod(rad2deg(theta_angles(i)) - 45, 360);
        if th_deg >= 315 || th_deg < 45
            Cp_polar_45(i) = Cp_windward;
        elseif th_deg >= 45 && th_deg < 135
            Cp_polar_45(i) = Cp_sidewall;
        elseif th_deg >= 135 && th_deg < 225
            Cp_polar_45(i) = Cp_leeward;
        else
            Cp_polar_45(i) = Cp_sidewall;
        end
    end
    polarplot(ax5, theta_angles, abs(Cp_polar_45) * 0.75, '--', 'Color', [1 0.4 0.4], 'LineWidth', 2);
    title(ax5, '|C_p| Distribution (Polar View)','Color','w','FontSize',12);
    set(ax5, 'Color', 'k', 'ThetaColor', 'w', 'RColor', 'w');
    legend(ax5, {'0° wind','45° wind (\times0.75)'},...
           'TextColor','w','Color',[0.15 0.15 0.15],'Location','southoutside');

    % --- Panel 6: Summary text ---
    ax6 = subplot(2,3,6);
    set(ax6,'Color','k','XColor','k','YColor','k','XTick',[],'YTick',[]);
    box off; axis off;
    summary_text = {
        '\bf{ASCE 7-22 Wind Analysis Summary}'
        ''
        sprintf('Design Wind Speed: %d mph (1970 code)', params.V_design)
        sprintf('Modern ASCE 7-22: %d mph (Risk Cat III)', params.V_asce7_catIII)
        sprintf('Exposure Category: %s (Urban)', params.exposure)
        sprintf('Directionality K_d: %.2f', params.Kd)
        ''
        sprintf('Gust Effect Factor G_f: %.3f', Gf)
        sprintf('  Background Q = %.3f', Q)
        sprintf('  Resonant   R = %.3f', R)
        sprintf('  Period T_1 = %.1f sec', T1)
        ''
        sprintf('\\bf{Base Shear (V = %d mph):}', V)
        sprintf('  Perpendicular: %.0f kips', V_base_perp/1000)
        sprintf('  Quartering:    %.0f kips', V_base_quar/1000)
        sprintf('  Ratio: %.2f', V_base_quar/V_base_perp)
        ''
        sprintf('\\bf{Overturning Moment:}')
        sprintf('  Perpendicular: %.2e ft-lbs', M_perp)
        sprintf('  Quartering:    %.2e ft-lbs', M_quar)
    };
    text(0.05, 0.95, summary_text, 'Color', 'w', 'FontSize', 10,...
         'VerticalAlignment', 'top', 'FontName', 'FixedWidth',...
         'Interpreter', 'tex');

    sgtitle('Module 2: ASCE 7-22 Wind Pressure Analysis — Citicorp Center',...
            'Color','w','FontSize',16,'FontWeight','bold');
end
