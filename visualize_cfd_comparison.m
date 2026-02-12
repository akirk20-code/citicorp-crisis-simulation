function fig = visualize_cfd_comparison(params, wind_results, cfd_results)
%VISUALIZE_CFD_COMPARISON  Side-by-side ASCE 7-22 vs CFD wind load comparison
%  Compares analytical wind pressure profile with CFD-derived forces.
%  Shows pressure distribution, base shear convergence, and summary.

    %% ===== ASCE 7-22 DATA =====
    z_ft = wind_results.z;                % height (ft)
    z_m  = z_ft * 0.3048;                 % convert to meters
    p_perp = wind_results.p_net_perp;     % net pressure (psf)
    p_Pa   = p_perp * 47.88;             % psf to Pa

    H_ft = params.height;
    Hs_ft = params.stilt_height;
    H_m = H_ft * 0.3048;
    Hs_m = Hs_ft * 0.3048;
    W_m = params.width * 0.3048;          % building width in meters

    % ASCE 7 Cp values
    Cp_ww = wind_results.Cp_windward;
    Cp_lw = wind_results.Cp_leeward;
    Cp_sw = wind_results.Cp_sidewall;

    %% ===== FIGURE =====
    fig = figure('Name','ASCE 7 vs CFD — Wind Load Comparison',...
                 'Position',[50 50 1600 900],'Color','k');

    % --- Panel 1: Pressure profile comparison ---
    ax1 = subplot(2,3,1); hold on;
    set(ax1,'Color','k','XColor','w','YColor','w');

    % ASCE 7 analytical profile
    plot(p_Pa, z_m, '-', 'Color', [0.3 0.7 1], 'LineWidth', 2);

    % Reference: dynamic pressure q at Uref
    rho = 1.225;
    Uref = 44.7;
    q_ref = 0.5 * rho * Uref^2;  % ~1223 Pa

    % Log-law velocity profile for CFD comparison
    z0 = 0.1; kappa = 0.41;
    ustar = Uref * kappa / log(10/z0);
    z_profile = linspace(1, H_m, 100);
    U_profile = (ustar / kappa) * log(z_profile / z0);
    q_profile = 0.5 * rho * U_profile.^2;
    % Net pressure ~ (Cp_front - Cp_back) * q  ≈ 1.3 * q (typical)
    p_cfd_est = 1.3 * q_profile;
    plot(p_cfd_est, z_profile, '--', 'Color', [1 0.4 0.4], 'LineWidth', 1.5);

    % Mark stilt height
    yline(Hs_m, ':', 'Color', [0.9 0.6 0.1], 'LineWidth', 1,...
          'Label', 'Stilt top', 'LabelColor', [0.9 0.6 0.1]);

    xlabel('Net Pressure (Pa)','Color','w');
    ylabel('Height (m)','Color','w');
    title('Wind Pressure Profile','Color','w','FontSize',11);
    legend({'ASCE 7-22', 'CFD (log-law est.)', 'Stilt top'},...
           'TextColor','w','Color',[0.15 0.15 0.15],'Location','southeast');
    grid on; set(ax1,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 2: Cp comparison ---
    ax2 = subplot(2,3,2); hold on;
    set(ax2,'Color','k','XColor','w','YColor','w');

    % ASCE 7 Cp values (code-specified)
    asce7_Cp = [Cp_ww, abs(Cp_lw), abs(Cp_sw)];
    face_labels = {'Windward', 'Leeward', 'Sidewall'};

    b = bar(1:3, asce7_Cp, 0.5, 'FaceColor', 'flat');
    b.CData = [0.3 0.7 1; 0.3 0.7 1; 0.3 0.7 1];

    % If CFD Cd is available, show it as reference
    if isfield(cfd_results, 'Cd') && cfd_results.Cd > 0
        yline(cfd_results.Cd, '--', 'Color', [1 0.4 0.4], 'LineWidth', 2,...
              'Label', sprintf('CFD Cd = %.3f', cfd_results.Cd),...
              'LabelColor', [1 0.4 0.4]);
    end

    set(ax2, 'XTick', 1:3, 'XTickLabel', face_labels);
    ylabel('|Cp|','Color','w');
    title('Pressure Coefficients','Color','w','FontSize',11);
    grid on; set(ax2,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 3: Base shear bar chart ---
    ax3 = subplot(2,3,3); hold on;
    set(ax3,'Color','k','XColor','w','YColor','w');

    V_asce7 = wind_results.V_base_perp;
    V_quar  = wind_results.V_base_quar;
    V_cfd   = cfd_results.Fx_kips;

    bar_data = [V_asce7; V_quar; V_cfd];
    b = bar(bar_data, 0.6, 'FaceColor', 'flat');
    b.CData = [0.3 0.7 1; 0.8 0.5 0.2; 1 0.3 0.3];
    set(ax3, 'XTick', 1:3, 'XTickLabel', {'ASCE 7 (0\circ)', 'ASCE 7 (45\circ)', 'CFD (0\circ)'});
    ylabel('Base Shear (kips)','Color','w');
    title('Base Shear Comparison','Color','w','FontSize',11);

    % Annotate values
    for i = 1:3
        text(i, bar_data(i) + max(bar_data)*0.03, sprintf('%.0f', bar_data(i)),...
             'Color','w','FontSize',10,'FontWeight','bold','HorizontalAlignment','center');
    end
    grid on; set(ax3,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 4: Cd convergence ---
    ax4 = subplot(2,3,4); hold on;
    set(ax4,'Color','k','XColor','w','YColor','w');

    if isfield(cfd_results, 'coeff_history') && ~isempty(cfd_results.coeff_history)
        ch = cfd_results.coeff_history;
        plot(ch(:,1), ch(:,2), '-', 'Color', [0.2 0.8 1], 'LineWidth', 1);
        plot(ch(:,1), ch(:,3), '-', 'Color', [1 0.4 0.4], 'LineWidth', 1);
        yline(cfd_results.Cd, '--', 'Color', [0.9 0.9 0.2], 'LineWidth', 1.5);
        legend({'Cd (drag)', 'Cl (lift)', 'Cd mean'},...
               'TextColor','w','Color',[0.15 0.15 0.15]);
    else
        text(0.5, 0.5, 'No CFD data', 'Color', 'w', 'FontSize', 14,...
             'HorizontalAlignment', 'center', 'Units', 'normalized');
    end
    xlabel('Iteration','Color','w');
    ylabel('Coefficient','Color','w');
    title('Force Coefficient History','Color','w','FontSize',11);
    grid on; set(ax4,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 5: Overturning moment ---
    ax5 = subplot(2,3,5); hold on;
    set(ax5,'Color','k','XColor','w','YColor','w');

    M_asce7 = wind_results.M_perp;
    M_quar  = wind_results.M_quar;
    M_cfd   = cfd_results.M_cfd_kipft;

    bar_m = [M_asce7; M_quar; M_cfd] / 1e6;  % to million kip-ft
    b = bar(bar_m, 0.6, 'FaceColor', 'flat');
    b.CData = [0.3 0.7 1; 0.8 0.5 0.2; 1 0.3 0.3];
    set(ax5, 'XTick', 1:3, 'XTickLabel', {'ASCE 7 (0\circ)', 'ASCE 7 (45\circ)', 'CFD (0\circ)'});
    ylabel('Overturning Moment (10^6 kip-ft)','Color','w');
    title('Overturning Moment Comparison','Color','w','FontSize',11);

    for i = 1:3
        text(i, bar_m(i) + max(bar_m)*0.03, sprintf('%.2f M', bar_m(i)),...
             'Color','w','FontSize',10,'FontWeight','bold','HorizontalAlignment','center');
    end
    grid on; set(ax5,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 6: Summary text ---
    ax6 = subplot(2,3,6);
    set(ax6,'Color','k','XColor','k','YColor','k');
    axis off;

    lines = {
        'ASCE 7-22 vs OpenFOAM RANS',                                    [1 1 1],     14, 'bold';
        '',                                                                [1 1 1],     10, 'normal';
        sprintf('Shear ratio (CFD/ASCE7): %.2f', cfd_results.ratio_V),   [0.3 1 0.3], 12, 'bold';
        sprintf('Moment ratio: %.2f', cfd_results.ratio_M),              [0.3 1 0.3], 12, 'bold';
        '',                                                                [1 1 1],     10, 'normal';
        'Wind: ABL inlet, Uref = 44.7 m/s',                              [0.6 0.6 0.6], 10, 'normal';
        'Turbulence: k-omega SST',                                        [0.6 0.6 0.6], 10, 'normal';
        sprintf('Mesh: ~2-3M cells, L/H ~ %.0f', 720/H_m),              [0.6 0.6 0.6], 10, 'normal';
        '',                                                                [1 1 1],     10, 'normal';
        'If ratio ~ 0.8-1.2: ASCE 7 is',                                 [0.8 0.8 0.2], 10, 'normal';
        'a reasonable approximation.',                                     [0.8 0.8 0.2], 10, 'normal';
    };

    y = 0.92;
    for i = 1:size(lines, 1)
        if ~isempty(lines{i,1})
            text(0.05, y, lines{i,1}, 'Color', lines{i,2}, 'FontSize', lines{i,3},...
                 'FontWeight', lines{i,4}, 'Units', 'normalized');
        end
        y = y - 0.085;
    end

    sgtitle('Module 9: ASCE 7-22 vs CFD Wind Load Validation',...
            'Color','w','FontSize',16,'FontWeight','bold');
end
