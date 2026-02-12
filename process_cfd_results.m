function [fig, cfd_results] = process_cfd_results(params, wind_results)
%PROCESS_CFD_RESULTS  Read OpenFOAM results and compare with ASCE 7-22
%  Reads force coefficients and force history from postProcessing/.
%  Computes base shear, overturning moment, and Cd/Cl.
%  Compares CFD-derived loads with the ASCE 7-22 analytical results.
%
%  Prerequisites: OpenFOAM case must have been run (postProcessing/ exists).

    cfd_dir = fullfile(pwd, 'citicorp_cfd', 'postProcessing');

    %% ===== CHECK FOR CFD OUTPUT =====
    if ~exist(cfd_dir, 'dir')
        warning('CFD:NoResults', ...
            'No postProcessing/ directory found. Run the OpenFOAM case first.');
        cfd_results = struct();
        fig = figure('Visible', 'off');
        return;
    end

    %% ===== READ FORCE COEFFICIENTS =====
    % forceCoeffs output: time Cd Cl CmRoll CmPitch CmYaw Cd(f) Cd(r) Cl(f) Cl(r)
    coeff_dir = fullfile(cfd_dir, 'towerForceCoeffs');
    coeff_data = read_openfoam_timeseries(coeff_dir);

    %% ===== READ FORCES =====
    % forces output: time (fx fy fz)(pressure) (fx fy fz)(viscous) (fx fy fz)(porous)
    force_dir = fullfile(cfd_dir, 'forces');
    force_data = read_openfoam_forces(force_dir);

    %% ===== EXTRACT CONVERGED VALUES =====
    % Use last 20% of iterations for averaging (steady-state solution)
    n_total = size(coeff_data, 1);
    n_avg = max(1, round(0.2 * n_total));
    idx_avg = (n_total - n_avg + 1):n_total;

    Cd_mean = mean(coeff_data(idx_avg, 2));
    Cl_mean = mean(coeff_data(idx_avg, 3));
    Cd_std  = std(coeff_data(idx_avg, 2));

    fprintf('    CFD Force Coefficients (last %d iterations):\n', n_avg);
    fprintf('      Cd = %.4f +/- %.4f\n', Cd_mean, Cd_std);
    fprintf('      Cl = %.4f\n', Cl_mean);

    % Base shear from forces file (sum of pressure + viscous)
    n_f = size(force_data, 1);
    n_avg_f = max(1, round(0.2 * n_f));
    idx_f = (n_f - n_avg_f + 1):n_f;

    % Columns: time, fpx, fpy, fpz, fvx, fvy, fvz
    Fx_mean = mean(force_data(idx_f, 2) + force_data(idx_f, 5));  % pressure + viscous x
    Fy_mean = mean(force_data(idx_f, 3) + force_data(idx_f, 6));
    Fz_mean = mean(force_data(idx_f, 4) + force_data(idx_f, 7));

    % Convert from N to kips (1 kip = 4448.22 N)
    N_to_kips = 1 / 4448.22;
    Fx_kips = Fx_mean * N_to_kips;
    Fy_kips = Fy_mean * N_to_kips;

    fprintf('    CFD Base Shear (converged):\n');
    fprintf('      Fx (drag)  = %.1f kN = %.1f kips\n', Fx_mean/1000, Fx_kips);
    fprintf('      Fy (cross) = %.1f kN = %.1f kips\n', Fy_mean/1000, Fy_kips);

    %% ===== COMPARE WITH ASCE 7-22 =====
    V_asce7 = wind_results.V_base_perp;   % kips (analytical base shear)
    M_asce7 = wind_results.M_perp;        % kip-ft (overturning moment)

    % CFD overturning moment estimate (force * arm at ~0.6*H)
    H_m = params.height * 0.3048;  % ft to m
    arm_m = 0.6 * H_m;            % approximate moment arm
    M_cfd_Nm = Fx_mean * arm_m;
    M_cfd_kipft = M_cfd_Nm * 0.000737562;  % N-m to kip-ft

    ratio_V = Fx_kips / max(V_asce7, 1);
    ratio_M = M_cfd_kipft / max(M_asce7, 1);

    fprintf('\n    ASCE 7-22 vs CFD Comparison:\n');
    fprintf('      Base shear:  ASCE 7 = %.0f kips, CFD = %.0f kips (ratio = %.2f)\n',...
            V_asce7, Fx_kips, ratio_V);
    fprintf('      Overturn. M: ASCE 7 = %.0f kip-ft, CFD = %.0f kip-ft (ratio = %.2f)\n',...
            M_asce7, M_cfd_kipft, ratio_M);

    %% ===== PACK RESULTS =====
    cfd_results = struct();
    cfd_results.Cd = Cd_mean;
    cfd_results.Cl = Cl_mean;
    cfd_results.Cd_std = Cd_std;
    cfd_results.Fx_kips = Fx_kips;
    cfd_results.Fy_kips = Fy_kips;
    cfd_results.M_cfd_kipft = M_cfd_kipft;
    cfd_results.V_asce7_kips = V_asce7;
    cfd_results.M_asce7_kipft = M_asce7;
    cfd_results.ratio_V = ratio_V;
    cfd_results.ratio_M = ratio_M;
    cfd_results.coeff_history = coeff_data;
    cfd_results.force_history = force_data;

    %% ===== PLOT =====
    fig = figure('Name','CFD Post-Processing — Force Comparison',...
                 'Position',[60 60 1500 800],'Color','k');

    % --- Panel 1: Cd convergence history ---
    ax1 = subplot(2,2,1); hold on;
    set(ax1,'Color','k','XColor','w','YColor','w');
    plot(coeff_data(:,1), coeff_data(:,2), '-', 'Color', [0.2 0.7 1], 'LineWidth', 1);
    yline(Cd_mean, '--', 'Color', [1 0.5 0.2], 'LineWidth', 1.5,...
          'Label', sprintf('Cd = %.4f', Cd_mean), 'LabelColor', [1 0.5 0.2]);
    xlabel('Iteration','Color','w'); ylabel('Cd','Color','w');
    title('Drag Coefficient Convergence','Color','w','FontSize',11);
    grid on; set(ax1,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 2: Base shear comparison ---
    ax2 = subplot(2,2,2); hold on;
    set(ax2,'Color','k','XColor','w','YColor','w');
    b = bar([V_asce7, Fx_kips], 'FaceColor', 'flat');
    b.CData = [0.2 0.6 1; 1 0.3 0.3];
    set(ax2, 'XTickLabel', {'ASCE 7-22', 'CFD (RANS)'});
    ylabel('Base Shear (kips)','Color','w');
    title('Base Shear Comparison','Color','w','FontSize',11);
    text(2, Fx_kips + V_asce7 * 0.03, sprintf('%.0f%%', ratio_V * 100),...
         'Color','w','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
    grid on; set(ax2,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 3: Force history (Fx over iterations) ---
    ax3 = subplot(2,2,3); hold on;
    set(ax3,'Color','k','XColor','w','YColor','w');
    Fx_hist = (force_data(:,2) + force_data(:,5)) * N_to_kips;
    plot(force_data(:,1), Fx_hist, '-', 'Color', [0.8 0.3 0.3], 'LineWidth', 1);
    yline(Fx_kips, '--', 'Color', [1 0.8 0.2], 'LineWidth', 1.5);
    yline(V_asce7, ':', 'Color', [0.3 0.8 0.3], 'LineWidth', 1.5,...
          'Label', 'ASCE 7', 'LabelColor', [0.3 0.8 0.3]);
    xlabel('Iteration','Color','w'); ylabel('Fx (kips)','Color','w');
    title('Drag Force Convergence','Color','w','FontSize',11);
    legend({'CFD Fx', 'CFD mean', 'ASCE 7-22'},...
           'TextColor','w','Color',[0.15 0.15 0.15]);
    grid on; set(ax3,'GridColor',[0.3 0.3 0.3]);

    % --- Panel 4: Summary table ---
    ax4 = subplot(2,2,4);
    set(ax4,'Color','k','XColor','k','YColor','k');
    axis off;

    summary = {
        'PARAMETER',         'ASCE 7-22',                     'CFD (RANS)';
        'Base shear',        sprintf('%.0f kips', V_asce7),   sprintf('%.0f kips', Fx_kips);
        'Overturn. M',       sprintf('%.0f kip-ft', M_asce7), sprintf('%.0f kip-ft', M_cfd_kipft);
        'Drag coeff (Cd)',   '—',                             sprintf('%.4f', Cd_mean);
        'Cross-force (Fy)',  '—',                             sprintf('%.0f kips', Fy_kips);
        'Shear ratio',       '1.00',                          sprintf('%.2f', ratio_V);
    };

    y_pos = 0.85;
    for row = 1:size(summary, 1)
        if row == 1
            col = [0.9 0.9 0.9]; fw = 'bold'; fs = 11;
        else
            col = [0.7 0.7 0.7]; fw = 'normal'; fs = 10;
        end
        text(0.05, y_pos, summary{row, 1}, 'Color', col, 'FontSize', fs,...
             'FontWeight', fw, 'Units', 'normalized');
        text(0.45, y_pos, summary{row, 2}, 'Color', [0.3 0.7 1], 'FontSize', fs,...
             'FontWeight', fw, 'Units', 'normalized');
        text(0.75, y_pos, summary{row, 3}, 'Color', [1 0.4 0.4], 'FontSize', fs,...
             'FontWeight', fw, 'Units', 'normalized');
        y_pos = y_pos - 0.14;
    end

    sgtitle('Module 8: CFD Results — ASCE 7-22 vs OpenFOAM RANS',...
            'Color','w','FontSize',16,'FontWeight','bold');
end

%% ===== HELPER: READ OPENFOAM FORCE COEFFICIENTS =====
function data = read_openfoam_timeseries(dir_path)
%READ_OPENFOAM_TIMESERIES  Read OpenFOAM postProcessing time-series data
    % Find the time directory (usually '0')
    subdirs = dir(dir_path);
    subdirs = subdirs([subdirs.isdir] & ~startsWith({subdirs.name}, '.'));
    if isempty(subdirs)
        error('No time directories found in %s', dir_path);
    end
    time_dir = fullfile(dir_path, subdirs(1).name);

    % Find the .dat file
    dat_files = dir(fullfile(time_dir, '*.dat'));
    if isempty(dat_files)
        dat_files = dir(fullfile(time_dir, 'coefficient*'));
    end
    if isempty(dat_files)
        error('No data files found in %s', time_dir);
    end

    fpath = fullfile(time_dir, dat_files(1).name);
    fid = fopen(fpath, 'r');
    raw = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = raw{1};

    % Skip comment lines (start with #)
    data_lines = lines(~startsWith(lines, '#'));
    n = length(data_lines);
    if n == 0
        data = zeros(0, 10);
        return;
    end

    % Parse numeric data
    first = str2num(data_lines{1}); %#ok<ST2NM>
    ncols = length(first);
    data = zeros(n, ncols);
    data(1, :) = first;
    for i = 2:n
        vals = str2num(data_lines{i}); %#ok<ST2NM>
        if length(vals) == ncols
            data(i, :) = vals;
        end
    end
end

%% ===== HELPER: READ OPENFOAM FORCES =====
function data = read_openfoam_forces(dir_path)
%READ_OPENFOAM_FORCES  Read OpenFOAM forces output
%  Returns: [time, fpx, fpy, fpz, fvx, fvy, fvz]
    subdirs = dir(dir_path);
    subdirs = subdirs([subdirs.isdir] & ~startsWith({subdirs.name}, '.'));
    if isempty(subdirs)
        error('No time directories found in %s', dir_path);
    end
    time_dir = fullfile(dir_path, subdirs(1).name);

    % Look for force.dat or force_0.dat
    dat_files = dir(fullfile(time_dir, 'force*'));
    if isempty(dat_files)
        error('No force files found in %s', time_dir);
    end

    fpath = fullfile(time_dir, dat_files(1).name);
    fid = fopen(fpath, 'r');
    raw = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = raw{1};

    % Skip comments, parse parenthesized format:
    % time ((fpx fpy fpz) (fvx fvy fvz) (fpox fpoy fpoz))
    data_lines = lines(~startsWith(lines, '#'));
    n = length(data_lines);
    data = zeros(n, 7);  % time, fpx, fpy, fpz, fvx, fvy, fvz

    for i = 1:n
        line = data_lines{i};
        % Extract all numbers from the line
        nums = regexp(line, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match');
        vals = str2double(nums);
        if length(vals) >= 7
            data(i, :) = vals(1:7);  % time + pressure(3) + viscous(3)
        end
    end
end
