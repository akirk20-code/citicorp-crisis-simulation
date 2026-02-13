function fig = draw_structural_elevation(params)
%DRAW_STRUCTURAL_ELEVATION  2D structural elevation drawings of Citicorp Center
%  Shows front elevation (0 deg) and quartering view (45 deg) side-by-side.
%  Clean engineering drawing style — structural frame only.

    W  = params.width;          % 157 ft
    H  = params.height;         % 915 ft
    Hs = params.stilt_height;   % 114 ft (stilt height / tower bottom)
    Sc = params.stilt_size;     % 24 ft stilt cross-section
    sh = (H - Hs) / params.stories;  % story height
    n_tiers    = params.n_tiers;      % 6
    tier_h     = params.tier_stories * sh;  % 8-story tier height
    brace_top  = Hs + n_tiers * tier_h;
    roof_peak  = H + W/4;

    % Colors
    col_brace  = [1.0 0.55 0.0];   % Orange — chevron braces
    col_stilt  = [0.9 0.6 0.1];    % Gold — stilts
    col_edge   = [0.6 0.4 0.15];   % Dark orange — building edges
    col_roof   = [0.7 0.7 0.7];    % Grey — roof
    col_dim    = [0.5 0.5 0.5];    % Grey — dimension lines
    col_wind   = [0.2 0.8 0.2];    % Green — 0 deg wind
    col_wind45 = [1.0 0.4 0.4];    % Red — 45 deg wind

    fig = figure('Name','Structural Elevation — 0° and 45°',...
                 'Position',[50 50 1400 800],'Color','k');

    %% ===== LEFT PANEL: FRONT ELEVATION (0 degrees) =====
    ax1 = subplot(1,2,1); hold on;
    set(ax1,'Color','k','XColor','w','YColor','w');

    % --- Ground line ---
    plot([-30 W+30], [0 0], '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);

    % --- Stilt (single stilt visible at face midpoint) ---
    sx = W/2;
    hw = Sc/2;
    % Stilt box
    plot([sx-hw sx-hw], [0 Hs], '-', 'Color', col_stilt, 'LineWidth', 4);
    plot([sx+hw sx+hw], [0 Hs], '-', 'Color', col_stilt, 'LineWidth', 4);
    plot([sx-hw sx+hw], [0 0], '-', 'Color', col_stilt, 'LineWidth', 4);
    plot([sx-hw sx+hw], [Hs Hs], '-', 'Color', col_stilt, 'LineWidth', 4);

    % --- Building outline (corner edges) ---
    plot([0 0], [Hs H], '-', 'Color', col_edge, 'LineWidth', 2);
    plot([W W], [Hs H], '-', 'Color', col_edge, 'LineWidth', 2);
    % Horizontal at stilt top
    plot([0 W], [Hs Hs], '-', 'Color', col_edge, 'LineWidth', 1.5);

    % --- Chevron V-braces (6 tiers) ---
    z = Hs;
    for t = 1:n_tiers
        z_top = z + tier_h;
        xm = W/2;
        % Left leg: top-left corner down to bottom midpoint
        plot([0 xm], [z_top z], '-', 'Color', col_brace, 'LineWidth', 2.5);
        % Right leg: top-right corner down to bottom midpoint
        plot([W xm], [z_top z], '-', 'Color', col_brace, 'LineWidth', 2.5);
        % Horizontal chord at tier top
        plot([0 W], [z_top z_top], '-', 'Color', col_edge*0.8, 'LineWidth', 1);
        z = z_top;
    end

    % --- Roof (45-degree crown) ---
    plot([0 W/2], [H roof_peak], '-', 'Color', col_roof, 'LineWidth', 2);
    plot([W W/2], [H roof_peak], '-', 'Color', col_roof, 'LineWidth', 2);
    plot([0 W], [H H], '-', 'Color', col_edge, 'LineWidth', 1.5);

    % --- TMD marker ---
    tmd_z = Hs + 54 * sh;
    plot(W/2, tmd_z, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    text(W/2 + 12, tmd_z, 'TMD', 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

    % --- Wind arrows (perpendicular, 0 degrees) ---
    arrow_x = -25;
    for az = linspace(Hs+50, H-50, 5)
        quiver(arrow_x, az, 18, 0, 0, 'Color', col_wind, 'LineWidth', 2, ...
               'MaxHeadSize', 1.5);
    end
    text(arrow_x-5, H/2, '0° Wind', 'Color', col_wind, 'FontSize', 11, ...
         'Rotation', 90, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    % --- Dimension annotations ---
    % Height
    dim_x = W + 20;
    plot([dim_x dim_x], [0 H], '-', 'Color', col_dim, 'LineWidth', 0.5);
    plot([dim_x-3 dim_x+3], [0 0], '-', 'Color', col_dim, 'LineWidth', 0.5);
    plot([dim_x-3 dim_x+3], [H H], '-', 'Color', col_dim, 'LineWidth', 0.5);
    text(dim_x+5, H/2, sprintf('%.0f ft', H), 'Color', col_dim, ...
         'FontSize', 9, 'Rotation', 90, 'HorizontalAlignment', 'center');

    % Stilt height
    dim_x2 = W + 40;
    plot([dim_x2 dim_x2], [0 Hs], '-', 'Color', col_dim, 'LineWidth', 0.5);
    plot([dim_x2-3 dim_x2+3], [0 0], '-', 'Color', col_dim, 'LineWidth', 0.5);
    plot([dim_x2-3 dim_x2+3], [Hs Hs], '-', 'Color', col_dim, 'LineWidth', 0.5);
    text(dim_x2+5, Hs/2, sprintf('%.0f ft', Hs), 'Color', col_dim, ...
         'FontSize', 8, 'Rotation', 90, 'HorizontalAlignment', 'center');

    % Width
    dim_y = -25;
    plot([0 W], [dim_y dim_y], '-', 'Color', col_dim, 'LineWidth', 0.5);
    plot([0 0], [dim_y-3 dim_y+3], '-', 'Color', col_dim, 'LineWidth', 0.5);
    plot([W W], [dim_y-3 dim_y+3], '-', 'Color', col_dim, 'LineWidth', 0.5);
    text(W/2, dim_y-12, sprintf('%.0f ft', W), 'Color', col_dim, ...
         'FontSize', 9, 'HorizontalAlignment', 'center');

    % Tier labels
    z = Hs;
    for t = 1:n_tiers
        z_mid = z + tier_h/2;
        text(-15, z_mid, sprintf('T%d', t), 'Color', col_brace*0.7, ...
             'FontSize', 7, 'HorizontalAlignment', 'center');
        z = z + tier_h;
    end

    title('Front Elevation (0° — Perpendicular Wind)', 'Color', 'w', 'FontSize', 13);
    xlabel('Width (ft)', 'Color', 'w');
    ylabel('Height (ft)', 'Color', 'w');
    axis equal;
    xlim([-50 W+60]);
    ylim([-40 roof_peak+30]);

    %% ===== RIGHT PANEL: QUARTERING VIEW (45 degrees) =====
    ax2 = subplot(1,2,2); hold on;
    set(ax2,'Color','k','XColor','w','YColor','w');

    % At 45 degrees, the building is viewed diagonally (corner-on).
    % Projected width = W * sqrt(2) = 222 ft
    % We see two faces simultaneously — project both onto the view plane.
    Wp = W * sqrt(2);   % projected width at 45 deg
    half_Wp = Wp / 2;

    % Corner positions projected at 45 degrees (from center):
    %   Left corner (0,0):   projected x = 0
    %   South stilt (W/2,0): projected x = W/2 * cos(45) = W/(2*sqrt(2))
    %   East stilt (W,W/2):  projected x = W*cos(45) + W/2*sin(45)
    %                       = W/sqrt(2) * (1 + 0.5) = 3W/(2*sqrt(2))
    %   Far corner (W,W):    projected x = W*sqrt(2)
    % Center shifted so building is centered on the panel.

    % Simplified: project all 4 corners and 4 stilt midpoints
    % View direction: (1,1,0)/sqrt(2), screen-x perpendicular to that
    % Screen-x = (-1,1,0)/sqrt(2), so projected x = (-px + py)/sqrt(2)

    % Corner (0,0) -> proj_x = 0
    % Corner (W,0) -> proj_x = (-W)/sqrt(2) = -W/sqrt(2)
    % Corner (W,W) -> proj_x = 0
    % Corner (0,W) -> proj_x = W/sqrt(2)
    % So visible range: [-W/sqrt(2), W/sqrt(2)] = [-111, 111]

    proj = @(px,py) (-px + py) / sqrt(2);  % project to screen x

    % Corners
    c1 = proj(0, 0);      % near-left corner  (SW) = 0
    c2 = proj(W, 0);      % near-right corner (SE) = -W/sqrt(2)
    c3 = proj(W, W);      % far corner        (NE) = 0
    c4 = proj(0, W);      % far-left corner   (NW) = +W/sqrt(2)

    % Stilt midpoints
    s_south = proj(W/2, 0);     % = -W/(2*sqrt(2))
    s_east  = proj(W, W/2);     % = -W/(2*sqrt(2))  same!
    s_north = proj(W/2, W);     % = +W/(2*sqrt(2))
    s_west  = proj(0, W/2);     % = +W/(2*sqrt(2))  same!

    % Note: at 45° view, South stilt and East stilt project to SAME x position!
    % This is correct — they overlap from this viewing angle.

    % Visible face edges: left silhouette (SE corner) and right silhouette (NW corner)
    left_edge = c2;   % -W/sqrt(2) ≈ -111
    right_edge = c4;  % +W/sqrt(2) ≈ +111

    % --- Ground line ---
    plot([left_edge-30, right_edge+30], [0 0], '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);

    % --- Stilts (at 45°, we see 2 stilts: South+East overlap, North+West overlap) ---
    % Near stilts (South/East overlap)
    stilt_x1 = s_south;
    plot([stilt_x1-hw stilt_x1-hw], [0 Hs], '-', 'Color', col_stilt, 'LineWidth', 4);
    plot([stilt_x1+hw stilt_x1+hw], [0 Hs], '-', 'Color', col_stilt, 'LineWidth', 4);
    plot([stilt_x1-hw stilt_x1+hw], [0 0], '-', 'Color', col_stilt, 'LineWidth', 4);
    plot([stilt_x1-hw stilt_x1+hw], [Hs Hs], '-', 'Color', col_stilt, 'LineWidth', 4);

    % Far stilts (North/West overlap)
    stilt_x2 = s_north;
    plot([stilt_x2-hw stilt_x2-hw], [0 Hs], '-', 'Color', col_stilt, 'LineWidth', 4);
    plot([stilt_x2+hw stilt_x2+hw], [0 Hs], '-', 'Color', col_stilt, 'LineWidth', 4);
    plot([stilt_x2-hw stilt_x2+hw], [0 0], '-', 'Color', col_stilt, 'LineWidth', 4);
    plot([stilt_x2-hw stilt_x2+hw], [Hs Hs], '-', 'Color', col_stilt, 'LineWidth', 4);

    % --- Building outline (silhouette edges at 45 deg) ---
    plot([left_edge left_edge], [Hs H], '-', 'Color', col_edge, 'LineWidth', 2);
    plot([right_edge right_edge], [Hs H], '-', 'Color', col_edge, 'LineWidth', 2);
    % Near corner (SW/NE overlap at center)
    plot([c1 c1], [Hs H], ':', 'Color', col_edge*0.6, 'LineWidth', 1);
    % Horizontal at stilt top
    plot([left_edge right_edge], [Hs Hs], '-', 'Color', col_edge, 'LineWidth', 1.5);

    % --- Chevron V-braces on BOTH visible faces ---
    % South face (y=0): corners at (0,0) and (W,0), midpoint at (W/2,0)
    % West face (x=0): corners at (0,0) and (0,W), midpoint at (0,W/2)

    z = Hs;
    for t = 1:n_tiers
        z_top = z + tier_h;

        % South face V-brace (projected)
        % Top corners: (0,0)->c1=0, (W,0)->c2=-W/sqrt(2)
        % Bottom midpoint: (W/2,0) -> s_south = -W/(2*sqrt(2))
        plot([c1 s_south], [z_top z], '-', 'Color', col_brace, 'LineWidth', 2.5);
        plot([c2 s_south], [z_top z], '-', 'Color', col_brace, 'LineWidth', 2.5);

        % West face V-brace (projected)
        % Top corners: (0,0)->c1=0, (0,W)->c4=+W/sqrt(2)
        % Bottom midpoint: (0,W/2) -> s_west = +W/(2*sqrt(2))
        plot([c1 s_west], [z_top z], '-', 'Color', col_brace*0.85, 'LineWidth', 2);
        plot([c4 s_west], [z_top z], '-', 'Color', col_brace*0.85, 'LineWidth', 2);

        % Horizontal chords
        plot([left_edge right_edge], [z_top z_top], '-', 'Color', col_edge*0.6, 'LineWidth', 0.8);

        z = z_top;
    end

    % --- Roof at 45 degrees ---
    % Ridge runs along (W/2, y) for y=0 to W, projected:
    ridge_left  = proj(W/2, 0);   % = -W/(2*sqrt(2))
    ridge_right = proj(W/2, W);   % = +W/(2*sqrt(2))
    % Roof peak is above ridge at H + W/4
    % From SE corner to ridge peak
    plot([left_edge ridge_left], [H roof_peak], '-', 'Color', col_roof, 'LineWidth', 2);
    % From NW corner to ridge peak
    plot([right_edge ridge_right], [H roof_peak], '-', 'Color', col_roof, 'LineWidth', 2);
    % Ridge line
    plot([ridge_left ridge_right], [roof_peak roof_peak], '-', 'Color', col_roof, 'LineWidth', 2);
    % Near-corner roof edges
    plot([c1 (ridge_left+ridge_right)/2], [H roof_peak], ':', 'Color', col_roof*0.7, 'LineWidth', 1);
    % Top edge
    plot([left_edge right_edge], [H H], '-', 'Color', col_edge, 'LineWidth', 1.5);

    % --- TMD marker ---
    plot(0, tmd_z, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    text(8, tmd_z, 'TMD', 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

    % --- Wind arrows (quartering, 45 degrees — comes toward near corner) ---
    arrow_x = left_edge - 35;
    for az = linspace(Hs+50, H-50, 5)
        quiver(arrow_x, az, 18, 0, 0, 'Color', col_wind45, 'LineWidth', 2, ...
               'MaxHeadSize', 1.5);
    end
    text(arrow_x-5, H/2, '45° Wind', 'Color', col_wind45, 'FontSize', 11, ...
         'Rotation', 90, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    % --- Dimension annotations ---
    dim_x = right_edge + 20;
    plot([dim_x dim_x], [0 H], '-', 'Color', col_dim, 'LineWidth', 0.5);
    plot([dim_x-3 dim_x+3], [0 0], '-', 'Color', col_dim, 'LineWidth', 0.5);
    plot([dim_x-3 dim_x+3], [H H], '-', 'Color', col_dim, 'LineWidth', 0.5);
    text(dim_x+5, H/2, sprintf('%.0f ft', H), 'Color', col_dim, ...
         'FontSize', 9, 'Rotation', 90, 'HorizontalAlignment', 'center');

    % Projected width
    dim_y = -25;
    plot([left_edge right_edge], [dim_y dim_y], '-', 'Color', col_dim, 'LineWidth', 0.5);
    plot([left_edge left_edge], [dim_y-3 dim_y+3], '-', 'Color', col_dim, 'LineWidth', 0.5);
    plot([right_edge right_edge], [dim_y-3 dim_y+3], '-', 'Color', col_dim, 'LineWidth', 0.5);
    text(0, dim_y-12, sprintf('%.0f ft (diagonal)', Wp), 'Color', col_dim, ...
         'FontSize', 9, 'HorizontalAlignment', 'center');

    % Face labels on the 45-degree view
    text((c1+c2)/2, Hs-15, 'South Face', 'Color', col_brace*0.7, 'FontSize', 8, ...
         'HorizontalAlignment', 'center');
    text((c1+c4)/2, Hs-15, 'West Face', 'Color', col_brace*0.6, 'FontSize', 8, ...
         'HorizontalAlignment', 'center');

    title('Quartering View (45° — Corner-On Wind)', 'Color', 'w', 'FontSize', 13);
    xlabel('Projected Width (ft)', 'Color', 'w');
    ylabel('Height (ft)', 'Color', 'w');
    axis equal;
    xlim([left_edge-55 right_edge+55]);
    ylim([-40 roof_peak+30]);

    sgtitle('Citicorp Center — Structural Elevation Drawings',...
            'Color','w','FontSize',16,'FontWeight','bold');
end
