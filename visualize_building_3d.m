function fig = visualize_building_3d(params)
%VISUALIZE_BUILDING_3D  3D rendering of the Citicorp Center structural system
%  Shows the unique column placement, chevron bracing, and TMD location.
%  Clean structural frame view — orange bracing, gold stilts, no floor plates.

    W = params.width;       % 157 ft square plan
    H = params.height;      % 915 ft total
    Hs = params.stilt_height; % 114 ft stilts
    Sc = params.stilt_size;   % 24 ft stilt cross-section
    Ca = params.cantilever;   % 72 ft cantilever

    fig = figure('Name','Citicorp Center — 3D Structural System',...
                 'Position',[50 50 1200 800],'Color','k');

    % --- Left panel: full 3D building ---
    ax1 = subplot(1,2,1); hold on;
    set(ax1,'Color','k','XColor','w','YColor','w','ZColor','w');

    story_h = (H - Hs) / params.stories;

    % Ground plane (subtle)
    fill3([0 W W 0], [0 0 W W], [0 0 0 0], [0.12 0.12 0.12],...
          'EdgeColor',[0.25 0.25 0.25], 'FaceAlpha', 0.3);

    % Corner vertical edges of the tower (subtle structural outline)
    edge_color = [0.6 0.4 0.15];  % Dark orange
    corners = [0 0; W 0; W W; 0 W];
    for c = 1:4
        plot3([corners(c,1) corners(c,1)], [corners(c,2) corners(c,2)],...
              [Hs H], '-', 'Color', edge_color, 'LineWidth', 1.5);
    end

    % Horizontal edges at stilt top and roof
    for z_lvl = [Hs, H]
        plot3([0 W W 0 0], [0 0 W W 0], z_lvl*ones(1,5), ...
              '-', 'Color', edge_color, 'LineWidth', 1.5);
    end

    % ===== STILTS (at midpoints of each face, NOT corners) =====
    stilt_positions = [W/2, 0;    % South face midpoint
                       W,   W/2;  % East face midpoint
                       W/2, W;    % North face midpoint
                       0,   W/2]; % West face midpoint
    stilt_colors = [0.9 0.6 0.1]; % Gold
    for s = 1:4
        cx = stilt_positions(s,1);
        cy = stilt_positions(s,2);
        hw = Sc/2;
        % Draw stilt as thick vertical lines (4 edges of box)
        sx = [cx-hw cx+hw cx+hw cx-hw cx-hw];
        sy = [cy-hw cy-hw cy+hw cy+hw cy-hw];
        % Bottom
        plot3(sx, sy, zeros(1,5), '-', 'Color', stilt_colors, 'LineWidth', 3);
        % Top
        plot3(sx, sy, Hs*ones(1,5), '-', 'Color', stilt_colors, 'LineWidth', 3);
        % Verticals
        for e = 1:4
            plot3([sx(e) sx(e)], [sy(e) sy(e)], [0 Hs],...
                  '-', 'Color', stilt_colors, 'LineWidth', 3);
        end
    end

    % ===== TRANSFER ZONE (base stories above stilts) =====
    if isfield(params, 'transfer_stories') && params.transfer_stories > 0
        n_transfer = params.transfer_stories;
    else
        n_transfer = 0;
    end
    n_tiers = params.n_tiers;          % 6
    tier_height = params.tier_stories * story_h;  % 8-story tiers
    brace_base = Hs + n_transfer * story_h;  % V-bracing starts above transfer
    brace_top = brace_base + n_tiers * tier_height;
    brace_color = [1.0 0.55 0.0]; % Orange
    xfer_color  = [1.0 0.4 0.8];  % Magenta

    if n_transfer > 0
        % Transfer V: stilt midpoint at Hs spreads to corners at brace_base
        % Horizontal at top of transfer zone
        plot3([0 W W 0 0], [0 0 W W 0], brace_base*ones(1,5), ...
              '-', 'Color', edge_color, 'LineWidth', 1.5);
        % South face: stilt center up to corners
        draw_transfer_v(0, W, 0, 0, Hs, brace_base, 'y', xfer_color);
        % North face
        draw_transfer_v(0, W, W, W, Hs, brace_base, 'y', xfer_color);
        % West face
        draw_transfer_v(0, 0, 0, W, Hs, brace_base, 'x', xfer_color);
        % East face
        draw_transfer_v(W, W, 0, W, Hs, brace_base, 'x', xfer_color);
    end

    % ===== CHEVRON BRACING (V pattern on each face) =====
    % "48 braces, in 6 tiers of 8" — V-braces start above transfer zone

    % South face (y=0): chevrons in x-z plane
    draw_chevrons(0, W, 0, 0, brace_base, brace_top, tier_height, 'y', brace_color);
    % North face (y=W)
    draw_chevrons(0, W, W, W, brace_base, brace_top, tier_height, 'y', brace_color);
    % West face (x=0)
    draw_chevrons(0, 0, 0, W, brace_base, brace_top, tier_height, 'x', brace_color);
    % East face (x=W)
    draw_chevrons(W, W, 0, W, brace_base, brace_top, tier_height, 'x', brace_color);

    % ===== TIER BOUNDARY RINGS (closed perimeter at each tier top) =====
    % Floor diaphragms connect chevron bracing on all 4 faces at corners,
    % forming closed rectangular rings — critical for quartering wind transfer.
    ring_color = brace_color * 0.6;  % Subdued orange
    z_tier = brace_base;
    while z_tier + tier_height <= brace_top + 0.1
        tz = min(z_tier + tier_height, brace_top);
        plot3([0 W W 0 0], [0 0 W W 0], tz*ones(1,5), ...
              '-', 'Color', ring_color, 'LineWidth', 1.3);
        z_tier = z_tier + tier_height;
    end

    % ===== SLOPED ROOF (45 degree crown) =====
    roof_peak = H + W/4; % approximate peak height
    roof_color = [0.7 0.7 0.7];
    % Triangular roof on south/north faces
    plot3([0 W/2 W], [0 W/2 0], [H roof_peak H], '-', 'Color', roof_color, 'LineWidth', 2);
    plot3([0 W/2 W], [W W/2 W], [H roof_peak H], '-', 'Color', roof_color, 'LineWidth', 2);
    % Ridge line
    plot3([W/2 W/2], [0 W], [roof_peak roof_peak], '-', 'Color', roof_color, 'LineWidth', 2);
    % Roof side edges
    plot3([0 0], [0 W], [H H], '-', 'Color', roof_color, 'LineWidth', 1);
    plot3([W W], [0 W], [H H], '-', 'Color', roof_color, 'LineWidth', 1);

    % ===== TMD MARKER (floor 63) =====
    tmd_z = Hs + 54 * story_h; % approximate floor 63
    scatter3(W/2, W/2, tmd_z, 200, 'r', 'filled', 'MarkerEdgeColor', 'w');
    text(W/2 + 15, W/2 + 15, tmd_z, 'TMD (400 tons)',...
         'Color','r','FontSize',10,'FontWeight','bold');

    % ===== WIND ARROWS =====
    % Perpendicular wind (0 degrees — hits south face)
    arrow_y = -60;
    for az = [40 W/2 W-40]
        quiver3(az, arrow_y, H/2, 0, 50, 0, 0,...
                'Color','g','LineWidth',2,'MaxHeadSize',2);
    end
    text(W/2, arrow_y-20, H/2, '0° Wind','Color','g',...
         'FontSize',11,'HorizontalAlignment','center','FontWeight','bold');

    % Quartering wind (45 degrees)
    arrow_start = -40;
    for az = [30 W/2 W-30]
        quiver3(arrow_start, arrow_start + az - W/2, H*0.7,...
                35, 35, 0, 0,...
                'Color',[1 0.4 0.4],'LineWidth',2,'MaxHeadSize',2);
    end
    text(arrow_start-25, arrow_start, H*0.7, '45° Wind',...
         'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');

    % Labels
    title('Citicorp Center — Structural System','Color','w','FontSize',14);
    xlabel('X (ft)','Color','w'); ylabel('Y (ft)','Color','w'); zlabel('Z (ft)','Color','w');
    view([-37 25]);
    axis equal; grid on;
    set(ax1,'GridColor',[0.3 0.3 0.3]);
    zlim([-20 roof_peak+50]);

    % --- Right panel: annotated cross-section ---
    ax2 = subplot(1,2,2); hold on;
    set(ax2,'Color','k','XColor','w','YColor','w');

    % Floor plan with stilt positions
    rectangle('Position',[0 0 W W],'EdgeColor','w','LineWidth',2);

    % Stilt positions (midpoints)
    for s = 1:4
        cx = stilt_positions(s,1);
        cy = stilt_positions(s,2);
        hw = Sc/2;
        rectangle('Position',[cx-hw cy-hw Sc Sc],...
                  'EdgeColor',stilt_colors,'FaceColor',stilt_colors,...
                  'LineWidth',2);
    end

    % Show where corners are NOT supported
    corner_pts = [0 0; W 0; W W; 0 W];
    for c = 1:4
        plot(corner_pts(c,1), corner_pts(c,2), 'rx', 'MarkerSize', 15, 'LineWidth', 3);
    end

    % Cantilever annotations
    annotation_color = [1 0.7 0.3];
    plot([0 W/2-Sc/2], [0 0], '--', 'Color', annotation_color, 'LineWidth', 1.5);
    text(20, -12, sprintf('%.0f ft cantilever', Ca), 'Color', annotation_color, 'FontSize', 9);

    % Church location (reason for column placement)
    rectangle('Position',[-30 -30 50 50],'EdgeColor',[0.5 0.5 0.5],...
              'LineStyle','--','LineWidth',1);
    text(-25, -45, 'St. Peter''s Church', 'Color', [0.5 0.5 0.5], 'FontSize', 8);

    % Wind direction indicators
    quiver(W/2, -40, 0, 30, 0, 'Color', 'g', 'LineWidth', 2, 'MaxHeadSize', 1);
    text(W/2+5, -35, '0°', 'Color', 'g', 'FontSize', 10);

    quiver(-30, -30, 25, 25, 0, 'Color', [1 0.4 0.4], 'LineWidth', 2, 'MaxHeadSize', 1);
    text(-45, -35, '45°', 'Color', [1 0.4 0.4], 'FontSize', 10);

    title('Plan View — Column Placement','Color','w','FontSize',14);
    xlabel('X (ft)','Color','w'); ylabel('Y (ft)','Color','w');
    axis equal;
    xlim([-60 W+30]); ylim([-60 W+30]);

    % Legend
    text(W+10, W-20, 'Columns at midpoints', 'Color', stilt_colors, 'FontSize', 9);
    text(W+10, W-40, 'No corner support', 'Color', 'r', 'FontSize', 9);
    text(W+10, W-60, '72 ft cantilevers', 'Color', annotation_color, 'FontSize', 9);
    text(W+10, W-80, 'Chevron V-bracing', 'Color', brace_color, 'FontSize', 9);
    text(W+10, W-100, 'Tier boundary rings', 'Color', brace_color*0.6, 'FontSize', 9);
    if n_transfer > 0
        text(W+10, W-120, 'Transfer zone', 'Color', xfer_color, 'FontSize', 9);
    end

    sgtitle('Module 1: Citicorp Center 3D Structural System','Color','w','FontSize',16,'FontWeight','bold');
end

function draw_transfer_v(x1, x2, y1, y2, z_bot, z_top, fixed_axis, col)
%DRAW_TRANSFER_V  Draw transfer V in the base zone (stilt center to corners)
%  Inverted V: apex at bottom-center (stilt midpoint), legs up to corners.
    if strcmp(fixed_axis, 'y')
        xm = (x1 + x2)/2;
        % From stilt center (bottom) up to left corner (top)
        plot3([xm x1], [y1 y1], [z_bot z_top], '-', 'Color', col, 'LineWidth', 2.5);
        % From stilt center (bottom) up to right corner (top)
        plot3([xm x2], [y2 y2], [z_bot z_top], '-', 'Color', col, 'LineWidth', 2.5);
    else
        ym = (y1 + y2)/2;
        plot3([x1 x1], [ym y1], [z_bot z_top], '-', 'Color', col, 'LineWidth', 2.5);
        plot3([x2 x2], [ym y2], [z_bot z_top], '-', 'Color', col, 'LineWidth', 2.5);
    end
end

function draw_chevrons(x1, x2, y1, y2, z_bot, z_top, tier_h, fixed_axis, col)
%DRAW_CHEVRONS  Draw V-bracing (chevron) on one face
%  Citicorp chevrons: apex at bottom center of each tier (V shape),
%  legs going from top corners DOWN to bottom midpoint.
    z = z_bot;
    while z + tier_h <= z_top + 0.1
        top_z = min(z + tier_h, z_top);

        if strcmp(fixed_axis, 'y')
            % Face is in x-z plane, y is fixed
            xm = (x1 + x2)/2;
            % V: from top corners down to bottom center
            plot3([x1 xm], [y1 y1], [top_z z], '-', 'Color', col, 'LineWidth', 2.5);
            plot3([x2 xm], [y2 y2], [top_z z], '-', 'Color', col, 'LineWidth', 2.5);
        else
            % Face is in y-z plane, x is fixed
            ym = (y1 + y2)/2;
            plot3([x1 x1], [y1 ym], [top_z z], '-', 'Color', col, 'LineWidth', 2.5);
            plot3([x2 x2], [y2 ym], [top_z z], '-', 'Color', col, 'LineWidth', 2.5);
        end
        % Note: horizontal chords drawn as full perimeter rings in main function

        z = z + tier_h;
    end
end
