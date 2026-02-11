function fig = visualize_building_3d(params)
%VISUALIZE_BUILDING_3D  3D rendering of the Citicorp Center structural system
%  Shows the unique column placement, chevron bracing, and TMD location.

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

    % Floor plates (every 8 stories for clarity)
    story_h = (H - Hs) / params.stories;
    for i = 0:8:params.stories
        z = Hs + i * story_h;
        if z > H, break; end
        floor_x = [0 W W 0 0];
        floor_y = [0 0 W W 0];
        floor_z = z * ones(1,5);
        plot3(floor_x, floor_y, floor_z, '-', 'Color', [0.3 0.3 0.5], 'LineWidth', 0.5);
    end

    % Ground plane
    fill3([0 W W 0], [0 0 W W], [0 0 0 0], [0.15 0.15 0.15],...
          'EdgeColor',[0.3 0.3 0.3], 'FaceAlpha', 0.5);

    % Corner vertical edges of the tower
    corners = [0 0; W 0; W W; 0 W];
    for c = 1:4
        plot3([corners(c,1) corners(c,1)], [corners(c,2) corners(c,2)],...
              [Hs H], '-', 'Color', [0.4 0.4 0.6], 'LineWidth', 1);
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

    % ===== CHEVRON BRACING (inverted V pattern on each face) =====
    tier_height = 8 * story_h; % 8-story tiers
    brace_color = [0.2 0.8 1.0]; % Cyan

    % South face (y=0): chevrons in x-z plane
    draw_chevrons(0, W, 0, 0, Hs, H, tier_height, 'y', brace_color);
    % North face (y=W)
    draw_chevrons(0, W, W, W, Hs, H, tier_height, 'y', brace_color);
    % West face (x=0)
    draw_chevrons(0, 0, 0, W, Hs, H, tier_height, 'x', brace_color);
    % East face (x=W)
    draw_chevrons(W, W, 0, W, Hs, H, tier_height, 'x', brace_color);

    % ===== SLOPED ROOF (45 degree crown) =====
    roof_peak = H + W/4; % approximate peak height
    roof_color = [0.7 0.7 0.7];
    % Triangular roof on south/north faces
    plot3([0 W/2 W], [0 W/2 0], [H roof_peak H], '-', 'Color', roof_color, 'LineWidth', 2);
    plot3([0 W/2 W], [W W/2 W], [H roof_peak H], '-', 'Color', roof_color, 'LineWidth', 2);
    % Ridge line
    plot3([W/2 W/2], [0 W], [roof_peak roof_peak], '-', 'Color', roof_color, 'LineWidth', 2);
    % Roof edges
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
    stilt_labels = {'South','East','North','West'};
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

    sgtitle('Module 1: Citicorp Center 3D Structural System','Color','w','FontSize',16,'FontWeight','bold');
end

function draw_chevrons(x1, x2, y1, y2, z_bot, z_top, tier_h, fixed_axis, col)
%DRAW_CHEVRONS  Draw V-bracing (chevron) on one face
%  Citicorp chevrons: apex at top center of each tier,
%  legs going DOWN to bottom corners. Each tier is one V.
    z = z_bot;
    while z + tier_h <= z_top
        top_z = z + tier_h;

        if strcmp(fixed_axis, 'y')
            % Face is in x-z plane, y is fixed
            xm = (x1 + x2)/2;
            % V: from bottom corners up to top center
            plot3([x1 xm], [y1 y1], [z top_z], '-', 'Color', col, 'LineWidth', 2);
            plot3([x2 xm], [y2 y2], [z top_z], '-', 'Color', col, 'LineWidth', 2);
        else
            % Face is in y-z plane, x is fixed
            ym = (y1 + y2)/2;
            plot3([x1 x1], [y1 ym], [z top_z], '-', 'Color', col, 'LineWidth', 2);
            plot3([x2 x2], [y2 ym], [z top_z], '-', 'Color', col, 'LineWidth', 2);
        end

        z = z + tier_h;
    end
end
