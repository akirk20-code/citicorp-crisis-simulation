function fig = visualize_cfd_geometry()
%VISUALIZE_CFD_GEOMETRY  3D rendering of the OpenFOAM CFD domain
%  Shows the Citicorp tower, stilts, 18 surrounding buildings,
%  domain boundaries, refinement zones, and ABL wind inlet.
%  All dimensions in meters (matching generate_stl.py).

    %% ===== GEOMETRY (from generate_stl.py) =====
    TOWER_W = 47.85;       % 157 ft — square plan
    TOWER_H = 278.9;       % 915 ft — roof height
    STILT_H = 34.75;       % 114 ft — stilt height
    STILT_W = 7.32;        % 24 ft — stilt cross-section
    HT = TOWER_W / 2;     % half-tower
    HS = STILT_W / 2;     % half-stilt

    % Stilt centers (face midpoints)
    stilts = [0, -HT; HT, 0; 0, HT; -HT, 0];

    % Surrounding buildings: [cx, cy, wx, wy, height]
    surr = {
        '399 Park Ave',      -200,  -50,  60, 45, 162;
        '280 Park Ave',      -200,   50,  55, 40, 135;
        '780 Third Ave',      200,  -25,  50, 50, 175;
        '919 Third Ave',      200,   60,  55, 45, 182;
        '885 Third Ave',      200,  180,  35, 35, 138;
        '599 Lexington',     -100,   80,  45, 40, 199;
        '731 Lexington',     -100, -100,  50, 45, 180;
        '135 E 54th',           0,  100,  40, 35, 150;
        '153 E 53rd',          80,  -80,  35, 30, 120;
        'Park Ave Tower S',  -400,  -50,  55, 40, 150;
        'Park Ave Tower N',  -400,   60,  50, 35, 120;
        'Second Ave S',       400,  -60,  50, 40, 130;
        'Second Ave N',       400,   70,  45, 35, 110;
        'E 55th & Lex',      -110,  170,  40, 30, 100;
        'E 52nd & Lex',      -110, -180,  45, 35,  90;
        'E 55th & Third',     110,  170,  35, 30,  85;
        'E 52nd & Third',     110, -170,  40, 35,  95;
        'E 53rd mid',         -50, -100,  30, 25,  70;
    };

    % Domain bounds (from blockMeshDict)
    xMin = -200;  xMax = 520;
    yMin = -360;  yMax = 360;
    zMax = 600;

    %% ===== FIGURE =====
    fig = figure('Name','CFD Domain — Citicorp Center',...
                 'Position',[40 40 1600 900],'Color','k');

    % --- Left: Full domain overview ---
    ax1 = subplot(1,2,1); hold on;
    set(ax1,'Color','k','XColor','w','YColor','w','ZColor','w');

    % Domain wireframe
    draw_wireframe_box(ax1, xMin, yMin, 0, xMax, yMax, zMax, [0.25 0.25 0.25], '--', 0.5);

    % Ground plane
    fill3([xMin xMax xMax xMin], [yMin yMin yMax yMax], [0 0 0 0],...
          [0.12 0.12 0.12], 'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.4);

    % Surrounding buildings (gray boxes)
    for i = 1:size(surr, 1)
        cx = surr{i,2}; cy = surr{i,3};
        wx = surr{i,4}; wy = surr{i,5}; h = surr{i,6};
        draw_solid_box(ax1, cx-wx/2, cy-wy/2, 0, cx+wx/2, cy+wy/2, h,...
                       [0.35 0.35 0.40], 0.6);
    end

    % Citicorp Tower (cyan)
    tower_col = [0.2 0.8 1.0];
    draw_solid_box(ax1, -HT, -HT, STILT_H, HT, HT, TOWER_H, tower_col, 0.7);

    % Stilts (gold)
    stilt_col = [0.9 0.6 0.1];
    for s = 1:4
        cx = stilts(s,1); cy = stilts(s,2);
        draw_solid_box(ax1, cx-HS, cy-HS, 0, cx+HS, cy+HS, STILT_H, stilt_col, 0.9);
    end

    % Refinement zones (dashed outlines)
    % Building region
    draw_wireframe_box(ax1, -60, -60, 0, 60, 60, 320, [0.3 0.8 0.3], ':', 1);
    text(62, 0, 160, 'Refine L2', 'Color', [0.3 0.8 0.3], 'FontSize', 8);

    % Wake region
    draw_wireframe_box(ax1, 24, -100, 0, 300, 100, 200, [0.8 0.5 0.3], ':', 1);
    text(160, 102, 100, 'Wake L1', 'Color', [0.8 0.5 0.3], 'FontSize', 8);

    % Stilt underpass region
    draw_wireframe_box(ax1, -30, -30, 0, 30, 30, 40, [1 0.8 0.2], ':', 1);
    text(32, 0, 20, 'Stilt L3', 'Color', [1 0.8 0.2], 'FontSize', 8);

    % Wind arrows (ABL inlet at x = xMin)
    arrow_col = [0.3 1.0 0.3];
    for yy = [-150 0 150]
        for zz = [50 150 300]
            % Arrow length proportional to log-law velocity
            z0 = 0.1; Uref = 44.7; Zref = 10;
            ustar = Uref * 0.41 / log(Zref/z0);
            U_z = (ustar / 0.41) * log(max(zz,1)/z0);
            len = 30 + 60 * (U_z / 60);  % scale
            quiver3(xMin, yy, zz, len, 0, 0, 0,...
                    'Color', arrow_col, 'LineWidth', 1.2, 'MaxHeadSize', 0.6);
        end
    end
    text(xMin - 15, 0, zMax * 0.85, 'ABL Inlet', 'Color', arrow_col,...
         'FontSize', 11, 'FontWeight', 'bold', 'Rotation', 90,...
         'HorizontalAlignment', 'center');

    % Outlet label
    text(xMax + 10, 0, zMax * 0.5, 'Outlet', 'Color', [0.7 0.3 0.3],...
         'FontSize', 10, 'HorizontalAlignment', 'left');

    % Labels
    title('CFD Domain — Full Overview (720m x 720m x 600m)',...
          'Color', 'w', 'FontSize', 13);
    xlabel('X [m]', 'Color', 'w');
    ylabel('Y [m]', 'Color', 'w');
    zlabel('Z [m]', 'Color', 'w');
    view([-42 22]);
    axis equal; grid on;
    set(ax1, 'GridColor', [0.25 0.25 0.25]);

    % --- Right: Close-up of building cluster ---
    ax2 = subplot(1,2,2); hold on;
    set(ax2,'Color','k','XColor','w','YColor','w','ZColor','w');

    % Ground
    fill3([-250 250 250 -250], [-220 -220 220 220], [0 0 0 0],...
          [0.12 0.12 0.12], 'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.4);

    % Surrounding buildings (labeled)
    for i = 1:size(surr, 1)
        cx = surr{i,2}; cy = surr{i,3};
        wx = surr{i,4}; wy = surr{i,5}; h = surr{i,6};
        if abs(cx) <= 250 && abs(cy) <= 220
            draw_solid_box(ax2, cx-wx/2, cy-wy/2, 0, cx+wx/2, cy+wy/2, h,...
                           [0.35 0.35 0.40], 0.6);
            name = surr{i,1};
            text(cx, cy, h + 8, name, 'Color', [0.6 0.6 0.6],...
                 'FontSize', 6, 'HorizontalAlignment', 'center');
        end
    end

    % Tower
    draw_solid_box(ax2, -HT, -HT, STILT_H, HT, HT, TOWER_H, tower_col, 0.7);
    text(0, 0, TOWER_H + 12, 'CITICORP CENTER', 'Color', tower_col,...
         'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % Stilts
    for s = 1:4
        cx = stilts(s,1); cy = stilts(s,2);
        draw_solid_box(ax2, cx-HS, cy-HS, 0, cx+HS, cy+HS, STILT_H, stilt_col, 0.9);
    end

    % Wind arrow
    quiver3(-200, 0, 150, 80, 0, 0, 0,...
            'Color', arrow_col, 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
    text(-210, 0, 160, 'Wind', 'Color', arrow_col, 'FontSize', 11,...
         'FontWeight', 'bold', 'HorizontalAlignment', 'right');

    % Compass rose
    text(220, -200, 0, 'N', 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
    text(220,  200, 0, 'S', 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');

    title('Building Cluster — Close-Up with Labels',...
          'Color', 'w', 'FontSize', 13);
    xlabel('X [m]', 'Color', 'w');
    ylabel('Y [m]', 'Color', 'w');
    zlabel('Z [m]', 'Color', 'w');
    view([-35 30]);
    axis equal; grid on;
    set(ax2, 'GridColor', [0.25 0.25 0.25]);
    xlim([-250 250]); ylim([-220 220]); zlim([0 320]);

    sgtitle('Module 7: CFD Domain & Geometry Validation',...
            'Color', 'w', 'FontSize', 16, 'FontWeight', 'bold');
end

%% ===== HELPER FUNCTIONS =====

function draw_solid_box(ax, x0, y0, z0, x1, y1, z1, col, alpha)
%DRAW_SOLID_BOX  Draw a filled 3D box
    vx = [x0 x1 x1 x0 x0 x1 x1 x0];
    vy = [y0 y0 y1 y1 y0 y0 y1 y1];
    vz = [z0 z0 z0 z0 z1 z1 z1 z1];

    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6];
    patch(ax, 'Vertices', [vx' vy' vz'], 'Faces', faces,...
          'FaceColor', col, 'FaceAlpha', alpha,...
          'EdgeColor', col * 0.7, 'LineWidth', 0.5);
end

function draw_wireframe_box(~, x0, y0, z0, x1, y1, z1, col, style, lw)
%DRAW_WIREFRAME_BOX  Draw wireframe edges of a box
    % Bottom rectangle
    plot3([x0 x1 x1 x0 x0], [y0 y0 y1 y1 y0], [z0 z0 z0 z0 z0],...
          style, 'Color', col, 'LineWidth', lw);
    % Top rectangle
    plot3([x0 x1 x1 x0 x0], [y0 y0 y1 y1 y0], [z1 z1 z1 z1 z1],...
          style, 'Color', col, 'LineWidth', lw);
    % Verticals
    for i = 1:4
        xx = [x0 x1 x1 x0]; yy = [y0 y0 y1 y1];
        plot3([xx(i) xx(i)], [yy(i) yy(i)], [z0 z1],...
              style, 'Color', col, 'LineWidth', lw);
    end
end
