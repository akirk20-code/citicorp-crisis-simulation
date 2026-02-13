function interactive_structure()
%INTERACTIVE_STRUCTURE  Node-based structural editor for Citicorp Center
%  Drag purple nodes to reshape the structure — all connected beams follow.
%  Like Desmos Geometry / GeoGebra, but for structural elevations.
%
%  Controls:
%    Drag node     — Move node (all connected beams follow)
%    Click beam    — Select it (white highlight)
%    Delete key    — Remove selected node (+ connected beams) or beam
%    Escape        — Cancel add mode
%    [Add Node]    — Click drawing to place a new node
%    [Add Beam]    — Click two nodes to connect them
%    [Generate]    — Regenerate default Citicorp geometry from parameters
%
%  Usage: interactive_structure()

    %% === BUILDING PARAMETERS ===
    p = struct('width', 157, 'height', 915, 'stilt_height', 114, ...
               'stilt_size', 24, 'stories', 59, 'n_tiers', 6, ...
               'tier_stories', 8, 'transfer_stories', 2);

    %% === GENERATE INITIAL GEOMETRY ===
    [nodes, beams] = generate_geometry(p);
    % nodes: Nx2 double [x, y]
    % beams: Mx4 cell {from_node_idx, to_node_idx, type_str, color_str}

    %% === EDITOR STATE ===
    sel_node = 0;
    sel_beam = 0;
    dragging = false;
    drag_idx = 0;
    mode = 'select';       % 'select', 'add_node', 'beam_from', 'beam_to'
    beam_from_node = 0;
    snap_val = 5;          % snap to 5 ft grid (0 = off)
    show_ids = true;
    show_floors = true;
    grab_r = 15;           % grab radius in ft

    %% === REFERENCE IMAGE ===
    has_img = false; ref_img = [];
    img_on = false; img_x = 0; img_y = 0; img_sc = 1.0; img_op = 0.3;
    img_path = fullfile(fileparts(mfilename('fullpath')), '..', 'Picture1.png');
    if ~exist(img_path, 'file')
        img_path = 'C:\Users\kirka\OneDrive\Documents\GMU PhD\OR 750 Reliability, Safety, and Risk\Picture1.png';
    end
    if exist(img_path, 'file')
        ref_img = imread(img_path);
        has_img = true;
        img_on = true;
    end

    %% === CREATE FIGURE ===
    bg = [0.06 0.06 0.06];
    fig = figure('Name','Citicorp — Node-Based Structural Editor', ...
        'Position',[30 30 1600 920], 'Color', bg, ...
        'NumberTitle','off', 'MenuBar','none', 'ToolBar','figure', ...
        'WindowButtonDownFcn', @mouse_down, ...
        'WindowButtonMotionFcn', @mouse_move, ...
        'WindowButtonUpFcn', @mouse_up, ...
        'KeyPressFcn', @key_press);

    % Drawing axes — left 60%
    ax = axes('Parent',fig,'Position',[0.01 0.07 0.57 0.89]);
    set(ax,'Color',[0.02 0.02 0.02],'XColor','w','YColor','w');
    hold(ax,'on'); grid(ax,'on');
    set(ax,'GridColor',[0.15 0.15 0.15],'GridAlpha',0.6);

    % Coordinate readout (below drawing)
    coord_txt = uicontrol(fig,'Style','text', ...
        'Units','normalized','Position',[0.01 0.005 0.57 0.05], ...
        'BackgroundColor',bg,'ForegroundColor',[0.3 1 0.3], ...
        'FontSize',11,'FontName','Consolas','HorizontalAlignment','center');

    %% === RIGHT PANEL ===
    rpx = 0.60;  rpw = 0.39;

    % ---- QUICK GENERATE ----
    mk_hdr(rpx, 0.965, 'QUICK GENERATE', [1 .55 0]);
    gy = 0.940;
    mk_lbl(rpx,gy,'Stilt'); e_stilt  = mk_ed(rpx+0.035,gy,0.03, num2str(p.stilt_height));
    mk_lbl(rpx+0.075,gy,'Width'); e_width  = mk_ed(rpx+0.11,gy,0.03, num2str(p.width));
    mk_lbl(rpx+0.15,gy,'Height'); e_height = mk_ed(rpx+0.19,gy,0.035, num2str(p.height));
    gy2 = 0.910;
    mk_lbl(rpx,gy2,'Tiers'); e_tiers = mk_ed(rpx+0.035,gy2,0.025, num2str(p.n_tiers));
    mk_lbl(rpx+0.065,gy2,'St/Tier'); e_stpt  = mk_ed(rpx+0.11,gy2,0.025, num2str(p.tier_stories));
    mk_lbl(rpx+0.14,gy2,'Xfer'); e_xfer  = mk_ed(rpx+0.17,gy2,0.025, num2str(p.transfer_stories));
    uicontrol(fig,'Style','pushbutton','String','GENERATE', ...
        'Units','normalized','Position',[rpx+0.21 gy2 0.065 0.028], ...
        'BackgroundColor',[.1 .5 .2],'ForegroundColor','w', ...
        'FontSize',9,'FontWeight','bold','Callback',@cb_generate);

    % ---- TOGGLES + IMAGE ----
    mk_hdr(rpx, 0.88, 'DISPLAY & IMAGE', [.5 .7 1]);
    ty = 0.855;
    chk_snap = uicontrol(fig,'Style','checkbox','String','Snap','Value',1, ...
        'Units','normalized','Position',[rpx ty 0.035 0.022], ...
        'BackgroundColor',bg,'ForegroundColor',[.8 .8 .8],'FontSize',8, ...
        'Callback',@(~,~) read_toggles());
    e_snap = mk_ed(rpx+0.035, ty, 0.02, '5');
    chk_ids = uicontrol(fig,'Style','checkbox','String','Node IDs','Value',1, ...
        'Units','normalized','Position',[rpx+0.06 ty 0.055 0.022], ...
        'BackgroundColor',bg,'ForegroundColor',[.8 .8 .8],'FontSize',8, ...
        'Callback',@(~,~) read_toggles_redraw());
    chk_flr = uicontrol(fig,'Style','checkbox','String','Floors','Value',1, ...
        'Units','normalized','Position',[rpx+0.115 ty 0.045 0.022], ...
        'BackgroundColor',bg,'ForegroundColor',[.8 .8 .8],'FontSize',8, ...
        'Callback',@(~,~) read_toggles_redraw());
    chk_img = uicontrol(fig,'Style','checkbox','String','Image','Value',has_img, ...
        'Units','normalized','Position',[rpx+0.16 ty 0.045 0.022], ...
        'BackgroundColor',bg,'ForegroundColor',[.8 .8 .8],'FontSize',8, ...
        'Callback',@(~,~) read_toggles_redraw());

    iy = 0.83;
    mk_lbl(rpx,iy,'ImgX'); e_imgx = mk_ed(rpx+0.03,iy,0.025,'0');
    mk_lbl(rpx+0.06,iy,'ImgY'); e_imgy = mk_ed(rpx+0.09,iy,0.025,'0');
    mk_lbl(rpx+0.12,iy,'Scale'); e_imgsc = mk_ed(rpx+0.155,iy,0.025,'1.0');
    mk_lbl(rpx+0.185,iy,'Opac'); e_imgop = mk_ed(rpx+0.22,iy,0.025,'0.3');
    uicontrol(fig,'Style','pushbutton','String','Load...', ...
        'Units','normalized','Position',[rpx+0.25 iy 0.04 0.022], ...
        'BackgroundColor',[.3 .3 .4],'ForegroundColor','w','FontSize',8, ...
        'Callback',@cb_load_img);

    % ---- ACTION BUTTONS ----
    ay = 0.795;
    bw = rpw / 5;
    btn_labels = {'Select','Add Node','Add Beam','Delete','Export'};
    btn_colors = {[.25 .25 .35],[.1 .4 .15],[.1 .35 .45],[.45 .12 .1],[.15 .25 .5]};
    btn_cbs = {@cb_sel, @cb_addnode, @cb_addbeam, @cb_delete, @cb_export};
    for bi = 1:5
        uicontrol(fig,'Style','pushbutton','String',btn_labels{bi}, ...
            'Units','normalized','Position',[rpx+(bi-1)*bw ay bw-0.003 0.028], ...
            'BackgroundColor',btn_colors{bi},'ForegroundColor','w', ...
            'FontSize',8,'Callback',btn_cbs{bi});
    end

    % ---- NODE TABLE (left half) ----
    mk_hdr(rpx, 0.77, 'NODES', [1 .55 0]);
    node_tbl = uitable(fig,'Units','normalized', ...
        'Position',[rpx 0.52 rpw*0.45 0.245], ...
        'ColumnName',{'#','X (ft)','Y (ft)'}, ...
        'ColumnWidth',{25, 55, 55}, ...
        'ColumnFormat',{'numeric','numeric','numeric'}, ...
        'ColumnEditable',[false true true], ...
        'BackgroundColor',[.12 .12 .12; .16 .16 .16], ...
        'ForegroundColor',[.9 .9 .9],'FontSize',9, ...
        'CellEditCallback',@cb_node_edit, ...
        'CellSelectionCallback',@cb_node_sel);

    % ---- BEAM TABLE (right half) ----
    mk_hdr(rpx+rpw*0.47, 0.77, 'BEAMS', [1 .55 0]);
    beam_tbl = uitable(fig,'Units','normalized', ...
        'Position',[rpx+rpw*0.47 0.52 rpw*0.53 0.245], ...
        'ColumnName',{'#','From','To','Type','Color'}, ...
        'ColumnWidth',{22, 30, 30, 52, 45}, ...
        'ColumnFormat',{'numeric','numeric','numeric','char','char'}, ...
        'ColumnEditable',[false true true true true], ...
        'BackgroundColor',[.12 .12 .12; .16 .16 .16], ...
        'ForegroundColor',[.9 .9 .9],'FontSize',9, ...
        'CellEditCallback',@cb_beam_edit, ...
        'CellSelectionCallback',@cb_beam_sel);

    % ---- FLOOR SCHEDULE ----
    mk_hdr(rpx, 0.495, 'FLOOR SCHEDULE  (Left & Right XY per floor)', [.3 .85 .3]);
    floor_tbl = uitable(fig,'Units','normalized', ...
        'Position',[rpx 0.06 rpw 0.43], ...
        'ColumnName',{'Floor','Left X','Left Y','Right X','Right Y','Ht (ft)'}, ...
        'ColumnWidth',{38, 48, 52, 48, 52, 52}, ...
        'ColumnEditable',false(1,6), ...
        'BackgroundColor',[.09 .09 .09; .12 .12 .12], ...
        'ForegroundColor',[.6 .6 .6],'FontSize',8);

    % Status bar
    status_txt = uicontrol(fig,'Style','text','Units','normalized', ...
        'Position',[rpx 0.01 rpw 0.04], ...
        'BackgroundColor',bg,'ForegroundColor',[.5 .5 .5], ...
        'FontSize',9,'FontName','Consolas','HorizontalAlignment','left');

    %% === INITIAL DRAW ===
    refresh_all();
    redraw();
    set_status(sprintf('%d nodes, %d beams. Drag nodes to reshape. [Delete] to remove.', ...
               size(nodes,1), size(beams,1)));

    %% =================================================================
    %%  GEOMETRY GENERATOR
    %% =================================================================
    function [N, B] = generate_geometry(pr)
        N = zeros(0, 2);
        B = cell(0, 4);

        function idx = add_n(x, y)
            for k = 1:size(N,1)
                if abs(N(k,1)-x) < 0.5 && abs(N(k,2)-y) < 0.5
                    idx = k; return;
                end
            end
            N(end+1,:) = [x, y];
            idx = size(N,1);
        end
        function add_b(i, j, typ, col)
            B(end+1,:) = {i, j, typ, col};
        end

        W = pr.width; H = pr.height; Hs = pr.stilt_height;
        Sc = pr.stilt_size; hw = Sc/2;
        sh = (H - Hs) / 50;  % 50 stories above stilts
        n_xfer = pr.transfer_stories;
        bb = Hs + n_xfer * sh;   % brace_base
        tier_h = pr.tier_stories * sh;

        % Ground
        ng1 = add_n(-30, 0); ng2 = add_n(W+30, 0);
        add_b(ng1, ng2, 'Ground', 'gray');

        % Stilt
        ns1 = add_n(W/2-hw, 0);  ns2 = add_n(W/2+hw, 0);
        ns3 = add_n(W/2-hw, Hs); ns4 = add_n(W/2+hw, Hs);
        add_b(ns1, ns3, 'Stilt', 'gold');
        add_b(ns2, ns4, 'Stilt', 'gold');
        add_b(ns1, ns2, 'Stilt', 'gold');
        add_b(ns3, ns4, 'Stilt', 'gold');

        % Building edges
        ne_bl = add_n(0, Hs);  ne_br = add_n(W, Hs);
        ne_tl = add_n(0, H);   ne_tr = add_n(W, H);
        add_b(ne_bl, ne_tl, 'Edge', 'edge');
        add_b(ne_br, ne_tr, 'Edge', 'edge');
        add_b(ne_bl, ne_br, 'Edge', 'edge');

        % Transfer zone
        if n_xfer > 0
            nx_l = add_n(0, bb);  nx_r = add_n(W, bb);
            nc   = add_n(W/2, Hs); % stilt center top
            add_b(nx_l, nx_r, 'Edge', 'edge');
            add_b(nc, nx_l, 'Transfer', 'magenta');
            add_b(nc, nx_r, 'Transfer', 'magenta');
        end

        % Chevron V-braces (6 tiers)
        prev_apex = add_n(W/2, bb);
        for t = 1:pr.n_tiers
            zt = bb + t * tier_h;
            if zt > H + 1, break; end
            ntl = add_n(0, zt);
            ntr = add_n(W, zt);
            add_b(ntl, prev_apex, 'V-Brace', 'orange');
            add_b(ntr, prev_apex, 'V-Brace', 'orange');
            add_b(ntl, ntr, 'Chord', 'chord');
            if t < pr.n_tiers
                prev_apex = add_n(W/2, zt);
            end
        end

        % Roof
        npeak = add_n(W/2, H + W/4);
        add_b(ne_tl, npeak, 'Roof', 'roof');
        add_b(ne_tr, npeak, 'Roof', 'roof');
        add_b(ne_tl, ne_tr, 'Edge', 'edge');
    end

    %% =================================================================
    %%  DRAWING
    %% =================================================================
    function redraw()
        cla(ax);

        % --- Reference image ---
        read_img_vals();
        if has_img && img_on
            ih = p.height * img_sc;
            asp = size(ref_img,2) / size(ref_img,1);
            iw = ih * asp;
            cx = p.width/2 + img_x;
            cy = p.height/2 + img_y;
            hI = image(ax,'XData',[cx-iw/2 cx+iw/2],'YData',[cy+ih/2 cy-ih/2],...
                       'CData',ref_img);
            set(hI,'AlphaData', img_op * ones(size(ref_img,1),size(ref_img,2)));
            set(hI,'HitTest','off','PickableParts','none');
        end

        % --- Floor lines (faint) ---
        if show_floors
            sh_f = (p.height - p.stilt_height) / 50;
            for f = 0:49
                fz = p.stilt_height + f * sh_f;
                h_fl = plot(ax,[0 p.width],[fz fz],'-','Color',[.13 .13 .13],'LineWidth',0.3);
                set(h_fl,'HitTest','off','PickableParts','none');
                if mod(f,5) == 0
                    text(ax, p.width+5, fz, sprintf('%d', f+10), ...
                         'Color',[.25 .25 .25],'FontSize',7,'HitTest','off');
                end
            end
        end

        % --- Draw beams ---
        for i = 1:size(beams,1)
            ni = beams{i,1}; nj = beams{i,2};
            if ni < 1 || ni > size(nodes,1) || nj < 1 || nj > size(nodes,1)
                continue;
            end
            x1 = nodes(ni,1); y1 = nodes(ni,2);
            x2 = nodes(nj,1); y2 = nodes(nj,2);
            rgb = get_rgb(beams{i,4});
            lw = get_lw(beams{i,3});

            % White highlight for selected beam
            if i == sel_beam
                hh = plot(ax,[x1 x2],[y1 y2],'-','Color',[1 1 1],'LineWidth',lw+3);
                set(hh,'HitTest','off','PickableParts','none');
            end
            hb = plot(ax,[x1 x2],[y1 y2],'-','Color',rgb,'LineWidth',lw);
            set(hb,'HitTest','off','PickableParts','none');
        end

        % --- Draw nodes (purple draggable dots) ---
        for i = 1:size(nodes,1)
            nx = nodes(i,1); ny = nodes(i,2);
            if i == sel_node
                % Selected node: larger, yellow ring
                plot(ax, nx, ny, 'o', 'Color',[1 1 0], 'MarkerSize',14, ...
                     'LineWidth',2.5, 'HitTest','off','PickableParts','none');
            end
            plot(ax, nx, ny, 'o', 'MarkerFaceColor',[.55 .4 .8], ...
                 'MarkerEdgeColor',[.3 .2 .5], 'MarkerSize',8, 'LineWidth',1.2, ...
                 'HitTest','off','PickableParts','none');

            % Node ID label
            if show_ids
                text(ax, nx+4, ny+8, sprintf('%d', i), ...
                     'Color',[.6 .5 .8],'FontSize',7,'HitTest','off');
            end
        end

        % --- Dimension lines ---
        W = p.width; H = p.height; Hs = p.stilt_height;
        dc = [.35 .35 .35];
        % Height
        plot(ax,[W+15 W+15],[0 H],'-','Color',dc,'LineWidth',.5,'HitTest','off','PickableParts','none');
        text(ax,W+20,H/2,sprintf('%.0f ft',H),'Color',dc,'FontSize',8,'Rotation',90,...
             'HorizontalAlignment','center','HitTest','off');
        % Stilt height
        plot(ax,[W+30 W+30],[0 Hs],'-','Color',dc,'LineWidth',.5,'HitTest','off','PickableParts','none');
        text(ax,W+35,Hs/2,sprintf('%.0f ft',Hs),'Color',dc,'FontSize',7,'Rotation',90,...
             'HorizontalAlignment','center','HitTest','off');
        % Width
        plot(ax,[0 W],[-20 -20],'-','Color',dc,'LineWidth',.5,'HitTest','off','PickableParts','none');
        text(ax,W/2,-30,sprintf('%.0f ft',W),'Color',dc,'FontSize',8,...
             'HorizontalAlignment','center','HitTest','off');

        % Axes formatting
        rp = H + W/4;
        axis(ax,'equal');
        xlim(ax,[-50 W+55]); ylim(ax,[-45 rp+25]);
        title(ax,'Citicorp Center — Drag Nodes to Edit','Color','w','FontSize',13);
        xlabel(ax,'X (ft)','Color','w'); ylabel(ax,'Y (ft)','Color','w');

        % Keep axes click working
        set(ax,'ButtonDownFcn',@axes_click);
    end

    %% =================================================================
    %%  REFRESH TABLES
    %% =================================================================
    function refresh_all()
        refresh_node_tbl();
        refresh_beam_tbl();
        refresh_floor_tbl();
    end

    function refresh_node_tbl()
        nn = size(nodes,1);
        D = cell(nn, 3);
        for i = 1:nn
            D{i,1} = i;
            D{i,2} = round(nodes(i,1), 1);
            D{i,3} = round(nodes(i,2), 1);
        end
        set(node_tbl, 'Data', D);
    end

    function refresh_beam_tbl()
        nb = size(beams,1);
        D = cell(nb, 5);
        for i = 1:nb
            D{i,1} = i;
            D{i,2} = beams{i,1};
            D{i,3} = beams{i,2};
            D{i,4} = beams{i,3};
            D{i,5} = beams{i,4};
        end
        set(beam_tbl, 'Data', D);
    end

    function refresh_floor_tbl()
        W = p.width; H = p.height; Hs = p.stilt_height;
        sh_f = (H - Hs) / 50;
        D = cell(50, 6);
        for f = 0:49
            fz = Hs + f * sh_f;
            D{f+1,1} = f + 10;       % Floor number (10-59)
            D{f+1,2} = 0;            % Left X
            D{f+1,3} = round(fz,1);  % Left Y
            D{f+1,4} = W;            % Right X
            D{f+1,5} = round(fz,1);  % Right Y
            D{f+1,6} = round(fz,1);  % Height in ft
        end
        set(floor_tbl, 'Data', D);
    end

    %% =================================================================
    %%  MOUSE CALLBACKS
    %% =================================================================
    function mouse_down(~, ~)
        cp = get(ax,'CurrentPoint');
        mx = cp(1,1); my = cp(1,2);
        if ~in_ax(mx, my), return; end

        switch mode
            case 'select'
                % Try to grab a node first
                [dn, ni] = nearest_node(mx, my);
                if dn < grab_r
                    sel_node = ni;
                    sel_beam = 0;
                    dragging = true;
                    drag_idx = ni;
                    set_status(sprintf('Node %d: (%.1f, %.1f)  — drag to move', ...
                               ni, nodes(ni,1), nodes(ni,2)));
                else
                    % Try to select a beam
                    [db, bi] = nearest_beam(mx, my);
                    if db < grab_r
                        sel_beam = bi;
                        sel_node = 0;
                        set_status(sprintf('Beam %d: %s [%s] nodes %d→%d', ...
                                   bi, beams{bi,3}, beams{bi,4}, beams{bi,1}, beams{bi,2}));
                    else
                        sel_node = 0;
                        sel_beam = 0;
                        set_status(sprintf('(%.1f, %.1f)', mx, my));
                    end
                end
                redraw();

            case 'add_node'
                [sx, sy] = do_snap(mx, my);
                nodes(end+1,:) = [sx, sy];
                sel_node = size(nodes,1);
                sel_beam = 0;
                mode = 'select';
                refresh_all(); redraw();
                set_status(sprintf('Node %d added at (%.1f, %.1f)', sel_node, sx, sy));

            case 'beam_from'
                [dn, ni] = nearest_node(mx, my);
                if dn < grab_r
                    beam_from_node = ni;
                    mode = 'beam_to';
                    set_status(sprintf('From node %d — now click destination node', ni));
                end

            case 'beam_to'
                [dn, ni] = nearest_node(mx, my);
                if dn < grab_r && ni ~= beam_from_node
                    beams(end+1,:) = {beam_from_node, ni, 'Custom', 'orange'};
                    sel_beam = size(beams,1);
                    sel_node = 0;
                    mode = 'select';
                    refresh_all(); redraw();
                    set_status(sprintf('Beam %d added: node %d → %d', ...
                               sel_beam, beam_from_node, ni));
                end
        end
    end

    function axes_click(~, ~)
        mouse_down([], []);
    end

    function mouse_move(~, ~)
        cp = get(ax,'CurrentPoint');
        mx = cp(1,1); my = cp(1,2);

        % Coordinate display with floor estimate
        sh_f = (p.height - p.stilt_height) / 50;
        if my >= p.stilt_height && my <= p.height
            fl = floor((my - p.stilt_height) / sh_f) + 10;
            fl = min(fl, 59);
            set(coord_txt,'String',sprintf('  X = %.1f ft   Y = %.1f ft   (~Floor %d)', mx, my, fl));
        else
            set(coord_txt,'String',sprintf('  X = %.1f ft   Y = %.1f ft', mx, my));
        end

        % Drag node
        if dragging && drag_idx > 0
            [sx, sy] = do_snap(mx, my);
            nodes(drag_idx,:) = [sx, sy];
            redraw();
        end
    end

    function mouse_up(~, ~)
        if dragging
            dragging = false;
            drag_idx = 0;
            refresh_all();
        end
    end

    function key_press(~, evt)
        switch evt.Key
            case 'delete'
                cb_delete([],[]);
            case 'escape'
                mode = 'select';
                sel_node = 0; sel_beam = 0;
                redraw();
                set_status('Selection cleared.');
        end
    end

    %% =================================================================
    %%  BUTTON CALLBACKS
    %% =================================================================
    function cb_sel(~,~)
        mode = 'select';
        set_status('SELECT mode: click to select, drag to move nodes.');
    end

    function cb_addnode(~,~)
        mode = 'add_node';
        set_status('ADD NODE: click on the drawing to place a new node.');
    end

    function cb_addbeam(~,~)
        mode = 'beam_from';
        beam_from_node = 0;
        set_status('ADD BEAM: click the FROM node...');
    end

    function cb_delete(~,~)
        if sel_node > 0 && sel_node <= size(nodes,1)
            % Remove beams connected to this node
            keep = true(size(beams,1),1);
            for b = 1:size(beams,1)
                if beams{b,1} == sel_node || beams{b,2} == sel_node
                    keep(b) = false;
                end
            end
            beams = beams(keep,:);
            % Decrement higher indices
            for b = 1:size(beams,1)
                if beams{b,1} > sel_node, beams{b,1} = beams{b,1} - 1; end
                if beams{b,2} > sel_node, beams{b,2} = beams{b,2} - 1; end
            end
            nodes(sel_node,:) = [];
            set_status(sprintf('Deleted node %d + connected beams. %d nodes, %d beams remain.', ...
                       sel_node, size(nodes,1), size(beams,1)));
            sel_node = 0; sel_beam = 0;
            refresh_all(); redraw();

        elseif sel_beam > 0 && sel_beam <= size(beams,1)
            beams(sel_beam,:) = [];
            set_status(sprintf('Deleted beam %d. %d beams remain.', sel_beam, size(beams,1)));
            sel_beam = 0;
            refresh_all(); redraw();

        else
            set_status('Nothing selected. Click a node or beam first.');
        end
    end

    function cb_generate(~,~)
        p.stilt_height = safe_num(get(e_stilt,'String'), 114);
        p.width = safe_num(get(e_width,'String'), 157);
        p.height = safe_num(get(e_height,'String'), 915);
        p.n_tiers = safe_num(get(e_tiers,'String'), 6);
        p.tier_stories = safe_num(get(e_stpt,'String'), 8);
        p.transfer_stories = safe_num(get(e_xfer,'String'), 2);
        [nodes, beams] = generate_geometry(p);
        sel_node = 0; sel_beam = 0;
        refresh_all(); redraw();
        set_status(sprintf('Generated: %d nodes, %d beams.', size(nodes,1), size(beams,1)));
    end

    function cb_export(~,~)
        fprintf('\n===== NODES (%d) =====\n', size(nodes,1));
        fprintf('%4s %10s %10s\n', '#', 'X (ft)', 'Y (ft)');
        for i = 1:size(nodes,1)
            fprintf('%4d %10.1f %10.1f\n', i, nodes(i,1), nodes(i,2));
        end
        fprintf('\n===== BEAMS (%d) =====\n', size(beams,1));
        fprintf('%4s %6s %6s %-10s %-8s\n', '#', 'From', 'To', 'Type', 'Color');
        for i = 1:size(beams,1)
            fprintf('%4d %6d %6d %-10s %-8s\n', i, beams{i,1}, beams{i,2}, beams{i,3}, beams{i,4});
        end
        fprintf('\n===== FLOOR SCHEDULE (50 floors) =====\n');
        sh_f = (p.height - p.stilt_height) / 50;
        fprintf('%6s %10s %10s %10s %10s\n', 'Floor', 'Left X', 'Left Y', 'Right X', 'Right Y');
        for f = 0:49
            fz = p.stilt_height + f * sh_f;
            fprintf('%6d %10.1f %10.1f %10.1f %10.1f\n', f+10, 0, fz, p.width, fz);
        end
        fprintf('\n--- Parameters for run_citicorp_simulation.m ---\n');
        fprintf('params.stilt_height     = %g;\n', p.stilt_height);
        fprintf('params.width            = %g;\n', p.width);
        fprintf('params.height           = %g;\n', p.height);
        fprintf('params.n_tiers          = %g;\n', p.n_tiers);
        fprintf('params.tier_stories     = %g;\n', p.tier_stories);
        fprintf('params.transfer_stories = %g;\n', p.transfer_stories);

        save_path = fullfile(pwd, 'custom_structure.mat');
        node_data = nodes; beam_data = beams; params_used = p; %#ok<NASGU>
        save(save_path, 'node_data', 'beam_data', 'params_used');
        fprintf('Saved to: %s\n\n', save_path);
        set_status(sprintf('Exported %d nodes + %d beams to console + custom_structure.mat', ...
                   size(nodes,1), size(beams,1)));
    end

    function cb_load_img(~,~)
        [fn, fp] = uigetfile({'*.png;*.jpg;*.jpeg;*.bmp','Images'},'Select image');
        if fn ~= 0
            ref_img = imread(fullfile(fp, fn));
            has_img = true; img_on = true;
            set(chk_img,'Value',1);
            redraw();
            set_status(sprintf('Loaded: %s', fn));
        end
    end

    %% =================================================================
    %%  TABLE CALLBACKS
    %% =================================================================
    function cb_node_edit(~, evt)
        r = evt.Indices(1); c = evt.Indices(2);
        v = evt.NewData;
        if c == 2 && isnumeric(v), nodes(r,1) = v; end
        if c == 3 && isnumeric(v), nodes(r,2) = v; end
        sel_node = r; sel_beam = 0;
        redraw();
    end

    function cb_node_sel(~, evt)
        if ~isempty(evt.Indices)
            sel_node = evt.Indices(1); sel_beam = 0;
            redraw();
            set_status(sprintf('Node %d: (%.1f, %.1f)', sel_node, ...
                       nodes(sel_node,1), nodes(sel_node,2)));
        end
    end

    function cb_beam_edit(~, evt)
        r = evt.Indices(1); c = evt.Indices(2);
        v = evt.NewData;
        if c == 2 && isnumeric(v), beams{r,1} = round(v); end
        if c == 3 && isnumeric(v), beams{r,2} = round(v); end
        if c == 4, beams{r,3} = v; end
        if c == 5, beams{r,4} = v; end
        sel_beam = r; sel_node = 0;
        redraw();
    end

    function cb_beam_sel(~, evt)
        if ~isempty(evt.Indices)
            sel_beam = evt.Indices(1); sel_node = 0;
            redraw();
            set_status(sprintf('Beam %d: %s [%s] nodes %d→%d', sel_beam, ...
                       beams{sel_beam,3}, beams{sel_beam,4}, beams{sel_beam,1}, beams{sel_beam,2}));
        end
    end

    %% =================================================================
    %%  TOGGLES
    %% =================================================================
    function read_toggles()
        show_ids = get(chk_ids,'Value');
        show_floors = get(chk_flr,'Value');
        img_on = get(chk_img,'Value');
        if get(chk_snap,'Value')
            snap_val = safe_num(get(e_snap,'String'), 5);
        else
            snap_val = 0;
        end
    end

    function read_toggles_redraw()
        read_toggles();
        redraw();
    end

    function read_img_vals()
        img_on = get(chk_img,'Value');
        img_x  = safe_num(get(e_imgx,'String'), 0);
        img_y  = safe_num(get(e_imgy,'String'), 0);
        img_sc = max(0.1, safe_num(get(e_imgsc,'String'), 1));
        img_op = max(0, min(1, safe_num(get(e_imgop,'String'), 0.3)));
    end

    %% =================================================================
    %%  UTILITIES
    %% =================================================================
    function [sx, sy] = do_snap(x, y)
        read_toggles();
        if snap_val > 0
            sx = round(x / snap_val) * snap_val;
            sy = round(y / snap_val) * snap_val;
        else
            sx = x; sy = y;
        end
    end

    function [d, idx] = nearest_node(mx, my)
        d = inf; idx = 0;
        for i = 1:size(nodes,1)
            dd = sqrt((nodes(i,1)-mx)^2 + (nodes(i,2)-my)^2);
            if dd < d, d = dd; idx = i; end
        end
    end

    function [d, idx] = nearest_beam(mx, my)
        d = inf; idx = 0;
        for i = 1:size(beams,1)
            ni = beams{i,1}; nj = beams{i,2};
            if ni < 1 || ni > size(nodes,1) || nj < 1 || nj > size(nodes,1)
                continue;
            end
            dd = pt2seg(mx, my, nodes(ni,1), nodes(ni,2), nodes(nj,1), nodes(nj,2));
            if dd < d, d = dd; idx = i; end
        end
    end

    function ok = in_ax(mx, my)
        xl = get(ax,'XLim'); yl = get(ax,'YLim');
        ok = mx >= xl(1) && mx <= xl(2) && my >= yl(1) && my <= yl(2);
    end

    function d = pt2seg(px, py, x1, y1, x2, y2)
        dx = x2-x1; dy = y2-y1;
        L2 = dx*dx + dy*dy;
        if L2 < 1e-10
            d = sqrt((px-x1)^2 + (py-y1)^2); return;
        end
        t = max(0, min(1, ((px-x1)*dx + (py-y1)*dy) / L2));
        d = sqrt((px - x1 - t*dx)^2 + (py - y1 - t*dy)^2);
    end

    function rgb = get_rgb(name)
        switch lower(name)
            case 'orange',  rgb = [1.0  0.55 0.0];
            case 'gold',    rgb = [0.9  0.6  0.1];
            case 'edge',    rgb = [0.6  0.4  0.15];
            case 'magenta', rgb = [1.0  0.4  0.8];
            case 'roof',    rgb = [0.7  0.7  0.7];
            case 'gray',    rgb = [0.35 0.35 0.35];
            case 'chord',   rgb = [0.42 0.28 0.0];
            case 'red',     rgb = [1.0  0.2  0.2];
            case 'green',   rgb = [0.2  1.0  0.3];
            case 'blue',    rgb = [0.3  0.5  1.0];
            case 'cyan',    rgb = [0.0  0.9  0.9];
            case 'white',   rgb = [1.0  1.0  1.0];
            case 'yellow',  rgb = [1.0  1.0  0.2];
            otherwise,      rgb = [1.0  0.55 0.0];
        end
    end

    function lw = get_lw(typ)
        switch lower(typ)
            case 'stilt',    lw = 4;
            case 'v-brace',  lw = 2.5;
            case 'transfer', lw = 2.5;
            case 'edge',     lw = 1.5;
            case 'chord',    lw = 0.8;
            case 'roof',     lw = 2;
            case 'ground',   lw = 1;
            otherwise,       lw = 2;
        end
    end

    function v = safe_num(s, default)
        v = str2double(s);
        if isnan(v), v = default; end
    end

    function set_status(msg)
        set(status_txt,'String',['  ' msg]);
    end

    %% === UI HELPERS ===
    function mk_hdr(x, y, str, col)
        uicontrol(fig,'Style','text','String',str, ...
            'Units','normalized','Position',[x y rpw 0.022], ...
            'BackgroundColor',bg,'ForegroundColor',col, ...
            'FontSize',10,'FontWeight','bold','HorizontalAlignment','left');
    end

    function mk_lbl(x, y, str)
        uicontrol(fig,'Style','text','String',str, ...
            'Units','normalized','Position',[x y 0.035 0.02], ...
            'BackgroundColor',bg,'ForegroundColor',[.6 .6 .6], ...
            'FontSize',8,'HorizontalAlignment','left');
    end

    function h = mk_ed(x, y, w, str)
        h = uicontrol(fig,'Style','edit','String',str, ...
            'Units','normalized','Position',[x y w 0.023], ...
            'BackgroundColor',[.15 .15 .15],'ForegroundColor',[1 .8 .3], ...
            'FontSize',9);
    end
end
