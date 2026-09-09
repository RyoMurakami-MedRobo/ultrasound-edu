function main_gui()
%MAIN_GUI  Interactive DAS beamforming explorer for ultrasound (uifigure based).
%
%   MAIN_GUI
%
%   An educational tool built on MUST (Matlab UltraSound Toolbox). It lets you
%   change the transducer, the transmit sequence and the phantom layout
%   interactively while visualising and verifying what a Delay-And-Sum
%   beamformer actually does. When MUST is not installed, SIM_ENGINE falls
%   back to an analytic mock backend so every feature stays usable.
%
%   Modules:
%     SIM_ENGINE          RF generation (SIMUS wrapper / analytic mock)
%     DAS_REFERENCE       Textbook DAS reference implementation
%     DAS_CUSTOM_TEMPLATE Skeleton for user-written algorithms
%     WAVE_ANIMATOR       Wave propagation, delay curve and alignment views
%
%   Unified interface for custom beamformers:
%     [bmode_img, delays] = my_das(rf_data, tx_info, rx_pos, ...
%                                  grid_x, grid_z, sound_speed, fs)
%
%   The control column has a Simple / Detailed mode. Simple mode hides the
%   advanced panels; the values they hold stay in effect.
%
%   See also SIM_ENGINE, DAS_REFERENCE, DAS_CUSTOM_TEMPLATE, WAVE_ANIMATOR.

%% ---------------- State ----------------
app = struct();
app.S        = [];      % output of SIM_ENGINE
app.imgA     = [];      % linear envelope of algorithm A
app.imgB     = [];
app.nameA    = '';
app.nameB    = '';
app.msA      = NaN;
app.msB      = NaN;
app.gx       = [];
app.gz       = [];
app.pt       = [0 20e-3];
app.cursorXZ = [0 20];  % last clicked phantom coordinate [mm]
app.playing  = false;
app.animT    = 0;

% Wall-clock duration of one full animation sweep at 1x speed [s].
% The playback loop is paced against real time rather than frame count,
% because drawnow('limitrate') drops frames and a frame-counted loop would
% run at whatever speed the machine happens to render at.
SWEEP_SECONDS = 20;

ui = struct();
buildUI();
onSchemeChanged();
applyMode();
applyPreset();
drawPhantom();
setStatus(['Ready. Press [1] Simulate, then [2] Beamform / Compare.' ...
           '   MUST: ' mustStatusText()]);

%% ======================================================================
%  UI construction
%% ======================================================================
    function buildUI()
        ui.fig = uifigure('Name', 'Ultrasound DAS Beamforming Explorer', ...
                          'Position', [40 40 1580 940]);
        ui.fig.CloseRequestFcn = @onClose;

        g = uigridlayout(ui.fig, [2 2]);
        g.RowHeight    = {'1x', 24};
        g.ColumnWidth  = {380, '1x'};
        g.Padding = [6 6 6 6]; g.RowSpacing = 4; g.ColumnSpacing = 6;

        % ---------- Left: controls ----------
        % Row order must match ui.rowH and ui.advRows below.
        ui.rowH    = {30, 96, 144, 152, 120, 284, 96, 194, 96, 96, 88};
        ui.advRows = [3 5 8 10];          % rows collapsed in Simple mode

        ui.ctrl = uigridlayout(g, [numel(ui.rowH) 1]);
        ui.ctrl.Layout.Row = 1; ui.ctrl.Layout.Column = 1;
        ui.ctrl.RowHeight = ui.rowH;
        ui.ctrl.Scrollable = 'on';
        ui.ctrl.Padding = [2 2 2 2]; ui.ctrl.RowSpacing = 4;

        buildModeHeader(ui.ctrl);
        ui.advPanels = gobjects(0);
        buildProbeBasic(ui.ctrl);
        buildProbeAdvanced(ui.ctrl);
        buildTxBasic(ui.ctrl);
        buildTxAdvanced(ui.ctrl);
        buildPhantomPanel(ui.ctrl);
        buildReconBasic(ui.ctrl);
        buildReconAdvanced(ui.ctrl);
        buildAlgoBasic(ui.ctrl);
        buildAlgoAdvanced(ui.ctrl);
        buildRunPanel(ui.ctrl);

        % ---------- Right: tabs ----------
        ui.tabs = uitabgroup(g);
        ui.tabs.Layout.Row = 1; ui.tabs.Layout.Column = 2;
        buildTabImage();
        buildTabCompare();
        buildTabAnim();
        buildTabDelay();

        % ---------- Status bar ----------
        ui.status = uilabel(g, 'Text', '', 'FontSize', 12, ...
                            'BackgroundColor', [0.94 0.94 0.96]);
        ui.status.Layout.Row = 2; ui.status.Layout.Column = [1 2];
    end

    function [gg, p] = section(parent, ttl, heights)
        p = uipanel(parent, 'Title', ttl, 'FontWeight', 'bold', 'FontSize', 12);
        gg = uigridlayout(p, [numel(heights) 2]);
        gg.ColumnWidth = {172, '1x'};
        gg.RowHeight   = heights;
        gg.Padding = [6 4 6 4]; gg.RowSpacing = 3; gg.ColumnSpacing = 6;
    end

    function lab(parent, txt, r)
        lb = uilabel(parent, 'Text', txt, 'HorizontalAlignment', 'right');
        lb.Layout.Row = r; lb.Layout.Column = 1;
    end

    function h = put(h, r, c)
        h.Layout.Row = r; h.Layout.Column = c;
    end

    function markAdvanced(p)
        ui.advPanels(end+1) = p;
    end

    %% ---------------- Mode header ----------------
    function buildModeHeader(parent)
        gg = uigridlayout(parent, [1 2]);
        gg.ColumnWidth = {172, '1x'};
        gg.Padding = [6 2 6 2]; gg.ColumnSpacing = 6;
        uilabel(gg, 'Text', 'Control panel mode', 'HorizontalAlignment', 'right', ...
                'FontWeight', 'bold');
        ui.mode = uidropdown(gg, 'Items', {'Simple', 'Detailed'}, 'Value', 'Simple', ...
            'ValueChangedFcn', @(s,e) applyMode(), ...
            'Tooltip', 'Simple hides the advanced panels; their values stay in effect.');
    end

    %% ---------------- Transducer ----------------
    function buildProbeBasic(parent)
        gg = section(parent, 'Transducer', repmat({22}, 1, 2));
        lab(gg, 'Elements', 1);
        ui.nel   = put(uidropdown(gg, 'Items', {'16','32','64','128','192'}, ...
                        'Value', '64', 'Editable', 'on', ...
                        'ValueChangedFcn', @(s,e) drawPhantom()), 1, 2);
        lab(gg, 'Centre frequency [MHz]', 2);
        ui.fc    = put(uieditfield(gg, 'numeric', 'Value', 5, 'Limits', [0.5 30]), 2, 2);
    end

    function buildProbeAdvanced(parent)
        [gg, p] = section(parent, 'Transducer (advanced)', repmat({22}, 1, 4));
        markAdvanced(p);
        lab(gg, 'Pitch [mm]', 1);
        ui.pitch = put(uieditfield(gg, 'numeric', 'Value', 0.30, ...
                        'Limits', [0.01 5], 'ValueDisplayFormat', '%.3f', ...
                        'ValueChangedFcn', @(s,e) drawPhantom()), 1, 2);
        lab(gg, 'Bandwidth -6 dB [%]', 2);
        ui.bw    = put(uieditfield(gg, 'numeric', 'Value', 75, 'Limits', [5 200]), 2, 2);
        lab(gg, 'Speed of sound [m/s]', 3);
        ui.c     = put(uieditfield(gg, 'numeric', 'Value', 1540, 'Limits', [300 4000]), 3, 2);
        lab(gg, 'Sampling  fs / fc', 4);
        ui.fsfac = put(uidropdown(gg, 'Items', {'4','6','8','12'}, 'Value', '4'), 4, 2);
    end

    %% ---------------- Transmit ----------------
    function buildTxBasic(parent)
        gg = section(parent, 'Transmit scheme', {22, 40, 22, 22});
        lab(gg, 'Transmit mode', 1);
        ui.scheme = put(uidropdown(gg, 'Items', ...
            {'Plane wave', 'Focused', 'Diverging'}, ...
            'Value', 'Plane wave', 'ValueChangedFcn', @(s,e) onSchemeChanged()), 1, 2);

        lab(gg, 'Steering angle [deg]', 2);
        ui.angleSld = put(uislider(gg, 'Limits', [-40 40], 'Value', 0, ...
            'MajorTicks', -40:20:40, ...
            'ValueChangedFcn', @(s,e) syncAngle('slider')), 2, 2);
        lab(gg, 'Angle (numeric) [deg]', 3);
        ui.angle = put(uieditfield(gg, 'numeric', 'Value', 0, 'Limits', [-40 40], ...
            'ValueChangedFcn', @(s,e) syncAngle('edit')), 3, 2);

        lab(gg, 'Focal depth F [mm]', 4);
        ui.focus = put(uieditfield(gg, 'text', 'Value', '20', ...
            'Tooltip', ['Comma-separated values request a multi-focus ' ...
                        'sequence (composited by depth zone).']), 4, 2);
    end

    function buildTxAdvanced(parent)
        [gg, p] = section(parent, 'Transmit (advanced)', repmat({22}, 1, 3));
        markAdvanced(p);
        lab(gg, 'Diverging source [mm]', 1);
        ui.srcdepth = put(uieditfield(gg, 'numeric', 'Value', 10, 'Limits', [0.5 100], ...
            'Tooltip', 'Places a virtual source behind the array at z = -value'), 1, 2);
        lab(gg, 'Single-element transmit', 2);
        ui.singleel = put(uicheckbox(gg, 'Text', 'Fire one element only', 'Value', false, ...
            'ValueChangedFcn', @(s,e) onSchemeChanged()), 2, 2);
        lab(gg, 'Element index', 3);
        ui.elidx = put(uieditfield(gg, 'numeric', 'Value', 32, 'Limits', [1 1024], ...
            'RoundFractionalValues', 'on'), 3, 2);
    end

    %% ---------------- Phantom ----------------
    function buildPhantomPanel(parent)
        gg = section(parent, 'Targets (phantom)', {22, 26, 130, 26, 26});
        lab(gg, 'Preset', 1);
        ui.preset = put(uidropdown(gg, 'Items', ...
            {'Single point (PSF)', 'Wire phantom', 'Anechoic cyst', 'Custom'}, ...
            'Value', 'Single point (PSF)'), 1, 2);

        put(uibutton(gg, 'Text', 'Apply preset', ...
            'ButtonPushedFcn', @(s,e) applyPresetAndRefresh()), 2, 1);
        ui.clickAdd = put(uicheckbox(gg, 'Text', 'Left-click to add', 'Value', true), 2, 2);

        ui.tbl = uitable(gg, 'ColumnName', {'X [mm]', 'Z [mm]', 'Amplitude'}, ...
            'ColumnEditable', [true true true], 'ColumnWidth', {80, 80, 90}, ...
            'Data', [0 20 1], 'CellEditCallback', @(s,e) onTableEdit());
        ui.tbl.Layout.Row = 3; ui.tbl.Layout.Column = [1 2];

        put(uibutton(gg, 'Text', 'Add row', ...
            'ButtonPushedFcn', @(s,e) addRow()), 4, 1);
        put(uibutton(gg, 'Text', 'Delete selected', ...
            'ButtonPushedFcn', @(s,e) delSelectedRows()), 4, 2);
        put(uibutton(gg, 'Text', 'Clear all', ...
            'ButtonPushedFcn', @(s,e) clearRows()), 5, 1);
        put(uilabel(gg, 'Text', 'Right-click: add / delete', 'FontSize', 11, ...
            'FontColor', [.4 .4 .4]), 5, 2);
    end

    %% ---------------- Reconstruction ----------------
    function buildReconBasic(parent)
        gg = section(parent, 'Reconstruction', repmat({22}, 1, 2));
        lab(gg, 'Receive f-number (0 = full)', 1);
        ui.fnum = put(uieditfield(gg, 'numeric', 'Value', 1.5, 'Limits', [0 8]), 1, 2);
        lab(gg, 'Dynamic range [dB]', 2);
        ui.dr = put(uieditfield(gg, 'numeric', 'Value', 50, 'Limits', [10 120], ...
            'ValueChangedFcn', @(s,e) refreshImages()), 2, 2);
    end

    function buildReconAdvanced(parent)
        [gg, p] = section(parent, 'Reconstruction grid (advanced)', repmat({22}, 1, 6));
        markAdvanced(p);
        lab(gg, 'Pixels X', 1);
        ui.nx = put(uieditfield(gg, 'numeric', 'Value', 161, 'Limits', [8 1201], ...
                    'RoundFractionalValues', 'on'), 1, 2);
        lab(gg, 'Pixels Z', 2);
        ui.nz = put(uieditfield(gg, 'numeric', 'Value', 221, 'Limits', [8 1201], ...
                    'RoundFractionalValues', 'on'), 2, 2);
        lab(gg, 'Lateral half-width [mm]', 3);
        ui.xhalf = put(uieditfield(gg, 'numeric', 'Value', 12, 'Limits', [1 100], ...
            'ValueChangedFcn', @(s,e) drawPhantom()), 3, 2);
        lab(gg, 'Depth min [mm]', 4);
        ui.zmin = put(uieditfield(gg, 'numeric', 'Value', 3, 'Limits', [0.1 300]), 4, 2);
        lab(gg, 'Depth max [mm]', 5);
        ui.zmax = put(uieditfield(gg, 'numeric', 'Value', 40, 'Limits', [1 400], ...
            'ValueChangedFcn', @(s,e) drawPhantom()), 5, 2);
        lab(gg, 'Receive apodisation', 6);
        ui.rxapod = put(uidropdown(gg, 'Items', {'rect', 'hann'}, 'Value', 'rect'), 6, 2);
    end

    %% ---------------- Algorithms ----------------
    function buildAlgoBasic(parent)
        gg = section(parent, 'DAS algorithm comparison', repmat({22}, 1, 2));
        algItems = {'Reference DAS', 'Custom DAS', 'MUST das()'};
        lab(gg, 'Algorithm A', 1);
        ui.algA = put(uidropdown(gg, 'Items', algItems, 'Value', algItems{1}), 1, 2);
        lab(gg, 'Algorithm B', 2);
        ui.algB = put(uidropdown(gg, 'Items', [{'None'}, algItems], ...
                      'Value', 'Custom DAS'), 2, 2);
    end

    function buildAlgoAdvanced(parent)
        [gg, p] = section(parent, 'Algorithm (advanced)', repmat({22}, 1, 2));
        markAdvanced(p);
        lab(gg, 'Custom function name', 1);
        ui.customfn = put(uieditfield(gg, 'text', 'Value', 'das_custom_template', ...
            'Tooltip', 'Name of your copy of das_custom_template'), 1, 2);
        lab(gg, 'Backend', 2);
        ui.forcemock = put(uicheckbox(gg, 'Text', 'Force the mock backend', ...
            'Value', false), 2, 2);
    end

    %% ---------------- Run ----------------
    function buildRunPanel(parent)
        p  = uipanel(parent, 'Title', 'Run', 'FontWeight', 'bold', 'FontSize', 12);
        gg = uigridlayout(p, [2 1]);
        gg.RowHeight = {28, 28}; gg.Padding = [6 4 6 4]; gg.RowSpacing = 4;
        ui.btnSim = uibutton(gg, 'Text', '[1] Simulate (generate RF)', ...
            'FontWeight', 'bold', 'ButtonPushedFcn', @(s,e) onSimulate());
        ui.btnBf  = uibutton(gg, 'Text', '[2] Beamform / Compare', ...
            'FontWeight', 'bold', 'ButtonPushedFcn', @(s,e) onBeamform());
    end

    %% ---------------- Tab 1: images ----------------
    function buildTabImage()
        t  = uitab(ui.tabs, 'Title', 'Phantom & B-mode');
        gg = uigridlayout(t, [1 3]);
        gg.ColumnWidth = {'1x','1x','1x'}; gg.Padding = [8 8 8 8];
        ui.axPh = uiaxes(gg);  title(ui.axPh, 'Phantom layout');
        ui.axA  = uiaxes(gg);  title(ui.axA,  'Algorithm A');
        ui.axB  = uiaxes(gg);  title(ui.axB,  'Algorithm B');

        % The pointer has already moved onto the menu by the time an item is
        % selected, so capture the coordinate when the menu opens.
        ui.cmPh = uicontextmenu(ui.fig, ...
            'ContextMenuOpeningFcn', @(s,e) captureCursor());
        uimenu(ui.cmPh, 'Text', 'Add scatterer here', ...
               'MenuSelectedFcn', @(s,e) menuAddScat());
        uimenu(ui.cmPh, 'Text', 'Delete nearest scatterer', ...
               'MenuSelectedFcn', @(s,e) menuDelScat());
    end

    %% ---------------- Tab 2: comparison ----------------
    function buildTabCompare()
        t  = uitab(ui.tabs, 'Title', 'Compare & metrics');
        gg = uigridlayout(t, [2 2]);
        gg.RowHeight = {'1x', 168}; gg.ColumnWidth = {'1x','1x'};
        gg.Padding = [8 8 8 8];
        ui.axDiff = uiaxes(gg);   ui.axDiff.Layout.Row = 1; ui.axDiff.Layout.Column = 1;
        ui.axProf = uiaxes(gg);   ui.axProf.Layout.Row = 1; ui.axProf.Layout.Column = 2;
        ui.tblMet = uitable(gg, 'ColumnName', {'Metric', 'A', 'B'}, ...
            'ColumnWidth', {280, 180, 180}, 'Data', cell(0,3));
        ui.tblMet.Layout.Row = 2; ui.tblMet.Layout.Column = [1 2];
        title(ui.axDiff, 'Difference image |A - B|');
        title(ui.axProf, 'Lateral profile / FWHM');
    end

    %% ---------------- Tab 3: wave animation ----------------
    function buildTabAnim()
        t  = uitab(ui.tabs, 'Title', 'Wave animation');
        gg = uigridlayout(t, [2 1]);
        gg.RowHeight = {'1x', 40}; gg.Padding = [8 8 8 8];
        ui.axAnim = uiaxes(gg);
        cg = uigridlayout(gg, [1 6]);
        cg.ColumnWidth = {90, '1x', 120, 60, 90, 150};
        cg.Padding = [0 0 0 0]; cg.ColumnSpacing = 6;
        ui.btnPlay = uibutton(cg, 'Text', 'Play', ...
            'ButtonPushedFcn', @(s,e) onPlayToggle());
        ui.tslider = uislider(cg, 'Limits', [0 1], 'Value', 0, ...
            'MajorTicks', [], 'ValueChangingFcn', @(s,e) onScrub(e.Value), ...
            'ValueChangedFcn', @(s,e) onScrub(s.Value));
        ui.lblT = uilabel(cg, 'Text', 't = 0.00 us');
        uilabel(cg, 'Text', 'Speed', 'HorizontalAlignment', 'right');
        ui.speed = uidropdown(cg, 'Items', {'0.25x','0.5x','1x','2x','4x'}, 'Value', '1x', ...
            'Tooltip', 'Playback speed. 1x takes about 20 s for one full sweep.');
        ui.txsel = uidropdown(cg, 'Items', {'Transmit 1'}, 'Value', 'Transmit 1', ...
            'ValueChangedFcn', @(s,e) onTxEventChanged());
    end

    %% ---------------- Tab 4: delay & alignment ----------------
    function buildTabDelay()
        t  = uitab(ui.tabs, 'Title', 'Delay curve & alignment');
        gg = uigridlayout(t, [2 1]);
        gg.RowHeight = {'1x', 36}; gg.Padding = [8 8 8 8];
        ag = uigridlayout(gg, [1 3]);
        ag.ColumnWidth = {'1x','1x','1x'}; ag.Padding = [0 0 0 0];
        ui.axRF   = uiaxes(ag);
        ui.axPre  = uiaxes(ag);
        ui.axPost = uiaxes(ag);
        title(ui.axRF, 'RF data + delay curve');

        cg = uigridlayout(gg, [1 6]);
        cg.ColumnWidth = {260, 60, 80, 60, 80, '1x'};
        cg.Padding = [0 0 0 0]; cg.ColumnSpacing = 6;
        uilabel(cg, 'Text', 'Reconstruction point (or click the B-mode image):');
        uilabel(cg, 'Text', 'X [mm]', 'HorizontalAlignment', 'right');
        ui.ptx = uieditfield(cg, 'numeric', 'Value', 0, ...
            'ValueChangedFcn', @(s,e) onPointEdited());
        uilabel(cg, 'Text', 'Z [mm]', 'HorizontalAlignment', 'right');
        ui.ptz = uieditfield(cg, 'numeric', 'Value', 20, ...
            'ValueChangedFcn', @(s,e) onPointEdited());
        ui.lblDelayInfo = uilabel(cg, 'Text', '', 'FontColor', [.3 .3 .3]);
    end

%% ======================================================================
%  UI events
%% ======================================================================
    function applyMode()
        simple = strcmp(ui.mode.Value, 'Simple');
        h = ui.rowH;
        if simple
            for r = ui.advRows, h{r} = 0; end
        end
        ui.ctrl.RowHeight = h;
        for k = 1:numel(ui.advPanels)
            ui.advPanels(k).Visible = onoff(~simple);
        end
    end

    function onSchemeChanged()
        k = schemeKey();
        ui.focus.Enable    = onoff(strcmp(k, 'focused') && ~ui.singleel.Value);
        ui.srcdepth.Enable = onoff(strcmp(k, 'diverging') && ~ui.singleel.Value);
        ui.elidx.Enable    = onoff(ui.singleel.Value);
        ui.angleSld.Enable = onoff(~ui.singleel.Value);
        ui.angle.Enable    = onoff(~ui.singleel.Value);
    end

    function syncAngle(src)
        if strcmp(src, 'slider')
            ui.angle.Value = round(ui.angleSld.Value, 1);
        else
            ui.angleSld.Value = max(min(ui.angle.Value, 40), -40);
        end
    end

    function onTableEdit()
        ui.preset.Value = 'Custom';
        drawPhantom();
    end

    function addRow()
        ui.tbl.Data = [ui.tbl.Data; 0, 20, 1];
        ui.preset.Value = 'Custom';
        drawPhantom();
    end

    function delSelectedRows()
        sel = ui.tbl.Selection;
        D = ui.tbl.Data;
        if isempty(sel) || isempty(D)
            setStatus('No row selected.'); return;
        end
        D(unique(sel(:,1)), :) = [];
        ui.tbl.Data = D;
        ui.preset.Value = 'Custom';
        drawPhantom();
    end

    function clearRows()
        ui.tbl.Data = zeros(0, 3);
        ui.preset.Value = 'Custom';
        drawPhantom();
    end

    function applyPresetAndRefresh()
        applyPreset();
        drawPhantom();
        setStatus('Preset applied. Press [1] Simulate.');
    end

    function applyPreset()
        switch ui.preset.Value
            case 'Single point (PSF)'
                D = [0 20 1];
            case 'Wire phantom'
                [XW, ZW] = meshgrid(-8:4:8, 8:6:38);
                D = [XW(:), ZW(:), ones(numel(XW), 1)];
            case 'Anechoic cyst'
                rng(0);                                  % fixed seed for repeatability
                n  = 500;
                xb = (rand(n,1)*2 - 1) * 12;
                zb = 5 + rand(n,1) * 33;
                rb = randn(n,1) * 0.35;                  % speckle background
                cx = 0; cz = 22; rad = 5;                % anechoic disc
                rb(hypot(xb - cx, zb - cz) < rad) = 0;
                D = [xb, zb, rb];
                D = [D; -8 12 1; 8 32 1];                % point targets for reference
            otherwise                                    % Custom
                return;
        end
        ui.tbl.Data = D;
    end

    function captureCursor()
        %CAPTURECURSOR  Store the coordinate at the moment the context menu opens.
        cp = ui.axPh.CurrentPoint;
        app.cursorXZ = cp(1, 1:2);
    end

    function menuAddScat()
        % Anonymous functions capture values at creation time, so app must be
        % read inside a nested function (otherwise the coordinate would be
        % frozen at its initial value).
        addScatAt(app.cursorXZ);
    end

    function menuDelScat()
        delScatAt(app.cursorXZ);
    end

    function onPhantomClick(ax, ~)
        cp = ax.CurrentPoint;
        app.cursorXZ = cp(1, 1:2);
        if strcmp(ui.fig.SelectionType, 'alt')
            % Rarely reached: an assigned ContextMenu consumes the right-click.
            delScatAt(app.cursorXZ);
        elseif ui.clickAdd.Value
            addScatAt(app.cursorXZ);
        end
    end

    function addScatAt(xz)
        x = xz(1); z = xz(2);
        if z <= 0, setStatus('Click inside the region z > 0.'); return; end
        ui.tbl.Data = [ui.tbl.Data; x, z, 1];
        ui.preset.Value = 'Custom';
        drawPhantom();
        setStatus(sprintf('Scatterer added at (%.2f, %.2f) mm', x, z));
    end

    function delScatAt(xz)
        D = ui.tbl.Data;
        if isempty(D), return; end
        [dm, k] = min(hypot(D(:,1) - xz(1), D(:,2) - xz(2)));
        if dm > 3
            setStatus('No scatterer within 3 mm of the cursor.'); return;
        end
        D(k,:) = [];
        ui.tbl.Data = D;
        ui.preset.Value = 'Custom';
        drawPhantom();
        setStatus('Nearest scatterer deleted.');
    end

    function onPickPoint(ax, ~)
        cp = ax.CurrentPoint;
        app.pt = [cp(1,1)*1e-3, cp(1,2)*1e-3];
        ui.ptx.Value = cp(1,1);
        ui.ptz.Value = cp(1,2);
        refreshImages();
        updateDelayViews();
    end

    function onPointEdited()
        app.pt = [ui.ptx.Value*1e-3, ui.ptz.Value*1e-3];
        refreshImages();
        updateDelayViews();
    end

    function onTxEventChanged()
        setupAnimation();
        updateDelayViews();
    end

    function onClose(~, ~)
        app.playing = false;
        drawnow;
        delete(ui.fig);
    end

%% ======================================================================
%  Reading the settings
%% ======================================================================
    function k = schemeKey()
        switch ui.scheme.Value
            case 'Focused',   k = 'focused';
            case 'Diverging', k = 'diverging';
            otherwise,        k = 'plane';
        end
    end

    function cfg = readCfg()
        Nel = round(str2double(ui.nel.Value));
        if ~isfinite(Nel) || Nel < 2
            Nel = 64; ui.nel.Value = '64';
        end
        cfg.probe  = struct('Nelements', Nel, 'pitch', ui.pitch.Value*1e-3, ...
                            'fc', ui.fc.Value*1e6, 'bandwidth', ui.bw.Value);
        cfg.medium = struct('c', ui.c.Value);
        cfg.acq    = struct('fs_factor', str2double(ui.fsfac.Value));

        fmm = parseList(ui.focus.Value);
        if isempty(fmm), fmm = 20; end
        cfg.tx = struct('scheme', schemeKey(), 'angle_deg', ui.angle.Value, ...
                        'focus_mm', fmm, 'src_mm', ui.srcdepth.Value, ...
                        'single_element', ui.singleel.Value, ...
                        'element_index', ui.elidx.Value);

        D = ui.tbl.Data;
        if isempty(D), D = zeros(0,3); end
        cfg.scat = struct('x', D(:,1).'*1e-3, 'z', D(:,2).'*1e-3, 'rc', D(:,3).');

        cfg.recon = struct('zmax', ui.zmax.Value*1e-3, 'fnumber', ui.fnum.Value, ...
                           'rx_apod', ui.rxapod.Value);
        cfg.options = struct('force_mock', ui.forcemock.Value);
    end

    function [gx, gz] = reconGrid()
        xh = ui.xhalf.Value * 1e-3;
        gx = linspace(-xh, xh, round(ui.nx.Value));
        z0 = min(ui.zmin.Value, ui.zmax.Value - 1) * 1e-3;
        gz = linspace(max(z0, 1e-4), ui.zmax.Value*1e-3, round(ui.nz.Value));
    end

%% ======================================================================
%  Run: simulation
%% ======================================================================
    function onSimulate()
        cfg = readCfg();
        setStatus('Simulating...'); drawnow;
        try
            app.S = sim_engine(cfg);
        catch ME
            uialert(ui.fig, ME.message, 'Simulation failed');
            setStatus(['Simulation failed: ' ME.message]);
            return;
        end
        app.imgA = []; app.imgB = [];
        updateTxEventList();
        drawPhantom();
        setupAnimation();
        updateDelayViews();
        setStatus(sprintf(['RF generated | backend = %s | %d transmit event(s) | ' ...
            'RF %d x %d | fs = %.1f MHz | %.2f s'], app.S.backend, numel(app.S.tx), ...
            size(app.S.RF,1), size(app.S.RF,2), app.S.fs/1e6, app.S.elapsed));
    end

    function updateTxEventList()
        n = numel(app.S.tx);
        items = arrayfun(@(k) sprintf('Transmit %d', k), 1:n, 'UniformOutput', false);
        ui.txsel.Items = items;
        ui.txsel.Value = items{1};
    end

    function k = currentTx()
        k = 1;
        if ~isempty(app.S)
            idx = find(strcmp(ui.txsel.Items, ui.txsel.Value), 1);
            if ~isempty(idx), k = min(idx, numel(app.S.tx)); end
        end
    end

%% ======================================================================
%  Run: beamforming and comparison
%% ======================================================================
    function onBeamform()
        if isempty(app.S)
            onSimulate();
            if isempty(app.S), return; end
        end
        [gx, gz] = reconGrid();
        app.gx = gx; app.gz = gz;

        % The f-number and the receive apodisation may have changed after the
        % simulation, so push the current values into every transmit event.
        for k = 1:numel(app.S.tx)
            app.S.tx(k).fnumber = ui.fnum.Value;
            app.S.tx(k).rx_apod = ui.rxapod.Value;
        end

        setStatus('Beamforming...'); drawnow;
        try
            [app.imgA, app.msA, app.nameA] = runAlgorithm(ui.algA.Value, gx, gz);
        catch ME
            uialert(ui.fig, ME.message, 'Algorithm A failed');
            setStatus(['A failed: ' ME.message]); return;
        end

        app.imgB = []; app.msB = NaN; app.nameB = '';
        if ~strcmp(ui.algB.Value, 'None')
            try
                [app.imgB, app.msB, app.nameB] = runAlgorithm(ui.algB.Value, gx, gz);
            catch ME
                uialert(ui.fig, ME.message, 'Algorithm B failed');
                app.imgB = [];
            end
        end

        refreshImages();
        updateCompare();
        updateDelayViews();
        msg = sprintf('A: %s = %.1f ms', app.nameA, app.msA);
        if ~isempty(app.imgB)
            msg = sprintf('%s   |   B: %s = %.1f ms', msg, app.nameB, app.msB);
        end
        setStatus(['Beamforming done.   ' msg]);
    end

    function [img, ms, name] = runAlgorithm(key, gx, gz)
        [fh, name] = algoHandle(key);
        S = app.S;

        % JIT warm-up on a tiny grid: one throwaway call keeps the measured
        % time free of first-call compilation overhead.
        try
            fh(S.RF(:,:,1), S.tx(1), S.rx_pos, gx(1:min(4,end)), gz(1:min(4,end)), S.c, S.fs);
        catch
            % Ignore warm-up failures; the real call reports the error.
        end

        tic;
        img = runAllTx(fh, gx, gz);
        ms  = toc * 1000;

        expected = [numel(gz) numel(gx)];
        if ~isequal(size(img), expected)
            error('main_gui:imgSize', ...
                ['%s returned an image of size %s but %s was expected.\n' ...
                 'The unified interface requires a linear envelope of size\n' ...
                 '[numel(grid_z) x numel(grid_x)].'], name, ...
                 mat2str(size(img)), mat2str(expected));
        end
        if ~isreal(img) || any(~isfinite(img(:)))
            img = abs(img);
            img(~isfinite(img)) = 0;
        end
    end

    function img = runAllTx(fh, gx, gz)
        %RUNALLTX  A multi-focus sequence is composited zone by zone in depth.
        S   = app.S;
        Ntx = numel(S.tx);
        if Ntx == 1
            img = fh(S.RF(:,:,1), S.tx(1), S.rx_pos, gx, gz, S.c, S.fs);
            return;
        end
        fz = arrayfun(@(t) t.focus_mm, S.tx) * 1e-3;
        [fzs, ord] = sort(fz);
        edges = [-inf, (fzs(1:end-1) + fzs(2:end))/2, inf];
        img = zeros(numel(gz), numel(gx));
        for k = 1:Ntx
            sub = gz >= edges(k) & gz < edges(k+1);
            if ~any(sub), continue; end
            img(sub, :) = fh(S.RF(:,:,ord(k)), S.tx(ord(k)), S.rx_pos, gx, gz(sub), S.c, S.fs);
        end
    end

    function [fh, name] = algoHandle(key)
        switch key
            case 'Reference DAS'
                fh = @das_reference; name = 'das_reference';
            case 'Custom DAS'
                fname = strtrim(ui.customfn.Value);
                if isempty(which(fname))
                    error('main_gui:noCustom', ...
                        ['Custom function "%s" was not found.\n' ...
                         'Copy das_custom_template.m, rename the function to match\n' ...
                         'the file, and put it on the MATLAB path.'], fname);
                end
                fh = str2func(fname); name = fname;
            case 'MUST das()'
                if isempty(which('das'))
                    error('main_gui:noMUST', ...
                        ['das() from MUST was not found.\n' ...
                         'Install MUST (Matlab UltraSound Toolbox) and add it\n' ...
                         'to the MATLAB path.']);
                end
                fh = @mustDAS; name = 'MUST das()';
            otherwise
                error('main_gui:algo', 'Unknown algorithm "%s".', key);
        end
    end

    function [img, delays] = mustDAS(rf, tx, rx_pos, gx, gz, c, fs)
        %MUSTDAS  Adapter that exposes das() from MUST through the unified
        %  interface. Verified syntax: BFSIG = DAS(SIG, X, Z, DELAYS, PARAM)
        if isvector(gx) && isvector(gz)
            [XI, ZI] = meshgrid(gx(:).', gz(:).');
        else
            XI = gx; ZI = gz;
        end
        param = struct();
        param.fs        = fs;
        param.fc        = tx.fc;
        param.c         = c;
        param.pitch     = mean(diff(rx_pos));
        param.Nelements = numel(rx_pos);
        param.t0        = tx.t0;
        param.fnumber   = tx.fnumber;
        d = tx.delays; d(~isfinite(d)) = 0;
        bf = das(rf, XI, ZI, d, param);
        if isreal(bf) && size(bf,1) >= 4
            bf = hilbert(bf);
        end
        img = abs(bf);
        if nargout > 1
            % das() from MUST does not return delays, so borrow the ones from
            % the reference implementation.
            [~, delays] = das_reference(rf, tx, rx_pos, gx, gz, c, fs);
        end
    end

%% ======================================================================
%  Drawing
%% ======================================================================
    function drawPhantom()
        ax = ui.axPh;
        D  = ui.tbl.Data;
        cla(ax, 'reset');
        hold(ax, 'on');

        % Element row
        Nel = round(str2double(ui.nel.Value));
        if ~isfinite(Nel) || Nel < 2, Nel = 64; end
        ex = ((0:Nel-1) - (Nel-1)/2) * ui.pitch.Value;
        plot(ax, ex, zeros(1, Nel), 's', 'MarkerSize', 4, ...
             'MarkerEdgeColor', [.3 .3 .3], 'MarkerFaceColor', [.8 .8 .8], ...
             'PickableParts', 'none');

        if ~isempty(D)
            pos = D(:,3) >= 0;
            sz  = 20 + 60*min(abs(D(:,3)), 2)/2;
            scatter(ax, D(pos,1), D(pos,2), sz(pos), 'o', ...
                'MarkerEdgeColor', [0 .45 0], 'MarkerFaceColor', [.5 .9 .5], ...
                'PickableParts', 'none');
            scatter(ax, D(~pos,1), D(~pos,2), sz(~pos), 'o', ...
                'MarkerEdgeColor', [.6 0 0], 'MarkerFaceColor', [1 .7 .7], ...
                'PickableParts', 'none');
        end

        xh = max(ui.xhalf.Value, max(abs(ex))*1.1);
        xlim(ax, [-xh, xh]);
        ylim(ax, [-2, ui.zmax.Value]);
        set(ax, 'YDir', 'reverse');
        xlabel(ax, 'x [mm]'); ylabel(ax, 'z [mm]');
        title(ax, sprintf('Phantom layout (%d scatterers)', size(D,1)));
        grid(ax, 'on');
        hold(ax, 'off');
        ax.ButtonDownFcn = @onPhantomClick;
        ax.ContextMenu   = ui.cmPh;
    end

    function refreshImages()
        showBmode(ui.axA, app.imgA, app.nameA, true);
        showBmode(ui.axB, app.imgB, app.nameB, false);
    end

    function showBmode(ax, img, name, pickable)
        cla(ax, 'reset');
        if isempty(img)
            title(ax, '(not computed)');
            axis(ax, 'off');
            return;
        end
        axis(ax, 'on');
        dr = ui.dr.Value;
        L  = 20*log10(img / max(img(:)) + 1e-12);
        im = imagesc(ax, app.gx*1e3, app.gz*1e3, L);
        set(ax, 'CLim', [-dr 0], 'YDir', 'reverse');
        colormap(ax, gray(256));
        im.PickableParts = 'none';
        hold(ax, 'on');
        D = ui.tbl.Data;
        if ~isempty(D) && size(D,1) <= 60
            plot(ax, D(:,1), D(:,2), 'o', 'MarkerSize', 7, 'LineWidth', 1, ...
                 'MarkerEdgeColor', [.2 1 .2], 'PickableParts', 'none');
        end
        plot(ax, app.pt(1)*1e3, app.pt(2)*1e3, '+', 'Color', [1 .4 0], ...
             'MarkerSize', 12, 'LineWidth', 1.5, 'PickableParts', 'none');
        hold(ax, 'off');
        axis(ax, 'image');
        xlim(ax, [app.gx(1) app.gx(end)]*1e3);
        ylim(ax, [app.gz(1) app.gz(end)]*1e3);
        xlabel(ax, 'x [mm]'); ylabel(ax, 'z [mm]');
        title(ax, sprintf('%s  (DR %d dB)', name, round(dr)), 'Interpreter', 'none');
        colorbar(ax);
        if pickable
            ax.ButtonDownFcn = @onPickPoint;
        end
    end

    function updateCompare()
        A = app.imgA; B = app.imgB;
        cla(ui.axDiff, 'reset'); cla(ui.axProf, 'reset');
        if isempty(A)
            ui.tblMet.Data = cell(0,3); return;
        end
        nA = max(A(:)); if nA <= 0, nA = 1; end
        An = A / nA;

        rows = cell(0,3);
        [fwA, profA, pkA] = lateralMetrics(An, app.gx, app.gz);

        if ~isempty(B)
            Bn = B / nA;                       % normalise BOTH by the maximum of A
            Dm = abs(An - Bn);
            mse = mean((An(:) - Bn(:)).^2);
            imagesc(ui.axDiff, app.gx*1e3, app.gz*1e3, 20*log10(Dm + 1e-12));
            set(ui.axDiff, 'CLim', [-80 0], 'YDir', 'reverse');
            colormap(ui.axDiff, parula(256));
            colorbar(ui.axDiff);
            axis(ui.axDiff, 'image');
            xlabel(ui.axDiff, 'x [mm]'); ylabel(ui.axDiff, 'z [mm]');
            title(ui.axDiff, sprintf(['|A - B| [dB]  (both normalised by max(A))' ...
                '   MSE = %.3e'], mse));
            [fwB, profB, pkB] = lateralMetrics(Bn, app.gx, app.gz, pkA(2));
        else
            text(ui.axDiff, 0.5, 0.5, 'Select algorithm B to see the difference', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ui.axDiff, 'off');
            fwB = NaN; profB = []; pkB = [NaN NaN]; mse = NaN; Dm = [];
        end

        % ---- Lateral profile ----
        hold(ui.axProf, 'on');
        plot(ui.axProf, app.gx*1e3, todB(profA), '-', 'LineWidth', 1.6, ...
             'Color', [0 .35 .8], 'DisplayName', sprintf('A: %s', app.nameA));
        if ~isempty(profB)
            plot(ui.axProf, app.gx*1e3, todB(profB), '--', 'LineWidth', 1.6, ...
                 'Color', [.85 .2 .1], 'DisplayName', sprintf('B: %s', app.nameB));
        end
        plot(ui.axProf, [app.gx(1) app.gx(end)]*1e3, [-6 -6], ':', ...
             'Color', [.4 .4 .4], 'DisplayName', '-6 dB (FWHM)');
        hold(ui.axProf, 'off');
        grid(ui.axProf, 'on');
        ylim(ui.axProf, [-max(ui.dr.Value, 40) 2]);
        xlabel(ui.axProf, 'x [mm]'); ylabel(ui.axProf, 'Normalised amplitude [dB]');
        title(ui.axProf, sprintf('Lateral profile at z = %.2f mm', pkA(2)*1e3));
        legend(ui.axProf, 'Location', 'southwest', 'FontSize', 9, 'Interpreter', 'none');

        % ---- Metrics table ----
        rows(end+1,:) = {'Run time [ms]',              num2str(app.msA, '%.1f'),    numOrDash(app.msB, '%.1f')};
        rows(end+1,:) = {'Peak position x [mm]',       num2str(pkA(1)*1e3, '%.3f'), numOrDash(pkB(1)*1e3, '%.3f')};
        rows(end+1,:) = {'Peak position z [mm]',       num2str(pkA(2)*1e3, '%.3f'), numOrDash(pkB(2)*1e3, '%.3f')};
        rows(end+1,:) = {'Lateral FWHM (-6 dB) [mm]',  num2str(fwA*1e3, '%.4f'),    numOrDash(fwB*1e3, '%.4f')};
        rows(end+1,:) = {'MSE (normalised by max(A))', '-',                         numOrDash(mse, '%.4e')};
        if ~isempty(Dm)
            rows(end+1,:) = {'Max absolute error (vs A)', '-', num2str(max(Dm(:)), '%.4e')};
        end
        rows(end+1,:) = {'Algorithm', app.nameA, dashIfEmpty(app.nameB)};
        rows(end+1,:) = {'Backend',   app.S.backend, ''};
        ui.tblMet.Data = rows;
    end

    function [fw, prof, pk] = lateralMetrics(img, gx, gz, forceZ)
        %LATERALMETRICS  Lateral profile through the peak row and its -6 dB
        %  width, evaluated on linear values.
        if nargin >= 4 && ~isempty(forceZ) && isfinite(forceZ)
            [~, iz] = min(abs(gz - forceZ));
            [~, ix] = max(img(iz, :));
        else
            [~, k] = max(img(:));
            [iz, ix] = ind2sub(size(img), k);
        end
        prof = img(iz, :);
        pk   = [gx(ix), gz(iz)];
        fw   = fwhm(prof, gx);
    end

    function w = fwhm(prof, gx)
        %FWHM  Full width at half maximum (= -6 dB) of a linear envelope.
        %  The crossings are located by linear interpolation.
        [pkv, ip] = max(prof);
        if pkv <= 0, w = NaN; return; end
        half = pkv / 2;
        il = find(prof(1:ip) <= half, 1, 'last');
        ir = find(prof(ip:end) <= half, 1, 'first');
        if isempty(il) || isempty(ir), w = NaN; return; end
        ir = ir + ip - 1;
        xl = crossing(gx(il), prof(il), gx(il+1), prof(il+1), half);
        xr = crossing(gx(ir-1), prof(ir-1), gx(ir), prof(ir), half);
        w  = xr - xl;
    end

    function x = crossing(x1, y1, x2, y2, yt)
        if y2 == y1, x = x1; else, x = x1 + (yt - y1) * (x2 - x1) / (y2 - y1); end
    end

%% ======================================================================
%  Wave animation
%% ======================================================================
    function setupAnimation()
        if isempty(app.S), return; end
        xh = max(ui.xhalf.Value * 1e-3, max(abs(app.S.rx_pos))*1.05);
        tmax = wave_animator('setup_propagation', ui.axAnim, app.S, currentTx(), ...
                             [-xh xh], [0, ui.zmax.Value*1e-3]);
        ui.tslider.Limits = [0 tmax];
        ui.tslider.Value  = 0;
        app.animT = 0;
        wave_animator('draw_propagation', ui.axAnim, 0);
        ui.lblT.Text = 't = 0.00 us';
    end

    function onScrub(v)
        if isempty(app.S), return; end
        app.animT = v;
        wave_animator('draw_propagation', ui.axAnim, v);
        ui.lblT.Text = sprintf('t = %.2f us', v*1e6);
    end

    function onPlayToggle()
        if app.playing
            app.playing = false;
            return;
        end
        if isempty(app.S)
            setStatus('Press [1] Simulate first.'); return;
        end
        app.playing = true;
        ui.btnPlay.Text = 'Stop';
        spd    = str2double(strrep(ui.speed.Value, 'x', ''));
        tmax   = ui.tslider.Limits(2);
        dur    = SWEEP_SECONDS / spd;      % wall-clock seconds for one sweep
        tStart = app.animT;
        clk    = tic;
        while app.playing && isvalid(ui.fig)
            % Real-time pacing: the simulated time is derived from the elapsed
            % wall-clock time, so dropped frames slow the frame rate but never
            % speed up the wave.
            t = mod(tStart + tmax * toc(clk) / dur, tmax);
            app.animT = t;
            ui.tslider.Value = t;
            wave_animator('draw_propagation', ui.axAnim, t);
            ui.lblT.Text = sprintf('t = %.2f us', t*1e6);
            drawnow limitrate;
            pause(0.01);                   % yield and keep the CPU load sane
        end
        if isvalid(ui.fig)
            app.playing = false;
            ui.btnPlay.Text = 'Play';
        end
    end

%% ======================================================================
%  Delay curve / alignment
%% ======================================================================
    function updateDelayViews()
        if isempty(app.S), return; end
        k = currentTx();
        note = '';
        try
            [fh, name] = algoHandle(ui.algA.Value);
        catch
            fh = @das_reference; name = 'das_reference (fallback)';
        end
        try
            % Calling the beamformer on a single-point grid returns exactly
            % the samples it sums, so the overlay can never drift away from
            % what is actually computed.
            [~, tau] = fh(app.S.RF(:,:,k), app.S.tx(k), app.S.rx_pos, ...
                          app.pt(1), app.pt(2), app.S.c, app.S.fs);
        catch ME
            [~, tau] = das_reference(app.S.RF(:,:,k), app.S.tx(k), app.S.rx_pos, ...
                          app.pt(1), app.pt(2), app.S.c, app.S.fs);
            name = 'das_reference (fallback)';
            note = ['  - ' ME.message];
        end
        tau = tau(1, :);
        wave_animator('delay_curve', ui.axRF, app.S, k, tau, app.pt);
        wave_animator('alignment', ui.axPre, ui.axPost, app.S, k, tau, app.pt);
        ui.lblDelayInfo.Text = sprintf('Delays from: %s   %d / %d elements used%s', ...
            name, sum(isfinite(tau)), numel(tau), note);
    end

%% ======================================================================
%  Small helpers
%% ======================================================================
    function setStatus(msg)
        ui.status.Text = ['  ' msg];
        drawnow limitrate;
    end

    function s = mustStatusText()
        if isempty(which('simus'))
            s = 'not found (running on the mock backend)';
        else
            s = ['found (' which('simus') ')'];
        end
    end

    function v = parseList(str)
        v = str2double(strsplit(strtrim(str), {',', ' ', ';'}));
        v = v(isfinite(v) & v ~= 0);
    end

    function s = onoff(tf)
        if tf, s = 'on'; else, s = 'off'; end
    end

    function s = numOrDash(v, fmt)
        if isempty(v) || ~isfinite(v), s = '-'; else, s = num2str(v, fmt); end
    end

    function s = dashIfEmpty(v)
        if isempty(v), s = '-'; else, s = v; end
    end

    function y = todB(p)
        m = max(p);
        if m <= 0, m = 1; end
        y = 20*log10(p/m + 1e-12);
    end
end
