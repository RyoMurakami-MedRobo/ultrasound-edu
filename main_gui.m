function main_gui()
%MAIN_GUI  超音波 DAS ビームフォーミング学習・検証 GUI（uifigure ベース）。
%
%   MAIN_GUI
%
%   MUST (Matlab UltraSound Toolbox) をバックエンドに、トランスデューサ構成・
%   送信シーケンス・ファントム配置を対話的に変更しながら、DAS ビームフォーミング
%   の中身を可視化・検証するための教育用ツール。
%   MUST が見つからない場合は SIM_ENGINE が解析的モックへ自動フォールバック
%   するため、MUST 未インストール環境でも全機能を試せる。
%
%   構成モジュール:
%     SIM_ENGINE          RF データ生成（MUST simus ラッパ / モック）
%     DAS_REFERENCE       教育用の標準 DAS 実装（参照コード）
%     DAS_CUSTOM_TEMPLATE 自作アルゴリズム用ひな形
%     WAVE_ANIMATOR       波面伝搬・遅延曲線・整相前後比較の描画
%
%   自作 DAS の統一インターフェース:
%     [bmode_img, delays] = my_das(rf_data, tx_info, rx_pos, ...
%                                  grid_x, grid_z, sound_speed, fs)
%
%   See also SIM_ENGINE, DAS_REFERENCE, DAS_CUSTOM_TEMPLATE, WAVE_ANIMATOR.

%% ---------------- 状態 ----------------
app = struct();
app.S       = [];      % SIM_ENGINE の出力
app.imgA    = [];      % アルゴリズム A の線形エンベロープ
app.imgB    = [];
app.nameA   = '';
app.nameB   = '';
app.msA     = NaN;
app.msB     = NaN;
app.gx      = [];
app.gz      = [];
app.pt      = [0 20e-3];
app.cursorXZ = [0 20];   % 直近にクリックされたファントム座標 [mm]
app.playing = false;
app.animT   = 0;

ui = struct();
buildUI();
applyPreset();
drawPhantom();
setStatus(['準備完了。 [1] シミュレーション実行 → [2] ビームフォーミング/比較 ' ...
           'の順に押してください。   MUST 検出: ' mustStatusText()]);

%% ======================================================================
%  UI 構築
%% ======================================================================
    function buildUI()
        ui.fig = uifigure('Name', '超音波 DAS ビームフォーミング学習ツール', ...
                          'Position', [40 40 1580 940]);
        ui.fig.CloseRequestFcn = @onClose;

        g = uigridlayout(ui.fig, [2 2]);
        g.RowHeight    = {'1x', 24};
        g.ColumnWidth  = {380, '1x'};
        g.Padding = [6 6 6 6]; g.RowSpacing = 4; g.ColumnSpacing = 6;

        % ---------- 左：コントロール ----------
        ctrl = uigridlayout(g, [6 1]);
        ctrl.Layout.Row = 1; ctrl.Layout.Column = 1;
        ctrl.RowHeight = {192, 238, 302, 248, 148, 88};
        ctrl.Scrollable = 'on';
        ctrl.Padding = [2 2 2 2]; ctrl.RowSpacing = 6;

        buildProbePanel(ctrl);
        buildTxPanel(ctrl);
        buildPhantomPanel(ctrl);
        buildReconPanel(ctrl);
        buildAlgoPanel(ctrl);
        buildRunPanel(ctrl);

        % ---------- 右：タブ ----------
        ui.tabs = uitabgroup(g);
        ui.tabs.Layout.Row = 1; ui.tabs.Layout.Column = 2;
        buildTabImage();
        buildTabCompare();
        buildTabAnim();
        buildTabDelay();

        % ---------- ステータスバー ----------
        ui.status = uilabel(g, 'Text', '', 'FontSize', 12, ...
                            'BackgroundColor', [0.94 0.94 0.96]);
        ui.status.Layout.Row = 2; ui.status.Layout.Column = [1 2];
    end

    function gg = section(parent, ttl, heights)
        p = uipanel(parent, 'Title', ttl, 'FontWeight', 'bold', 'FontSize', 12);
        gg = uigridlayout(p, [numel(heights) 2]);
        gg.ColumnWidth = {160, '1x'};
        gg.RowHeight   = heights;
        gg.Padding = [6 4 6 4]; gg.RowSpacing = 3; gg.ColumnSpacing = 6;
    end

    function lb = lab(parent, txt, r)
        lb = uilabel(parent, 'Text', txt, 'HorizontalAlignment', 'right');
        lb.Layout.Row = r; lb.Layout.Column = 1;
    end

    function h = put(h, r, c)
        h.Layout.Row = r; h.Layout.Column = c;
    end

    %% ---------------- トランスデューサ ----------------
    function buildProbePanel(parent)
        gg = section(parent, 'トランスデューサ', repmat({22}, 1, 6));
        lab(gg, '素子数',            1);
        ui.nel   = put(uidropdown(gg, 'Items', {'16','32','64','128','192'}, ...
                        'Value', '64', 'Editable', 'on', ...
                        'ValueChangedFcn', @(s,e) drawPhantom()), 1, 2);
        lab(gg, 'ピッチ [mm]',       2);
        ui.pitch = put(uieditfield(gg, 'numeric', 'Value', 0.30, ...
                        'Limits', [0.01 5], 'ValueDisplayFormat', '%.3f', ...
                        'ValueChangedFcn', @(s,e) drawPhantom()), 2, 2);
        lab(gg, '中心周波数 [MHz]',  3);
        ui.fc    = put(uieditfield(gg, 'numeric', 'Value', 5, 'Limits', [0.5 30]), 3, 2);
        lab(gg, '比帯域 -6dB [%]',   4);
        ui.bw    = put(uieditfield(gg, 'numeric', 'Value', 75, 'Limits', [5 200]), 4, 2);
        lab(gg, '音速 [m/s]',        5);
        ui.c     = put(uieditfield(gg, 'numeric', 'Value', 1540, 'Limits', [300 4000]), 5, 2);
        lab(gg, 'サンプリング fs/fc',6);
        ui.fsfac = put(uidropdown(gg, 'Items', {'4','6','8','12'}, 'Value', '4'), 6, 2);
    end

    %% ---------------- 送信 ----------------
    function buildTxPanel(parent)
        gg = section(parent, '送信シーケンス (Tx Scheme)', {22, 40, 22, 22, 22, 22, 22});
        lab(gg, '送信モード', 1);
        ui.scheme = put(uidropdown(gg, 'Items', ...
            {'平面波 (Plane wave)', '集束波 (Focused)', '発散波 (Diverging)'}, ...
            'Value', '平面波 (Plane wave)', 'ValueChangedFcn', @(s,e) onSchemeChanged()), 1, 2);

        lab(gg, 'ステアリング角 \theta', 2);
        ui.angleSld = put(uislider(gg, 'Limits', [-40 40], 'Value', 0, ...
            'MajorTicks', -40:20:40, ...
            'ValueChangedFcn', @(s,e) syncAngle('slider')), 2, 2);
        lab(gg, '角度 [deg]', 3);
        ui.angle = put(uieditfield(gg, 'numeric', 'Value', 0, 'Limits', [-40 40], ...
            'ValueChangedFcn', @(s,e) syncAngle('edit')), 3, 2);

        lab(gg, '集束深度 F [mm]', 4);
        ui.focus = put(uieditfield(gg, 'text', 'Value', '20', ...
            'Tooltip', 'カンマ区切りで複数指定するとマルチフォーカス送信（深度ゾーン合成）になります'), 4, 2);

        lab(gg, '発散: 仮想音源 [mm]', 5);
        ui.srcdepth = put(uieditfield(gg, 'numeric', 'Value', 10, 'Limits', [0.5 100], ...
            'Tooltip', 'アレイ背後 z = -この値 に仮想音源を置きます'), 5, 2);

        lab(gg, '単一素子送信', 6);
        ui.singleel = put(uicheckbox(gg, 'Text', '1 素子だけを励振する', 'Value', false, ...
            'ValueChangedFcn', @(s,e) onSchemeChanged()), 6, 2);
        lab(gg, '励振素子番号', 7);
        ui.elidx = put(uieditfield(gg, 'numeric', 'Value', 32, 'Limits', [1 1024], ...
            'RoundFractionalValues', 'on'), 7, 2);
        onSchemeChanged();
    end

    %% ---------------- ファントム ----------------
    function buildPhantomPanel(parent)
        gg = section(parent, 'ターゲット（ファントム）配置', {22, 26, 148, 26, 26});
        lab(gg, 'プリセット', 1);
        ui.preset = put(uidropdown(gg, 'Items', ...
            {'単一点散乱体 (PSF)', 'ワイヤーファントム', '無エコー領域 (シスト)', 'カスタム'}, ...
            'Value', '単一点散乱体 (PSF)'), 1, 2);

        ui.btnPreset = put(uibutton(gg, 'Text', 'プリセット適用', ...
            'ButtonPushedFcn', @(s,e) applyPresetAndRefresh()), 2, 1);
        ui.clickAdd  = put(uicheckbox(gg, 'Text', '左クリックで追加', 'Value', true), 2, 2);

        ui.tbl = uitable(gg, 'ColumnName', {'X [mm]', 'Z [mm]', '反射強度'}, ...
            'ColumnEditable', [true true true], 'ColumnWidth', {80, 80, 80}, ...
            'Data', [0 20 1], 'CellEditCallback', @(s,e) onTableEdit());
        ui.tbl.Layout.Row = 3; ui.tbl.Layout.Column = [1 2];

        put(uibutton(gg, 'Text', '行を追加', ...
            'ButtonPushedFcn', @(s,e) addRow()), 4, 1);
        put(uibutton(gg, 'Text', '選択行を削除', ...
            'ButtonPushedFcn', @(s,e) delSelectedRows()), 4, 2);
        put(uibutton(gg, 'Text', '全消去', ...
            'ButtonPushedFcn', @(s,e) clearRows()), 5, 1);
        put(uilabel(gg, 'Text', '右クリック → 追加/削除', 'FontSize', 11, ...
            'FontColor', [.4 .4 .4]), 5, 2);
    end

    %% ---------------- 再構成 ----------------
    function buildReconPanel(parent)
        gg = section(parent, '再構成グリッド / 受信', repmat({22}, 1, 8));
        lab(gg, '画素数 X', 1);
        ui.nx = put(uieditfield(gg, 'numeric', 'Value', 161, 'Limits', [8 1201], ...
                    'RoundFractionalValues', 'on'), 1, 2);
        lab(gg, '画素数 Z', 2);
        ui.nz = put(uieditfield(gg, 'numeric', 'Value', 221, 'Limits', [8 1201], ...
                    'RoundFractionalValues', 'on'), 2, 2);
        lab(gg, '横方向範囲 ± [mm]', 3);
        ui.xhalf = put(uieditfield(gg, 'numeric', 'Value', 12, 'Limits', [1 100], ...
            'ValueChangedFcn', @(s,e) drawPhantom()), 3, 2);
        lab(gg, '深さ 最小 [mm]', 4);
        ui.zmin = put(uieditfield(gg, 'numeric', 'Value', 3, 'Limits', [0.1 300]), 4, 2);
        lab(gg, '深さ 最大 [mm]', 5);
        ui.zmax = put(uieditfield(gg, 'numeric', 'Value', 40, 'Limits', [1 400], ...
            'ValueChangedFcn', @(s,e) drawPhantom()), 5, 2);
        lab(gg, '受信 f 値 (0=全開口)', 6);
        ui.fnum = put(uieditfield(gg, 'numeric', 'Value', 1.5, 'Limits', [0 8]), 6, 2);
        lab(gg, '受信アポダイゼーション', 7);
        ui.rxapod = put(uidropdown(gg, 'Items', {'rect', 'hann'}, 'Value', 'rect'), 7, 2);
        lab(gg, 'ダイナミックレンジ [dB]', 8);
        ui.dr = put(uieditfield(gg, 'numeric', 'Value', 50, 'Limits', [10 120], ...
            'ValueChangedFcn', @(s,e) refreshImages()), 8, 2);
    end

    %% ---------------- アルゴリズム ----------------
    function buildAlgoPanel(parent)
        gg = section(parent, 'DAS アルゴリズム比較', repmat({22}, 1, 4));
        algItems = {'参照DAS (das_reference)', '自作DAS (Custom)', 'MUST das()'};
        lab(gg, 'アルゴリズム A', 1);
        ui.algA = put(uidropdown(gg, 'Items', algItems, 'Value', algItems{1}), 1, 2);
        lab(gg, 'アルゴリズム B', 2);
        ui.algB = put(uidropdown(gg, 'Items', [{'なし'}, algItems], ...
                      'Value', '自作DAS (Custom)'), 2, 2);
        lab(gg, '自作関数名', 3);
        ui.customfn = put(uieditfield(gg, 'text', 'Value', 'das_custom_template', ...
            'Tooltip', 'das_custom_template をコピーして作った関数名を入力'), 3, 2);
        lab(gg, 'バックエンド', 4);
        ui.forcemock = put(uicheckbox(gg, 'Text', 'MUST があってもモックを使う', ...
            'Value', false), 4, 2);
    end

    %% ---------------- 実行 ----------------
    function buildRunPanel(parent)
        p  = uipanel(parent, 'Title', '実行', 'FontWeight', 'bold', 'FontSize', 12);
        gg = uigridlayout(p, [2 1]);
        gg.RowHeight = {28, 28}; gg.Padding = [6 4 6 4]; gg.RowSpacing = 4;
        ui.btnSim = uibutton(gg, 'Text', '[1] シミュレーション実行 (RF 生成)', ...
            'FontWeight', 'bold', 'ButtonPushedFcn', @(s,e) onSimulate());
        ui.btnBf  = uibutton(gg, 'Text', '[2] ビームフォーミング / 比較', ...
            'FontWeight', 'bold', 'ButtonPushedFcn', @(s,e) onBeamform());
    end

    %% ---------------- タブ 1: 画像 ----------------
    function buildTabImage()
        t  = uitab(ui.tabs, 'Title', 'ファントム & B モード');
        gg = uigridlayout(t, [1 3]);
        gg.ColumnWidth = {'1x','1x','1x'}; gg.Padding = [8 8 8 8];
        ui.axPh = uiaxes(gg);  title(ui.axPh, 'ファントム配置');
        ui.axA  = uiaxes(gg);  title(ui.axA,  'アルゴリズム A');
        ui.axB  = uiaxes(gg);  title(ui.axB,  'アルゴリズム B');

        % メニュー項目を選ぶ頃にはポインタが移動してしまうため、
        % 右クリックでメニューが開いた瞬間の座標を確定して保持する。
        ui.cmPh = uicontextmenu(ui.fig, ...
            'ContextMenuOpeningFcn', @(s,e) captureCursor());
        uimenu(ui.cmPh, 'Text', 'ここに散乱体を追加', ...
               'MenuSelectedFcn', @(s,e) menuAddScat());
        uimenu(ui.cmPh, 'Text', '最も近い散乱体を削除', ...
               'MenuSelectedFcn', @(s,e) menuDelScat());
    end

    %% ---------------- タブ 2: 比較 ----------------
    function buildTabCompare()
        t  = uitab(ui.tabs, 'Title', '比較・評価');
        gg = uigridlayout(t, [2 2]);
        gg.RowHeight = {'1x', 168}; gg.ColumnWidth = {'1x','1x'};
        gg.Padding = [8 8 8 8];
        ui.axDiff = uiaxes(gg);   ui.axDiff.Layout.Row = 1; ui.axDiff.Layout.Column = 1;
        ui.axProf = uiaxes(gg);   ui.axProf.Layout.Row = 1; ui.axProf.Layout.Column = 2;
        ui.tblMet = uitable(gg, 'ColumnName', {'指標', 'A', 'B'}, ...
            'ColumnWidth', {260, 170, 170}, 'Data', cell(0,3));
        ui.tblMet.Layout.Row = 2; ui.tblMet.Layout.Column = [1 2];
        title(ui.axDiff, '差分画像 |A - B|');
        title(ui.axProf, 'ラテラルプロファイル / FWHM');
    end

    %% ---------------- タブ 3: 波面アニメーション ----------------
    function buildTabAnim()
        t  = uitab(ui.tabs, 'Title', '波面アニメーション');
        gg = uigridlayout(t, [2 1]);
        gg.RowHeight = {'1x', 40}; gg.Padding = [8 8 8 8];
        ui.axAnim = uiaxes(gg);
        cg = uigridlayout(gg, [1 6]);
        cg.ColumnWidth = {90, '1x', 120, 70, 90, 150};
        cg.Padding = [0 0 0 0]; cg.ColumnSpacing = 6;
        ui.btnPlay = uibutton(cg, 'Text', '▶ 再生', ...
            'ButtonPushedFcn', @(s,e) onPlayToggle());
        ui.tslider = uislider(cg, 'Limits', [0 1], 'Value', 0, ...
            'MajorTicks', [], 'ValueChangingFcn', @(s,e) onScrub(e.Value), ...
            'ValueChangedFcn', @(s,e) onScrub(s.Value));
        ui.lblT = uilabel(cg, 'Text', 't = 0.00 us');
        uilabel(cg, 'Text', '速度', 'HorizontalAlignment', 'right');
        ui.speed = uidropdown(cg, 'Items', {'x0.25','x0.5','x1','x2'}, 'Value', 'x1');
        ui.txsel = uidropdown(cg, 'Items', {'送信 1'}, 'Value', '送信 1', ...
            'ValueChangedFcn', @(s,e) onTxEventChanged());
    end

    %% ---------------- タブ 4: 遅延 & 整相 ----------------
    function buildTabDelay()
        t  = uitab(ui.tabs, 'Title', '遅延カーブ & 整相');
        gg = uigridlayout(t, [2 1]);
        gg.RowHeight = {'1x', 36}; gg.Padding = [8 8 8 8];
        ag = uigridlayout(gg, [1 3]);
        ag.ColumnWidth = {'1x','1x','1x'}; ag.Padding = [0 0 0 0];
        ui.axRF   = uiaxes(ag);
        ui.axPre  = uiaxes(ag);
        ui.axPost = uiaxes(ag);
        title(ui.axRF, 'RF データ + 遅延曲線');

        cg = uigridlayout(gg, [1 6]);
        cg.ColumnWidth = {190, 60, 80, 60, 80, '1x'};
        cg.Padding = [0 0 0 0]; cg.ColumnSpacing = 6;
        uilabel(cg, 'Text', '再構成点（B モード画像をクリックでも選択）:');
        uilabel(cg, 'Text', 'X [mm]', 'HorizontalAlignment', 'right');
        ui.ptx = uieditfield(cg, 'numeric', 'Value', 0, ...
            'ValueChangedFcn', @(s,e) onPointEdited());
        uilabel(cg, 'Text', 'Z [mm]', 'HorizontalAlignment', 'right');
        ui.ptz = uieditfield(cg, 'numeric', 'Value', 20, ...
            'ValueChangedFcn', @(s,e) onPointEdited());
        ui.lblDelayInfo = uilabel(cg, 'Text', '', 'FontColor', [.3 .3 .3]);
    end

%% ======================================================================
%  UI イベント
%% ======================================================================
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
        ui.preset.Value = 'カスタム';
        drawPhantom();
    end

    function addRow()
        D = ui.tbl.Data;
        ui.tbl.Data = [D; 0, 20, 1];
        ui.preset.Value = 'カスタム';
        drawPhantom();
    end

    function delSelectedRows()
        sel = ui.tbl.Selection;
        D = ui.tbl.Data;
        if isempty(sel) || isempty(D)
            setStatus('削除する行が選択されていません。'); return;
        end
        rows = unique(sel(:,1));
        D(rows, :) = [];
        ui.tbl.Data = D;
        ui.preset.Value = 'カスタム';
        drawPhantom();
    end

    function clearRows()
        ui.tbl.Data = zeros(0, 3);
        ui.preset.Value = 'カスタム';
        drawPhantom();
    end

    function applyPresetAndRefresh()
        applyPreset();
        drawPhantom();
        setStatus('プリセットを適用しました。[1] シミュレーション実行 を押してください。');
    end

    function applyPreset()
        switch ui.preset.Value
            case '単一点散乱体 (PSF)'
                D = [0 20 1];
            case 'ワイヤーファントム'
                xw = -8:4:8;
                zw = 8:6:38;
                [XW, ZW] = meshgrid(xw, zw);
                D = [XW(:), ZW(:), ones(numel(XW), 1)];
            case '無エコー領域 (シスト)'
                rng(0);                                  % 再現性のため固定シード
                n  = 500;
                xb = (rand(n,1)*2 - 1) * 12;
                zb = 5 + rand(n,1) * 33;
                rb = randn(n,1) * 0.35;                  % スペックル背景
                cx = 0; cz = 22; rad = 5;                % 無エコー円
                inCyst = hypot(xb - cx, zb - cz) < rad;
                rb(inCyst) = 0;
                D = [xb, zb, rb];
                D = [D; -8 12 1; 8 32 1];                % 位置基準の点散乱体
            otherwise                                    % カスタム
                return;
        end
        ui.tbl.Data = D;
    end

    function menuAddScat()
        % 匿名関数は生成時の値をキャプチャするため、app の読み出しは
        % 必ずこの入れ子関数の中で行う（座標が [0 20] に固定されるのを防ぐ）。
        addScatAt(app.cursorXZ);
    end

    function menuDelScat()
        delScatAt(app.cursorXZ);
    end

    function captureCursor()
        %CAPTURECURSOR  右クリックでコンテキストメニューが開いた瞬間の座標を保持。
        cp = ui.axPh.CurrentPoint;
        app.cursorXZ = cp(1, 1:2);
    end

    function onPhantomClick(ax, ~)
        cp = ax.CurrentPoint;
        app.cursorXZ = cp(1, 1:2);
        if strcmp(ui.fig.SelectionType, 'alt')
            % ContextMenu を割り当てている場合は通常こちらへは来ない
            % （右クリックはコンテキストメニューに消費される）
            delScatAt(app.cursorXZ);
        elseif ui.clickAdd.Value
            addScatAt(app.cursorXZ);
        end
    end

    function addScatAt(xz)
        x = xz(1); z = xz(2);
        if z <= 0, setStatus('z > 0 の領域をクリックしてください。'); return; end
        ui.tbl.Data = [ui.tbl.Data; x, z, 1];
        ui.preset.Value = 'カスタム';
        drawPhantom();
        setStatus(sprintf('散乱体を追加: (%.2f, %.2f) mm', x, z));
    end

    function delScatAt(xz)
        D = ui.tbl.Data;
        if isempty(D), return; end
        d  = hypot(D(:,1) - xz(1), D(:,2) - xz(2));
        [dm, k] = min(d);
        if dm > 3
            setStatus('近くに散乱体がありません（3 mm 以内が対象）。'); return;
        end
        D(k,:) = [];
        ui.tbl.Data = D;
        ui.preset.Value = 'カスタム';
        drawPhantom();
        setStatus('最も近い散乱体を削除しました。');
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
%  設定の読み出し
%% ======================================================================
    function k = schemeKey()
        switch ui.scheme.Value
            case '集束波 (Focused)',   k = 'focused';
            case '発散波 (Diverging)', k = 'diverging';
            otherwise,                 k = 'plane';
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
%  実行：シミュレーション
%% ======================================================================
    function onSimulate()
        cfg = readCfg();
        setStatus('シミュレーション実行中...'); drawnow;
        try
            app.S = sim_engine(cfg);
        catch ME
            uialert(ui.fig, ME.message, 'シミュレーション失敗');
            setStatus(['シミュレーション失敗: ' ME.message]);
            return;
        end
        app.imgA = []; app.imgB = [];
        updateTxEventList();
        drawPhantom();
        setupAnimation();
        updateDelayViews();
        setStatus(sprintf(['RF 生成完了 | backend = %s | 送信イベント %d | ' ...
            'RF %d x %d | fs = %.1f MHz | %.2f 秒'], app.S.backend, numel(app.S.tx), ...
            size(app.S.RF,1), size(app.S.RF,2), app.S.fs/1e6, app.S.elapsed));
    end

    function updateTxEventList()
        n = numel(app.S.tx);
        items = arrayfun(@(k) sprintf('送信 %d', k), 1:n, 'UniformOutput', false);
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
%  実行：ビームフォーミングと比較
%% ======================================================================
    function onBeamform()
        if isempty(app.S)
            onSimulate();
            if isempty(app.S), return; end
        end
        [gx, gz] = reconGrid();
        app.gx = gx; app.gz = gz;

        % f 値・受信アポダイゼーションは送信後に変えられるので毎回反映する
        for k = 1:numel(app.S.tx)
            app.S.tx(k).fnumber = ui.fnum.Value;
            app.S.tx(k).rx_apod = ui.rxapod.Value;
        end

        setStatus('ビームフォーミング中...'); drawnow;
        try
            [app.imgA, app.msA, app.nameA] = runAlgorithm(ui.algA.Value, gx, gz);
        catch ME
            uialert(ui.fig, ME.message, 'アルゴリズム A の実行に失敗');
            setStatus(['A 失敗: ' ME.message]); return;
        end

        app.imgB = []; app.msB = NaN; app.nameB = '';
        if ~strcmp(ui.algB.Value, 'なし')
            try
                [app.imgB, app.msB, app.nameB] = runAlgorithm(ui.algB.Value, gx, gz);
            catch ME
                uialert(ui.fig, ME.message, 'アルゴリズム B の実行に失敗');
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
        setStatus(['ビームフォーミング完了   ' msg]);
    end

    function [img, ms, name] = runAlgorithm(key, gx, gz)
        [fh, name] = algoHandle(key);
        S = app.S;

        % JIT ウォームアップ（極小グリッドで 1 回捨て実行 → 計測値を汚さない）
        try
            fh(S.RF(:,:,1), S.tx(1), S.rx_pos, gx(1:min(4,end)), gz(1:min(4,end)), S.c, S.fs);
        catch
            % ウォームアップ失敗は無視して本実行のエラーで報告する
        end

        tic;
        img = runAllTx(fh, gx, gz);
        ms  = toc * 1000;

        expected = [numel(gz) numel(gx)];
        if ~isequal(size(img), expected)
            error('main_gui:imgSize', ...
                ['%s が返した画像サイズ %s が想定 %s と一致しません。\n' ...
                 '統一インターフェースでは [numel(grid_z) x numel(grid_x)] の\n' ...
                 '線形エンベロープを返す必要があります。'], name, ...
                 mat2str(size(img)), mat2str(expected));
        end
        if ~isreal(img) || any(~isfinite(img(:)))
            img = abs(img);
            img(~isfinite(img)) = 0;
        end
    end

    function img = runAllTx(fh, gx, gz)
        %RUNALLTX  マルチフォーカス送信は深度ゾーンごとに画像を切り替えて合成する。
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
            tmp = fh(S.RF(:,:,ord(k)), S.tx(ord(k)), S.rx_pos, gx, gz(sub), S.c, S.fs);
            img(sub, :) = tmp;
        end
    end

    function [fh, name] = algoHandle(key)
        switch key
            case '参照DAS (das_reference)'
                fh = @das_reference; name = 'das_reference';
            case '自作DAS (Custom)'
                fname = strtrim(ui.customfn.Value);
                if isempty(which(fname))
                    error('main_gui:noCustom', ...
                        ['自作関数 "%s" が見つかりません。\n' ...
                         'das_custom_template.m をコピーして関数名を合わせ、\n' ...
                         'MATLAB のパス（またはカレントフォルダ）に置いてください。'], fname);
                end
                fh = str2func(fname); name = fname;
            case 'MUST das()'
                if isempty(which('das'))
                    error('main_gui:noMUST', ...
                        ['MUST の das() が見つかりません。\n' ...
                         'MUST (Matlab UltraSound Toolbox) をインストールし、\n' ...
                         'addpath でパスを通してください。']);
                end
                fh = @mustDAS; name = 'MUST das()';
            otherwise
                error('main_gui:algo', '未知のアルゴリズム "%s" です。', key);
        end
    end

    function [img, delays] = mustDAS(rf, tx, rx_pos, gx, gz, c, fs)
        %MUSTDAS  MUST 組み込みの das() を統一インターフェースに合わせるアダプタ。
        %  検証済み構文: BFSIG = DAS(SIG, X, Z, DELAYS, PARAM)
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
            % MUST の das() は遅延を返さないため、参照実装の遅延で代用する
            [~, delays] = das_reference(rf, tx, rx_pos, gx, gz, c, fs);
        end
    end

%% ======================================================================
%  描画
%% ======================================================================
    function drawPhantom()
        ax = ui.axPh;
        D  = ui.tbl.Data;
        cla(ax, 'reset');
        hold(ax, 'on');

        % 素子列
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

        xh = ui.xhalf.Value;
        xlim(ax, [-max(xh, max(abs(ex))*1.1), max(xh, max(abs(ex))*1.1)]);
        ylim(ax, [-2, ui.zmax.Value]);
        set(ax, 'YDir', 'reverse');
        xlabel(ax, 'x [mm]'); ylabel(ax, 'z [mm]');
        title(ax, sprintf('ファントム配置（散乱体 %d 個）', size(D,1)));
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
            title(ax, '（未計算）');
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
            Bn = B / nA;                       % ★ A の最大値で共通正規化する
            Dm = abs(An - Bn);
            mse = mean((An(:) - Bn(:)).^2);
            imagesc(ui.axDiff, app.gx*1e3, app.gz*1e3, 20*log10(Dm + 1e-12));
            set(ui.axDiff, 'CLim', [-80 0], 'YDir', 'reverse');
            colormap(ui.axDiff, parula(256));
            colorbar(ui.axDiff);
            axis(ui.axDiff, 'image');
            xlabel(ui.axDiff, 'x [mm]'); ylabel(ui.axDiff, 'z [mm]');
            title(ui.axDiff, sprintf('|A - B| [dB]  (A の最大値で共通正規化)  MSE = %.3e', mse));
            [fwB, profB, pkB] = lateralMetrics(Bn, app.gx, app.gz, pkA(2));
        else
            text(ui.axDiff, 0.5, 0.5, 'アルゴリズム B を選ぶと差分が表示されます', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
            axis(ui.axDiff, 'off');
            fwB = NaN; profB = []; pkB = [NaN NaN]; mse = NaN; Dm = [];
        end

        % ---- ラテラルプロファイル ----
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
        xlabel(ui.axProf, 'x [mm]'); ylabel(ui.axProf, '正規化振幅 [dB]');
        title(ui.axProf, sprintf('深さ z = %.2f mm のラテラルプロファイル', pkA(2)*1e3));
        legend(ui.axProf, 'Location', 'southwest', 'FontSize', 9, 'Interpreter', 'none');

        % ---- 指標テーブル ----
        rows(end+1,:) = {'実行時間 [ms]',            num2str(app.msA, '%.1f'),  numOrDash(app.msB, '%.1f')};
        rows(end+1,:) = {'ピーク位置 x [mm]',        num2str(pkA(1)*1e3, '%.3f'), numOrDash(pkB(1)*1e3, '%.3f')};
        rows(end+1,:) = {'ピーク位置 z [mm]',        num2str(pkA(2)*1e3, '%.3f'), numOrDash(pkB(2)*1e3, '%.3f')};
        rows(end+1,:) = {'ラテラル FWHM (-6dB) [mm]',num2str(fwA*1e3, '%.4f'),   numOrDash(fwB*1e3, '%.4f')};
        rows(end+1,:) = {'MSE (A 基準正規化)',       '-',                        numOrDash(mse, '%.4e')};
        if ~isempty(Dm)
            rows(end+1,:) = {'最大絶対誤差 (A 基準)', '-', num2str(max(Dm(:)), '%.4e')};
        end
        rows(end+1,:) = {'アルゴリズム名', app.nameA, dashIfEmpty(app.nameB)};
        rows(end+1,:) = {'バックエンド',   app.S.backend, ''};
        ui.tblMet.Data = rows;
    end

    function [fw, prof, pk] = lateralMetrics(img, gx, gz, forceZ)
        %LATERALMETRICS  ピーク行のラテラルプロファイルと -6dB 全幅（線形値で評価）。
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
        %FWHM  線形エンベロープの半値（= -6 dB）全幅。交点は線形補間で求める。
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
%  波面アニメーション
%% ======================================================================
    function setupAnimation()
        if isempty(app.S), return; end
        xh = ui.xhalf.Value * 1e-3;
        xr = [-max(xh, max(abs(app.S.rx_pos))*1.05), max(xh, max(abs(app.S.rx_pos))*1.05)];
        zr = [0, ui.zmax.Value*1e-3];
        tmax = wave_animator('setup_propagation', ui.axAnim, app.S, currentTx(), xr, zr);
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
            setStatus('先に [1] シミュレーション実行 を押してください。'); return;
        end
        app.playing = true;
        ui.btnPlay.Text = '■ 停止';
        spd  = str2double(strrep(ui.speed.Value, 'x', ''));
        tmax = ui.tslider.Limits(2);
        dt   = tmax / 200 * spd;
        t    = app.animT;
        while app.playing && isvalid(ui.fig)
            t = t + dt;
            if t > tmax, t = 0; end
            app.animT = t;
            ui.tslider.Value = t;
            wave_animator('draw_propagation', ui.axAnim, t);
            ui.lblT.Text = sprintf('t = %.2f us', t*1e6);
            drawnow limitrate;
        end
        if isvalid(ui.fig)
            app.playing = false;
            ui.btnPlay.Text = '▶ 再生';
        end
    end

%% ======================================================================
%  遅延曲線 / 整相前後
%% ======================================================================
    function updateDelayViews()
        if isempty(app.S), return; end
        k = currentTx();
        note = '';
        try
            [fh, name] = algoHandle(ui.algA.Value);
        catch
            fh = @das_reference; name = 'das_reference (代替)';
        end
        try
            % 1 点だけのグリッドで呼ぶと、ビームフォーマが実際に加算している
            % サンプル位置がそのまま得られる（可視化と計算の乖離が起きない）
            [~, tau] = fh(app.S.RF(:,:,k), app.S.tx(k), app.S.rx_pos, ...
                          app.pt(1), app.pt(2), app.S.c, app.S.fs);
        catch ME
            [~, tau] = das_reference(app.S.RF(:,:,k), app.S.tx(k), app.S.rx_pos, ...
                          app.pt(1), app.pt(2), app.S.c, app.S.fs);
            name = 'das_reference (代替)';
            note = ['  ※' ME.message];
        end
        tau = tau(1, :);
        wave_animator('delay_curve', ui.axRF, app.S, k, tau, app.pt);
        wave_animator('alignment', ui.axPre, ui.axPost, app.S, k, tau, app.pt);
        ui.lblDelayInfo.Text = sprintf('遅延の出所: %s   有効素子 %d / %d%s', ...
            name, sum(isfinite(tau)), numel(tau), note);
    end

%% ======================================================================
%  小物
%% ======================================================================
    function setStatus(msg)
        ui.status.Text = ['  ' msg];
        drawnow limitrate;
    end

    function s = mustStatusText()
        if isempty(which('simus'))
            s = '未検出（モックで動作します）';
        else
            s = ['検出済み (' which('simus') ')'];
        end
    end

    function v = parseList(str)
        parts = strsplit(strtrim(str), {',', ' ', ';'});
        v = str2double(parts);
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
