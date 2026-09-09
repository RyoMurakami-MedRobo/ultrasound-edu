function varargout = wave_animator(action, varargin)
%WAVE_ANIMATOR  波面伝搬・遅延曲線・整相前後比較の可視化ロジック。
%
%   MAIN_GUI から呼ばれる描画専用モジュール。状態は対象 axes の UserData に
%   保持するため、呼び出し側は axes ハンドルだけを持っていればよい。
%
%   ---- 波面アニメーション ---------------------------------------------
%   tmax = WAVE_ANIMATOR('setup_propagation', ax, S, txIdx, xlim_m, zlim_m)
%       送信波面と散乱エコーの場を計算する準備を行い、アニメーションの
%       最終時刻 [s] を返す。S は SIM_ENGINE の出力。
%   WAVE_ANIMATOR('draw_propagation', ax, t)
%       時刻 t [s] のスナップショットを描画する（CData の更新のみ）。
%
%   ---- 遅延曲線オーバーレイ -------------------------------------------
%   WAVE_ANIMATOR('delay_curve', ax, S, txIdx, tau, pt)
%       RF データを画像表示し、その上に「再構成点 pt に対して各素子の
%       どのサンプルが加算されるか」を遅延曲線として重ね描きする。
%       tau は [1 x Nel] の往復合計時間 [s]（DAS 実装の delays 出力の 1 行）。
%
%   ---- 整相前後の比較 -------------------------------------------------
%   WAVE_ANIMATOR('alignment', axPre, axPost, S, txIdx, tau, pt)
%       整相前のばらついた波形束（左）と、遅延補正で同相に揃った波形束＋
%       積算結果（右）を並べて描画する。
%
%   See also SIM_ENGINE, DAS_REFERENCE, MAIN_GUI.

switch lower(action)
    case 'setup_propagation'
        varargout{1} = setupPropagation(varargin{:});
    case 'draw_propagation'
        drawPropagation(varargin{:});
    case 'delay_curve'
        drawDelayCurve(varargin{:});
    case 'alignment'
        drawAlignment(varargin{:});
    otherwise
        error('wave_animator:action', '未知のアクション "%s" です。', action);
end
end

%% ======================================================================
function tmax = setupPropagation(ax, S, txIdx, xrange, zrange)
%SETUPPROPAGATION  波面／エコーの到達時刻マップを事前計算して描画枠を作る。
tx = S.tx(txIdx);
c  = S.c;

NX = 241; NZ = 241;
xv = linspace(xrange(1), xrange(2), NX);
zv = linspace(max(zrange(1), 0), zrange(2), NZ);
[X, Z] = meshgrid(xv, zv);

% 送信波面が各点に到達する時刻（DAS_REFERENCE と同一モデル）
TAU = txArrivalTime(X, Z, tx, c);
TAU = reshape(TAU, size(X));

% 散乱体ごとのエコー：到達時刻 + 散乱体からの距離 / c
%   スペックルファントムのように散乱体が多い場合は描画が破綻するので、
%   反射強度の大きい順に MAXSCAT 個だけをアニメーション対象にする。
MAXSCAT = 30;
xs = S.scat.x; zs = S.scat.z; rc = S.scat.rc;
keep = rc ~= 0;
xs = xs(keep); zs = zs(keep); rc = rc(keep);
nDropped = 0;
if numel(xs) > MAXSCAT
    [~, ord] = sort(abs(rc), 'descend');
    ord = ord(1:MAXSCAT);
    nDropped = numel(xs) - MAXSCAT;
    xs = xs(ord); zs = zs(ord); rc = rc(ord);
end
Ns = numel(xs);
tauS = zeros(1, Ns);
Rmat = zeros(NZ, NX, Ns);          % 散乱体からの距離を前計算（毎フレーム再利用）
for k = 1:Ns
    tauS(k) = txArrivalTime(xs(k), zs(k), tx, c);
    Rmat(:,:,k) = hypot(X - xs(k), Z - zs(k));
end

% パルス幅（可視化用の波面の厚み）
sigma = sqrt(2*log(2)) / (pi * S.probe.fc * max(S.probe.bandwidth,1)/100);
sigma = max(sigma, 0.3/S.probe.fc);

cla(ax, 'reset');
hold(ax, 'on');
st = struct();
st.X = X; st.Z = Z; st.TAU = TAU; st.c = c; st.sigma = sigma;
st.xs = xs; st.zs = zs; st.rc = rc; st.tauS = tauS;
st.Rmat = Rmat; st.nDropped = nDropped;
st.xv = xv; st.zv = zv;

% 場を表示する画像（赤 = 送信波面、青 = 散乱エコー）
st.hImg = imagesc(ax, xv*1e3, zv*1e3, zeros(NZ, NX));
colormap(ax, blueWhiteRed());
set(ax, 'CLim', [-1 1]);

% 素子・散乱体・焦点マーカー
plot(ax, tx.elem_x*1e3, tx.elem_z*1e3, 's', 'MarkerSize', 4, ...
     'MarkerEdgeColor', [.25 .25 .25], 'MarkerFaceColor', [.75 .75 .75], ...
     'HitTest','off');
act = isfinite(tx.delays) & tx.apod ~= 0;
plot(ax, tx.elem_x(act)*1e3, tx.elem_z(act)*1e3, 's', 'MarkerSize', 4, ...
     'MarkerEdgeColor', [.6 0 0], 'MarkerFaceColor', [1 .4 .4], 'HitTest','off');
st.hHit = plot(ax, nan, nan, 's', 'MarkerSize', 7, 'MarkerEdgeColor', [0 0 .6], ...
     'MarkerFaceColor', [.3 .6 1], 'HitTest','off');
if Ns > 0
    plot(ax, xs*1e3, zs*1e3, 'o', 'MarkerSize', 6, 'LineWidth', 1.2, ...
         'MarkerEdgeColor', [0 .5 0], 'MarkerFaceColor', 'none', 'HitTest','off');
end
if all(isfinite(tx.src_xz)) && ~strcmpi(tx.scheme,'plane')
    plot(ax, tx.src_xz(1)*1e3, tx.src_xz(2)*1e3, 'p', 'MarkerSize', 12, ...
         'MarkerEdgeColor', [.5 .3 0], 'MarkerFaceColor', [1 .8 .2], 'HitTest','off');
end

st.hTitle = title(ax, '');
axis(ax, 'image');
set(ax, 'YDir', 'reverse');
xlim(ax, xrange*1e3);
ylim(ax, [min(zv) max(zv)]*1e3);
xlabel(ax, 'x [mm]'); ylabel(ax, 'z [mm]');
grid(ax, 'off');
hold(ax, 'off');

% アニメーション終了時刻：最も遠い素子までエコーが戻りきる時刻
if Ns > 0
    dmax = max(hypot(xs(:) - tx.elem_x(:).', zs(:) - tx.elem_z(:).'), [], 2);
    tmax = max(tauS(:) + dmax/c) + 6*sigma;
else
    tmax = max(TAU(:)) + 6*sigma;
end
tmax = max(tmax, 6*sigma);
st.tmax = tmax;
st.elem_x = tx.elem_x; st.elem_z = tx.elem_z;
ax.UserData = st;
end

%% ======================================================================
function drawPropagation(ax, t)
%DRAWPROPAGATION  時刻 t のスナップショットに更新する。
st = ax.UserData;
if isempty(st) || ~isfield(st, 'hImg') || ~isvalid(st.hImg), return; end

sig = st.sigma;

% 送信波面（正 = 赤）
F = exp(-0.5*((st.TAU - t)/sig).^2);

% 散乱エコー（負 = 青）: 散乱体が照射された後に同心円状に広がる
%   重なった円弧が加算で飽和して送信波面を覆い隠さないよう、和ではなく
%   最大値で合成する。
Fe   = zeros(size(F));
hitx = [];
for k = 1:numel(st.xs)
    dt = t - st.tauS(k);
    if dt <= 0, continue; end
    amp = min(abs(st.rc(k)), 1) * 0.9;
    Fe  = max(Fe, amp * exp(-0.5*((st.Rmat(:,:,k)/st.c - dt)/sig).^2));

    % エコーが到達した受信素子をハイライト
    re = hypot(st.elem_x - st.xs(k), st.elem_z - st.zs(k));
    hitx = [hitx, st.elem_x(abs(re/st.c - dt) < 2*sig)]; %#ok<AGROW>
end

st.hImg.CData = max(min(F - Fe, 1), -1);
if isempty(hitx)
    set(st.hHit, 'XData', nan, 'YData', nan);
else
    set(st.hHit, 'XData', hitx*1e3, 'YData', zeros(size(hitx)));
end
if st.nDropped > 0
    extra = sprintf('  [表示は強い散乱体 %d 個のみ / 他 %d 個は非表示]', ...
                    numel(st.xs), st.nDropped);
else
    extra = '';
end
st.hTitle.String = sprintf('t = %.2f \\mus   （赤: 送信波面 / 青: 散乱エコー）%s', ...
                           t*1e6, extra);
end

%% ======================================================================
function drawDelayCurve(ax, S, txIdx, tau, pt)
%DRAWDELAYCURVE  RF 画像の上に、加算対象サンプルを示す遅延曲線を重ねる。
RF = S.RF(:,:,txIdx);
Nel = size(RF, 2);
tus = (S.t0 + (0:size(RF,1)-1)/S.fs) * 1e6;

cla(ax, 'reset');
% 深さ方向のゲイン補償をかけて浅部の飽和を防ぐ（表示のみ）
disp_rf = RF ./ max(max(abs(RF), [], 2), 1e-12).^0.7;
imagesc(ax, 1:Nel, tus, disp_rf);
colormap(ax, gray(256));
set(ax, 'CLim', [-1 1]*max(abs(disp_rf(:)))*0.6 + [-eps eps]);
hold(ax, 'on');

tauUs = tau(:).' * 1e6;
ok = isfinite(tauUs);
plot(ax, find(ok), tauUs(ok), '-', 'Color', [1 .25 .1], 'LineWidth', 2);
plot(ax, find(ok), tauUs(ok), '.', 'Color', [1 .8 0], 'MarkerSize', 8);
if any(~ok)
    yl = mean(tauUs(ok), 'omitnan');
    plot(ax, find(~ok), repmat(yl, 1, sum(~ok)), 'x', 'Color', [.4 .4 .4], ...
         'MarkerSize', 5);
end
hold(ax, 'off');
set(ax, 'YDir', 'normal');
xlabel(ax, '受信素子番号');
ylabel(ax, '時間 [\mus]');
if any(ok)
    span = max(max(tauUs(ok)) - min(tauUs(ok)), 0.5);
    ylim(ax, [min(tauUs(ok)) - 0.8*span - 1, max(tauUs(ok)) + 0.8*span + 1]);
end
xlim(ax, [0.5 Nel+0.5]);
title(ax, sprintf('遅延曲線  再構成点 (x, z) = (%.2f, %.2f) mm   有効素子 %d/%d', ...
      pt(1)*1e3, pt(2)*1e3, sum(ok), Nel));
end

%% ======================================================================
function drawAlignment(axPre, axPost, S, txIdx, tau, pt)
%DRAWALIGNMENT  整相前（左）と整相後＋積算（右）の波形束を描画する。
RF  = S.RF(:,:,txIdx);
Nel = size(RF, 2);
fc  = S.probe.fc;

tau = tau(:).';
ok  = isfinite(tau);
if ~any(ok)
    cla(axPre,'reset');  text(axPre, .5,.5,'有効な素子がありません','Units','normalized','HorizontalAlignment','center');
    cla(axPost,'reset'); text(axPost,.5,.5,'有効な素子がありません','Units','normalized','HorizontalAlignment','center');
    return;
end

trel = linspace(-5/fc, 5/fc, 401);          % 中心周波数の ±5 周期分
tc   = mean(tau(ok));

pre  = zeros(numel(trel), Nel);
post = zeros(numel(trel), Nel);
for e = 1:Nel
    pre(:,e) = sampleRF(RF(:,e), S.t0, S.fs, tc + trel);
    if ok(e)
        post(:,e) = sampleRF(RF(:,e), S.t0, S.fs, tau(e) + trel);
    end
end
scale = max([max(abs(pre(:))), max(abs(post(:))), eps]);
gain  = 3.0;

drawBundle(axPre,  trel*1e6, pre,  ok, scale, gain, tau, tc);
title(axPre, sprintf('整相前（共通時間窓 t = %.2f \\mus 中心）', tc*1e6));

drawBundle(axPost, trel*1e6, post, ok, scale, gain, [], []);
% 積算結果を波形束の下に重ねる
sumTrace = sum(post(:, ok), 2);
sScale = max(abs(sumTrace)) + eps;
hold(axPost, 'on');
yBase = -0.10*Nel - 2;
plot(axPost, trel*1e6, yBase + 0.9*(0.06*Nel+2)*sumTrace/sScale, ...
     '-', 'Color', [.85 .1 .1], 'LineWidth', 1.8);
plot(axPost, [trel(1) trel(end)]*1e6, [yBase yBase], ':', 'Color', [.6 .6 .6]);
text(axPost, trel(1)*1e6, yBase, sprintf('  \\Sigma  (peak %.3g)', max(abs(sumTrace))), ...
     'Color', [.85 .1 .1], 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
ylim(axPost, [yBase - (0.06*Nel+3), Nel + gain + 1]);
hold(axPost, 'off');
title(axPost, sprintf('整相後 + 積算   点 (%.2f, %.2f) mm', pt(1)*1e3, pt(2)*1e3));
end

%% ----------------------------------------------------------------------
function drawBundle(ax, tus, W, ok, scale, gain, tau, tc)
cla(ax, 'reset');
hold(ax, 'on');
Nel = size(W, 2);
for e = 1:Nel
    if ok(e)
        col = [0 .35 .75];
    else
        col = [.75 .75 .75];
    end
    plot(ax, tus, e + gain*W(:,e)/scale, '-', 'Color', col, 'LineWidth', 0.6);
end
if ~isempty(tau)
    % 整相前は各素子の「加算すべき時刻」がばらついていることを示す
    d = (tau - tc) * 1e6;
    plot(ax, d(ok), find(ok), '.', 'Color', [1 .3 0], 'MarkerSize', 9);
    plot(ax, d(ok), find(ok), '-', 'Color', [1 .3 0], 'LineWidth', 1.2);
else
    plot(ax, [0 0], [0 Nel+gain+1], '-', 'Color', [1 .3 0], 'LineWidth', 1.2);
end
hold(ax, 'off');
xlabel(ax, '相対時間 [\mus]');
ylabel(ax, '受信素子番号');
xlim(ax, [tus(1) tus(end)]);
ylim(ax, [0, Nel + gain + 1]);
grid(ax, 'on');
end

%% ----------------------------------------------------------------------
function v = sampleRF(col, t0, fs, times)
%SAMPLERF  RF の 1 チャンネルを任意時刻で線形補間する（範囲外は 0）。
n   = numel(col);
pos = (times(:) - t0)*fs + 1;
i0  = floor(pos);
fr  = pos - i0;
ok  = i0 >= 1 & i0 < n;
i0c = min(max(i0,1), n-1);
v   = col(i0c).*(1-fr) + col(i0c+1).*fr;
v(~ok) = 0;
end

%% ----------------------------------------------------------------------
function tau = txArrivalTime(X, Z, tx, c)
%TXARRIVALTIME  送信波面の到達時刻 [s]。
%  ※ DAS_REFERENCE の同名ローカル関数と同一モデル。可視化モジュールを
%    ビームフォーマから独立させるため意図的に複製している。
%    （「ビームフォーマが仮定している波面」をそのまま描くことが目的）
ex = tx.elem_x; ez = tx.elem_z; td = tx.delays; ap = tx.apod;
if isempty(ez), ez = zeros(size(ex)); end
act = isfinite(td) & ap(:).' ~= 0;
ex = ex(act); ez = ez(act); td = td(act);
xg = X(:); zg = Z(:);

src = tx.src_xz;
isFocused = strcmpi(tx.scheme,'focused') && numel(src)==2 && all(isfinite(src)) && ...
            src(2) > 0 && ~tx.single_element;
if isFocused
    Rf   = hypot(ex - src(1), ez - src(2));
    Tf   = mean(td + Rf/c);
    D    = hypot(src(1), src(2));
    u    = [src(1) src(2)]/D;
    proj = xg*u(1) + zg*u(2);
    dF   = hypot(xg - src(1), zg - src(2));
    sgn  = ones(numel(xg),1); sgn(proj <= D) = -1;
    tau  = Tf + sgn.*dF/c;
else
    tau = inf(numel(xg), 1);
    for e = 1:numel(ex)
        tau = min(tau, td(e) + hypot(xg - ex(e), zg - ez(e))/c);
    end
end
end

%% ----------------------------------------------------------------------
function cmap = blueWhiteRed()
%BLUEWHITERED  青-白-赤のダイバージングカラーマップ（ツールボックス不要）。
n = 128;
u = linspace(0, 1, n).';
lower_ = [u, u, ones(n,1)];                 % 青 -> 白
upper_ = [ones(n,1), flipud(u), flipud(u)]; % 白 -> 赤
cmap = [lower_; upper_];
end
