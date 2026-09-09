function [bmode_img, delays] = das_reference(rf_data, tx_info, rx_pos, grid_x, grid_z, sound_speed, fs)
%DAS_REFERENCE  教育用の標準 Delay-And-Sum ビームフォーマ（参照実装）。
%
%   [bmode_img, delays] = DAS_REFERENCE(rf_data, tx_info, rx_pos, ...
%                                       grid_x, grid_z, sound_speed, fs)
%
%   本ツールの「自作 DAS 検証フレームワーク」が要求する統一インターフェース。
%   自作アルゴリズム（DAS_CUSTOM_TEMPLATE 参照）もこの引数・戻り値に従うこと。
%
%   ---- 入力 -----------------------------------------------------------
%   rf_data     [Nt x Nel] RF データ。**時間が第1次元、列が受信素子**。
%               3 次元配列が渡された場合は第1送信イベントのみを使用する。
%   tx_info     送信条件 struct。以下を参照する:
%                 .delays  [1 x Nel] 送信遅延 [s]（非励振素子は NaN）
%                 .apod    [1 x Nel] 送信アポダイゼーション（0 で非励振）
%                 .elem_x  [1 x Nel] 送信素子 x 座標 [m]
%                 .elem_z  [1 x Nel] 送信素子 z 座標 [m]
%                 .t0      RF 第1サンプルの時刻 [s]
%                 .fnumber 受信 f 値（0 または非有限で全開口）
%                 .rx_apod 'rect' | 'hann'
%   rx_pos      [1 x Nel] または [2 x Nel] 受信素子座標 [m]。
%               1 行なら x 座標のみとみなし z = 0 とする。受信ジオメトリは
%               この引数を正とする（tx_info.elem_* は送信側の定義）。
%   grid_x,grid_z 再構成グリッド [m]。ベクトルを渡すと meshgrid され、
%               画像サイズは [numel(grid_z) x numel(grid_x)] になる。
%               同サイズの行列を渡した場合はそのまま点群として扱う。
%   sound_speed 音速 [m/s]
%   fs          サンプリング周波数 [Hz]
%
%   ---- 出力 -----------------------------------------------------------
%   bmode_img   [Nz x Nx] **線形エンベロープ**（対数圧縮前）。
%               表示側で 20*log10(img/max(img)) して用いる。
%               MSE や FWHM といった指標は線形値で評価する方が正しいため、
%               本 I/F では線形エンベロープを返す契約とする。
%   delays      [Nz*Nx x Nel] 各画素・各素子の **往復合計時間 [s]**。
%               f 値や非励振によって加算対象外の要素は NaN。
%               nargout < 2 のときは計算しない（ベンチマーク時の負荷を避ける）。
%
%   ---- アルゴリズム ---------------------------------------------------
%   1) 送信到達時間 tau_tx(p)
%      ・平面波 / 発散波 / 単一素子 : 初到達（Huygens）モデル
%          tau_tx(p) = min_e ( txdelay(e) + |p - e| / c )
%        これらの送信では閉形式（x*sin+z*cos や仮想音源モデル）と厳密に
%        一致するため、分岐なしに扱える。
%      ・集束波 : 仮想音源モデル
%          tau_tx(p) = T_F + s * |p - F| / c   (s = -1 焦点手前 / +1 以遠)
%        初到達モデルは焦点以遠で開口端の弱い端部波を拾ってしまうため。
%   2) 受信到達時間  tau_rx(p,e) = |p - e_rx| / c
%   3) 合計時間 tau = tau_tx + tau_rx をサンプル位置へ換算し線形補間
%   4) f 値による受信開口制限とアポダイゼーションを掛けて素子方向に総和
%   5) ヒルベルト変換で深さ方向にエンベロープ検波
%
%   See also SIM_ENGINE, DAS_CUSTOM_TEMPLATE, WAVE_ANIMATOR.

narginchk(7, 7);

%% ---- 入力整形 -------------------------------------------------------
if ndims(rf_data) > 2 %#ok<ISMAT>
    rf_data = rf_data(:,:,1);
end
[Nt, Nel] = size(rf_data);

[rxx, rxz] = normalizeRxPos(rx_pos, Nel);
[XI, ZI, imgSize] = normalizeGrid(grid_x, grid_z);
Npix = numel(XI);

c  = sound_speed;
t0 = getdef(tx_info, 't0', 0);
fnum    = getdef(tx_info, 'fnumber', 0);
rxApod  = lower(getdef(tx_info, 'rx_apod', 'rect'));
if ~isfinite(fnum) || fnum <= 0, fnum = 0; end   % 0 = 全開口

%% ---- 1) 送信到達時間（初到達 / Huygens モデル）----------------------
tau_tx = txArrivalTime(XI, ZI, tx_info, c);

%% ---- 準備 ----------------------------------------------------------
bf = zeros(Npix, 1);
if nargout > 1
    delays = nan(Npix, Nel);
end

% f 値から決まる受信開口の半幅（深さに比例）
if fnum > 0
    halfAp = ZI(:) ./ (2*fnum);
else
    halfAp = inf(Npix, 1);
end

xg = XI(:); zg = ZI(:);

%% ---- 2)-4) 素子ごとの遅延加算 --------------------------------------
for e = 1:Nel
    dx  = xg - rxx(e);
    rrx = hypot(dx, zg - rxz(e));
    tau = tau_tx + rrx / c;                       % 往復合計時間 [s]

    % 受信開口制限
    inAp = abs(dx) <= halfAp;
    switch rxApod
        case 'hann'
            u = zeros(Npix,1);
            u(inAp) = dx(inAp) ./ max(halfAp(inAp), eps);   % -1..1
            w = inAp .* (0.5 * (1 + cos(pi*u)));
        otherwise
            w = double(inAp);
    end

    % 線形補間で RF をサンプリング
    pos  = (tau - t0) * fs + 1;
    i0   = floor(pos);
    frac = pos - i0;
    ok   = w > 0 & i0 >= 1 & i0 < Nt;
    i0c  = min(max(i0, 1), Nt-1);
    val  = rf_data(i0c, e) .* (1-frac) + rf_data(i0c+1, e) .* frac;
    val(~ok) = 0;

    bf = bf + w .* val;

    if nargout > 1
        d = tau;
        d(~ok) = NaN;
        delays(:, e) = d;
    end
end

%% ---- 5) エンベロープ検波 -------------------------------------------
bmode_img = envelopeZ(reshape(bf, imgSize));
end % das_reference

%% ======================================================================
function env = envelopeZ(x)
%ENVELOPEZ  第1次元（深さ方向）のエンベロープ検波。
if size(x,1) < 4
    env = abs(x);                       % 1 点／極小グリッドの縮退ケース
elseif ~isempty(which('hilbert'))
    env = abs(hilbert(x));              % Signal Processing Toolbox
else
    % ツールボックス無しでも動くよう FFT で解析信号を作る
    n = size(x,1);
    X = fft(x, n, 1);
    h = zeros(n,1);
    h(1) = 1;
    if mod(n,2) == 0
        h(n/2+1) = 1;  h(2:n/2) = 2;
    else
        h(2:(n+1)/2) = 2;
    end
    env = abs(ifft(X .* h, n, 1));
end
end

%% ======================================================================
function tau = txArrivalTime(XI, ZI, tx_info, c)
%TXARRIVALTIME  送信波面が各画素に到達する時刻（初到達）を返す [Npix x 1]。
ex = getdef(tx_info, 'elem_x', []);
ez = getdef(tx_info, 'elem_z', []);
td = getdef(tx_info, 'delays', []);
ap = getdef(tx_info, 'apod', []);

if isempty(ex) || isempty(td)
    error('das_reference:txinfo', 'tx_info に elem_x と delays が必要です。');
end
if isempty(ez), ez = zeros(size(ex)); end
if isempty(ap), ap = ones(size(ex)); end

act = isfinite(td) & ap(:).' ~= 0;
if ~any(act)
    error('das_reference:noActiveTx', '送信に有効な素子がありません。');
end
ex = ex(act); ez = ez(act); td = td(act);
xg = XI(:); zg = ZI(:);

src    = getdef(tx_info, 'src_xz', [NaN NaN]);
scheme = getdef(tx_info, 'scheme', '');
single = getdef(tx_info, 'single_element', false);
isFocused = strcmpi(scheme,'focused') && numel(src)==2 && all(isfinite(src)) && ...
            src(2) > 0 && ~single;

if isFocused
    % ---- 集束送信: 仮想音源モデル -----------------------------------
    % 初到達（min）モデルは焦点より手前では厳密だが、焦点以遠では開口端
    % から出る弱い端部波を拾ってしまい、主エネルギーの到達時刻より早く
    % なる。そこで集束送信だけは焦点を通過する球面波として扱う。
    %   tau_tx(p) = T_F + s * |p - F| / c   ( s = -1 焦点手前 / +1 焦点以遠 )
    % T_F（波面が焦点に収束する時刻）は遅延則の定義から全素子で一定なので
    % 遅延そのものから求められる（遅延の正規化方法に依存しない）。
    Rf   = hypot(ex - src(1), ez - src(2));
    Tf   = mean(td + Rf / c);
    D    = hypot(src(1), src(2));
    u    = [src(1) src(2)] / D;                 % ビーム軸の単位ベクトル
    proj = xg*u(1) + zg*u(2);                   % ビーム軸方向の投影距離
    dF   = hypot(xg - src(1), zg - src(2));
    sgn  = ones(numel(xg), 1);
    sgn(proj <= D) = -1;
    tau  = Tf + sgn .* dF / c;
else
    % ---- 平面波 / 発散波 / 単一素子: 初到達（Huygens）モデル ---------
    % これらの送信では初到達モデルが閉形式と厳密に一致する。
    tau = inf(numel(xg), 1);
    for e = 1:numel(ex)
        tau = min(tau, td(e) + hypot(xg - ex(e), zg - ez(e)) / c);
    end
end
end

%% ======================================================================
function [rxx, rxz] = normalizeRxPos(rx_pos, Nel)
if isvector(rx_pos)
    rxx = rx_pos(:).';
    rxz = zeros(1, numel(rxx));
elseif size(rx_pos,1) == 2
    rxx = rx_pos(1,:);
    rxz = rx_pos(2,:);
elseif size(rx_pos,2) == 2
    rxx = rx_pos(:,1).';
    rxz = rx_pos(:,2).';
else
    error('das_reference:rxpos', 'rx_pos は [1xNel] か [2xNel] で与えてください。');
end
assert(numel(rxx) == Nel, 'das_reference:rxpos', ...
    'rx_pos の要素数 (%d) が RF の列数 (%d) と一致しません。', numel(rxx), Nel);
end

%% ======================================================================
function [XI, ZI, imgSize] = normalizeGrid(grid_x, grid_z)
if isvector(grid_x) && isvector(grid_z)
    % 軸ベクトル指定 -> meshgrid（画像サイズ = [numel(grid_z) x numel(grid_x)]）
    [XI, ZI] = meshgrid(grid_x(:).', grid_z(:).');
elseif isequal(size(grid_x), size(grid_z))
    XI = grid_x; ZI = grid_z;
else
    error('das_reference:grid', 'grid_x と grid_z のサイズが不整合です。');
end
imgSize = size(XI);
end

%% ======================================================================
function v = getdef(s, name, def)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end
