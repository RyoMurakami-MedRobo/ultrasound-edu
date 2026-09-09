function out = sim_engine(cfg)
%SIM_ENGINE  トランスデューサ／送信シーケンス／散乱体設定から RF データを生成する。
%
%   out = SIM_ENGINE(cfg)
%
%   MUST (Matlab UltraSound Toolbox) の SIMUS をバックエンドとして使用する。
%   MUST が検出できない場合は、解析的な円筒波モデル（モック）に自動的に
%   フォールバックするため、MUST 未インストール環境でも GUI の全機能を
%   動作確認できる。
%
%   ---- 座標系 ---------------------------------------------------------
%   x : アレイに平行（左端素子 → 右端素子の向きが +x）、アレイ中心が x = 0
%   z : アレイに垂直で下向き（深さ方向）。素子は z = 0 上に並ぶ。
%   ステアリング角 theta は z 軸から測り、+x 側へ振る向きを正とする。
%
%   ---- 入力 cfg -------------------------------------------------------
%   cfg.probe.Nelements   素子数                          (既定 64)
%   cfg.probe.pitch       ピッチ [m]                      (既定 0.30e-3)
%   cfg.probe.fc          中心周波数 [Hz]                 (既定 5e6)
%   cfg.probe.bandwidth   -6dB 比帯域 [%]                 (既定 75)
%   cfg.probe.kerf        素子間隙 [m]                    (既定 pitch*0.1)
%   cfg.probe.height      素子高さ [m]                    (既定 5e-3)
%   cfg.probe.elevfocus   エレベーション焦点 [m]          (既定 20e-3)
%   cfg.medium.c          音速 [m/s]                      (既定 1540)
%   cfg.acq.fs_factor     fs = fs_factor * fc             (既定 4、>=4 を推奨)
%   cfg.tx.scheme         'plane' | 'focused' | 'diverging'
%   cfg.tx.angle_deg      ステアリング角 [deg]            (既定 0)
%   cfg.tx.focus_mm       集束深度 [mm]。複数指定でマルチフォーカス送信
%                         （送信イベントを深度ごとに分けて生成する）
%   cfg.tx.src_mm         発散波の仮想音源深度 [mm]（正値、実際は z<0 に配置）
%   cfg.tx.single_element 単一素子送信なら true           (既定 false)
%   cfg.tx.element_index  単一素子送信で励振する素子番号  (既定 中央)
%   cfg.scat.x/.z/.rc     散乱体の x,z 座標 [m] と反射強度（同サイズのベクトル）
%   cfg.recon.zmax        記録したい最大深度 [m]          (既定 40e-3)
%   cfg.options.force_mock  true で MUST があってもモックを使う (既定 false)
%
%   ---- 出力 out -------------------------------------------------------
%   out.RF      [Nt x Nelements x Ntx] RF データ。**時間は列方向（第1次元）**、
%               列インデックスが受信素子番号。第3次元が送信イベント。
%   out.t0      RF 第1サンプルの時刻 [s]（本エンジンでは常に 0）
%   out.fs      サンプリング周波数 [Hz]
%   out.c       音速 [m/s]
%   out.time    [Nt x 1] 時間軸 = t0 + (0:Nt-1)/fs
%   out.rx_pos  [1 x Nelements] 受信素子の x 座標 [m]（z = 0）
%   out.tx      [1 x Ntx] struct 配列。ビームフォーマへ渡す tx_info。
%   out.backend 'MUST' もしくは 'mock'
%   out.elapsed シミュレーション所要時間 [s]
%
%   tx_info の各フィールド（DAS 実装が参照する契約）:
%     .delays   [1 x Nel] 送信遅延 [s]（非励振素子は NaN、min(有効)=0 に正規化）
%     .apod     [1 x Nel] 送信アポダイゼーション（非励振素子は 0）
%     .elem_x   [1 x Nel] 送信素子 x 座標 [m]
%     .elem_z   [1 x Nel] 送信素子 z 座標 [m]（線形アレイなので全て 0）
%     .t0       RF 第1サンプルの時刻 [s]
%     .c,.fc,.fs 音速・中心周波数・サンプリング周波数
%     .fnumber  受信 f 値（0 で全開口）。※統一 I/F に引数が無いためここで渡す
%     .rx_apod  受信アポダイゼーション 'rect' | 'hann'
%     .scheme,.angle_deg,.focus_mm  送信条件（表示・注釈用）
%     .src_xz   [x z] 集束点／仮想音源の座標 [m]（教育用オーバーレイに使用）
%
%   See also DAS_REFERENCE, DAS_CUSTOM_TEMPLATE, WAVE_ANIMATOR, MAIN_GUI.

%% ---------------- 既定値の展開 ----------------
if nargin < 1 || isempty(cfg), cfg = struct(); end

probe   = getsub(cfg, 'probe');
medium  = getsub(cfg, 'medium');
acq     = getsub(cfg, 'acq');
tx      = getsub(cfg, 'tx');
scat    = getsub(cfg, 'scat');
recon   = getsub(cfg, 'recon');
opts    = getsub(cfg, 'options');

Nel        = getdef(probe, 'Nelements', 64);
pitch      = getdef(probe, 'pitch',     0.30e-3);
fc         = getdef(probe, 'fc',        5e6);
bandwidth  = getdef(probe, 'bandwidth', 75);
kerf       = getdef(probe, 'kerf',      pitch*0.1);
height     = getdef(probe, 'height',    5e-3);
elevfocus  = getdef(probe, 'elevfocus', 20e-3);
c          = getdef(medium,'c',         1540);
fsFactor   = getdef(acq,   'fs_factor', 4);
zmax       = getdef(recon, 'zmax',      40e-3);
forceMock  = getdef(opts,  'force_mock',false);

validateattributes(Nel,   {'numeric'},{'scalar','integer','>=',2,'<=',1024},mfilename,'Nelements');
validateattributes(pitch, {'numeric'},{'scalar','positive','finite'},mfilename,'pitch');
validateattributes(fc,    {'numeric'},{'scalar','positive','finite'},mfilename,'fc');
if fsFactor < 4
    warning('sim_engine:lowFs', ...
        ['fs_factor = %.2f は小さすぎます（fs = %.2f*fc）。補間ベースの DAS で ' ...
         'エイリアシングが生じます。4 以上を推奨します。'], fsFactor, fsFactor);
end
fs = fsFactor * fc;

% 素子座標（アレイ中心が x = 0）
elem_x = ((0:Nel-1) - (Nel-1)/2) * pitch;
elem_z = zeros(1, Nel);

%% ---------------- 送信シーケンスの構築 ----------------
scheme   = lower(getdef(tx, 'scheme', 'plane'));
angleDeg = getdef(tx, 'angle_deg', 0);
theta    = angleDeg * pi/180;
singleEl = getdef(tx, 'single_element', false);
elIdx    = getdef(tx, 'element_index', round((Nel+1)/2));

switch scheme
    case 'plane'
        srcDepths = NaN;                                   % 焦点なし
    case 'focused'
        fmm = getdef(tx, 'focus_mm', 20);
        srcDepths = fmm(:).' * 1e-3;                       % 正 = 集束
        srcDepths(srcDepths <= 0) = 1e-3;
    case 'diverging'
        smm = getdef(tx, 'src_mm', 10);
        srcDepths = -abs(smm(1)) * 1e-3;                   % 負 = 仮想音源
    otherwise
        error('sim_engine:scheme', '未知の送信モード "%s" です。', scheme);
end

Ntx = numel(srcDepths);
txArr = repmat(makeEmptyTx(), 1, Ntx);
for k = 1:Ntx
    [d, a, srcxz] = txLaw(elem_x, c, theta, srcDepths(k), scheme, singleEl, elIdx);
    txArr(k).delays   = d;
    txArr(k).apod     = a;
    txArr(k).elem_x   = elem_x;
    txArr(k).elem_z   = elem_z;
    txArr(k).t0       = 0;
    txArr(k).c        = c;
    txArr(k).fc       = fc;
    txArr(k).fs       = fs;
    txArr(k).fnumber  = getdef(recon, 'fnumber', 1.5);
    txArr(k).rx_apod  = getdef(recon, 'rx_apod', 'rect');
    txArr(k).scheme   = scheme;
    txArr(k).angle_deg= angleDeg;
    txArr(k).focus_mm = srcDepths(k)*1e3;
    txArr(k).src_xz   = srcxz;
    txArr(k).single_element = singleEl;
end

%% ---------------- 散乱体 ----------------
xs = getdef(scat, 'x', 0);   xs = xs(:).';
zs = getdef(scat, 'z', 20e-3); zs = zs(:).';
rc = getdef(scat, 'rc', 1);  rc = rc(:).';
if isscalar(rc) && numel(xs) > 1, rc = repmat(rc, 1, numel(xs)); end
assert(isequal(numel(xs), numel(zs), numel(rc)), ...
    'sim_engine:scat', '散乱体の x, z, rc は同じ要素数である必要があります。');
keep = isfinite(xs) & isfinite(zs) & isfinite(rc) & zs > 0;
xs = xs(keep); zs = zs(keep); rc = rc(keep);
if isempty(xs)
    xs = 0; zs = max(zmax/2, 1e-3); rc = 0;      % 空ファントムでも RF 長を確保
end

%% ---------------- 記録長の決定 ----------------
apHalf  = max(abs(elem_x)) + pitch;
zrec    = max([zmax, max(zs)]) * 1.05;
sigma   = pulseSigma(fc, bandwidth);
maxDel  = max([0, max(cellfun(@(d) max(d(isfinite(d))), {txArr.delays}))]);
Tend    = maxDel + 2*hypot(zrec, 2*apHalf)/c + 8*sigma;
Nt      = max(64, ceil(Tend*fs));
tvec    = (0:Nt-1).'/fs;

%% ---------------- バックエンド選択 ----------------
hasMUST = ~isempty(which('simus')) && ~isempty(which('pfield'));
useMUST = hasMUST && ~forceMock;

tic;
if useMUST
    try
        RF = runMUST(xs, zs, rc, txArr, Nel, pitch, kerf, fc, bandwidth, ...
                     height, elevfocus, c, fs, Nt);
        backend = 'MUST';
    catch ME
        warning('sim_engine:mustFailed', ...
            'MUST(simus) の実行に失敗したためモックに切り替えます: %s', ME.message);
        RF = runMock(xs, zs, rc, txArr, elem_x, elem_z, c, fc, sigma, tvec);
        backend = 'mock (MUST failed)';
    end
else
    RF = runMock(xs, zs, rc, txArr, elem_x, elem_z, c, fc, sigma, tvec);
    if hasMUST
        backend = 'mock (forced)';
    else
        backend = 'mock (MUST not found)';
    end
end
elapsed = toc;

% 振幅を正規化（バックエンド間で表示スケールを揃える）
pk = max(abs(RF(:)));
if pk > 0, RF = RF / pk; end

%% ---------------- 出力 ----------------
out = struct();
out.RF      = RF;
out.t0      = 0;
out.fs      = fs;
out.c       = c;
out.time    = (0:size(RF,1)-1).'/fs;
out.rx_pos  = elem_x;
out.rx_pos_z= elem_z;
out.tx      = txArr;
out.backend = backend;
out.elapsed = elapsed;
out.probe   = struct('Nelements',Nel,'pitch',pitch,'fc',fc,'bandwidth',bandwidth, ...
                     'kerf',kerf,'height',height,'elevfocus',elevfocus);
out.scat    = struct('x',xs,'z',zs,'rc',rc);
end % sim_engine

%% ======================================================================
%  MUST アダプタ（MUST 依存コードはこの関数の中だけに閉じ込める）
%  ----------------------------------------------------------------------
%  検証済みの呼び出し規約（biomecardio.com / MUST 公式ドキュメント）:
%    RF = SIMUS(X,Z,RC,DELAYS,PARAM)   % 2-D 構文。RF の列数 = 素子数
%    PARAM.fc, PARAM.pitch, PARAM.width または kerf は必須
%    PARAM.fs 既定 = 4*fc、PARAM.bandwidth はパルスエコー -6dB 比帯域 [%]
%    PARAM.c 既定 1540、PARAM.TXapodization、PARAM.RXdelay 既定 0
%    時間原点は t = 0（例示コードが t = (0:size(RF,1)-1)/param.fs を使用）
%  MUST のバージョン差で引数が変わった場合はここだけを直せばよい。
%% ======================================================================
function RF = runMUST(xs, zs, rc, txArr, Nel, pitch, kerf, fc, bandwidth, ...
                      height, elevfocus, c, fs, Nt)
param = struct();
param.Nelements    = Nel;
param.pitch        = pitch;
param.kerf         = kerf;
param.width        = pitch - kerf;
param.fc           = fc;
param.bandwidth    = bandwidth;
param.height       = height;
param.focus        = elevfocus;
param.c            = c;
param.fs           = fs;
param.RXdelay      = zeros(1, Nel);
param.t0           = 0;

opt = struct('WaitBar', false, 'ParPool', false);

Ntx = numel(txArr);
RFc = cell(1, Ntx);
for k = 1:Ntx
    d = txArr(k).delays;
    a = txArr(k).apod;
    d(~isfinite(d)) = 0;                 % 非励振素子は遅延 0 + apod 0 で無効化
    param.TXapodization = a;
    RFc{k} = simus(xs, zs, rc, d, param, opt);
end

% 送信イベント間で行数が揃わないことがあるためゼロ詰めして連結する
Nrow = max([Nt, cellfun(@(r) size(r,1), RFc)]);
RF = zeros(Nrow, Nel, Ntx);
for k = 1:Ntx
    RF(1:size(RFc{k},1), :, k) = RFc{k};
end
end

%% ======================================================================
%  モックバックエンド（MUST 非依存の解析的 2-D 円筒波モデル）
%  各散乱体について
%     (1) 全送信素子から届く球面波を重ね合わせて散乱体位置の場を作り
%     (2) その場を受信素子まで再度伝搬させて RF に加算する
%  回折・素子指向性の厳密モデルは持たないが、遅延構造は物理的に正しいので
%  DAS の整相を検証する目的には十分。
%% ======================================================================
function RF = runMock(xs, zs, rc, txArr, elem_x, elem_z, c, fc, sigma, tvec)
Nel  = numel(elem_x);
Nt   = numel(tvec);
Ntx  = numel(txArr);
fs   = 1/(tvec(2)-tvec(1));
RF   = zeros(Nt, Nel, Ntx);

for k = 1:Ntx
    txd = txArr(k).delays;
    apo = txArr(k).apod;
    act = isfinite(txd) & apo ~= 0;
    if ~any(act), continue; end
    acc = zeros(Nt, Nel);

    for is = 1:numel(xs)
        if rc(is) == 0, continue; end

        % --- (1) 送信: 散乱体位置での場 ------------------------------
        rt   = hypot(xs(is) - elem_x(act), zs(is) - elem_z(act));   % 1 x Nact
        taut = txd(act) + rt/c;
        ampt = apo(act) .* (zs(is)./rt) ./ sqrt(rt);                % 斜入射 + 円筒拡散
        sfield = gpulse(tvec - taut, fc, sigma) * ampt(:);          % Nt x 1

        % --- (2) 受信: 各素子への再伝搬 ------------------------------
        rr    = hypot(xs(is) - elem_x, zs(is) - elem_z);            % 1 x Nel
        ampr  = (zs(is)./rr) ./ sqrt(rr);
        pos   = ((tvec - rr/c) - tvec(1)) * fs + 1;                 % Nt x Nel
        i0    = floor(pos);
        frac  = pos - i0;
        ok    = i0 >= 1 & i0 < Nt;
        i0c   = min(max(i0, 1), Nt-1);
        v     = sfield(i0c).*(1-frac) + sfield(i0c+1).*frac;
        v(~ok)= 0;

        acc = acc + rc(is) * (v .* ampr);
    end
    RF(:,:,k) = acc;
end
end

%% ======================================================================
%  補助関数
%% ======================================================================
function [d, a, srcxz] = txLaw(elem_x, c, theta, srcDepth, scheme, singleEl, elIdx)
%TXLAW  送信遅延則。焦点／仮想音源はビーム軸上 D*(sin th, cos th) に置く。
Nel = numel(elem_x);
a = ones(1, Nel);
if strcmp(scheme, 'plane') || ~isfinite(srcDepth)
    d = elem_x * sin(theta) / c;                 % 平面波
    srcxz = [NaN NaN];
else
    Px = srcDepth * sin(theta);
    Pz = srcDepth * cos(theta);
    R  = hypot(elem_x - Px, Pz);
    if srcDepth > 0
        d = (max(R) - R) / c;                    % 集束波（外側が先に発射）
    else
        d = (R - min(R)) / c;                    % 発散波（仮想音源から等距離）
    end
    srcxz = [Px Pz];
end

if singleEl
    elIdx = min(max(round(elIdx), 1), Nel);
    a(:) = 0; a(elIdx) = 1;
    d(:) = NaN; d(elIdx) = 0;
    srcxz = [elem_x(elIdx) 0];
end

d = d - min(d(isfinite(d)));                     % min(有効遅延) = 0 に正規化
d(a == 0) = NaN;                                 % 非励振素子は NaN
end

function s = makeEmptyTx()
s = struct('delays',[],'apod',[],'elem_x',[],'elem_z',[],'t0',0,'c',1540, ...
           'fc',5e6,'fs',20e6,'fnumber',1.5,'rx_apod','rect','scheme','plane', ...
           'angle_deg',0,'focus_mm',NaN,'src_xz',[NaN NaN],'single_element',false);
end

function sigma = pulseSigma(fc, bwPercent)
%PULSESIGMA  比帯域 [%] からガウス包絡の標準偏差 [s] を求める。
bw = max(bwPercent, 1)/100;
sigma = sqrt(2*log(2)) / (pi * fc * bw);
end

function y = gpulse(t, fc, sigma)
%GPULSE  ガウス変調正弦パルス。
y = exp(-0.5*(t/sigma).^2) .* cos(2*pi*fc*t);
end

function s = getsub(c, name)
if isstruct(c) && isfield(c, name) && isstruct(c.(name)), s = c.(name); else, s = struct(); end
end

function v = getdef(s, name, def)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end
