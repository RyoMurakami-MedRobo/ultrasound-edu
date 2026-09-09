function [bmode_img, delays] = das_custom_template(rf_data, tx_info, rx_pos, grid_x, grid_z, sound_speed, fs)
%DAS_CUSTOM_TEMPLATE  ユーザーが自作 DAS を書くためのひな形。
%
%   [bmode_img, delays] = DAS_CUSTOM_TEMPLATE(rf_data, tx_info, rx_pos, ...
%                                             grid_x, grid_z, sound_speed, fs)
%
%   ■ 使い方
%     1) このファイルを好きな名前（例: my_das.m）でコピーする。
%     2) 関数名をファイル名に合わせて変更する。
%     3) 下の【STEP 1】〜【STEP 5】を書き換える。
%     4) MAIN_GUI の「Custom function」欄にその関数名を入力し、
%        アルゴリズム選択で "Custom" を選ぶと参照実装と比較できる。
%
%   ■ 守るべき契約（これを外すと比較・評価が破綻する）
%     - 引数の順序と意味は DAS_REFERENCE と完全に同一にすること。
%     - rf_data は [Nt x Nel]（時間が第1次元、列が受信素子）。
%     - bmode_img は [numel(grid_z) x numel(grid_x)] の **線形エンベロープ**
%       （対数圧縮しない）。サイズが違うと差分画像・MSE がエラーになる。
%     - delays は [Npix x Nel] の **往復合計時間 [s]**。加算対象外は NaN。
%       nargout < 2 のときは計算をスキップしてよい（実行時間計測を汚さない）。
%     - グリッドは列優先（画素 k = reshape 後の linear index）で並べる。
%
%   ■ 初期状態
%     すぐ動くように「最近傍サンプリング + 矩形アポダイゼーション」の
%     素朴な DAS を実装してある。参照実装（線形補間）との差分画像を見ると
%     補間誤差がサイドローブ状の残差として現れることが確認できる。
%
%   See also DAS_REFERENCE, SIM_ENGINE, MAIN_GUI.

narginchk(7, 7);

%% =====================================================================
%  【STEP 0】入力の整形（通常は触らなくてよい）
%% =====================================================================
if ndims(rf_data) > 2 %#ok<ISMAT>
    rf_data = rf_data(:,:,1);
end
[Nt, Nel] = size(rf_data);

rxx = rx_pos(:).';                       % 受信素子 x 座標 [m]
rxz = zeros(1, Nel);                     % 線形アレイなので z = 0
assert(numel(rxx) == Nel, 'rx_pos の要素数が RF の列数と一致しません。');

if isvector(grid_x) && isvector(grid_z)
    [XI, ZI] = meshgrid(grid_x(:).', grid_z(:).');
else
    XI = grid_x; ZI = grid_z;
end
imgSize = size(XI);
xg = XI(:); zg = ZI(:); Npix = numel(xg);

c  = sound_speed;
t0 = tx_info.t0;                         % RF 第1サンプルの時刻 [s]

%% =====================================================================
%  【STEP 1】送信到達時間 tau_tx(p) を決める  <<< ここを書き換える
%  ---------------------------------------------------------------------
%  参照実装と同じ 2 通りのモデルを実装してある。
%   (a) 初到達（Huygens）モデル ... 平面波・発散波・単一素子で厳密
%         tau_tx(p) = min_e ( delays(e) + |p - elem_e| / c )
%   (b) 仮想音源モデル ......... 集束波用
%         tau_tx(p) = T_F + s * |p - F| / c   (s = -1 焦点手前 / +1 以遠)
%  (a) を集束波にも使うと焦点以遠で開口端の弱い端部波を拾い、深さ方向に
%  1 mm 程度ずれる。この切り分けを外すとどうなるかを試すのも良い教材。
%% =====================================================================
ex  = tx_info.elem_x;
ez  = tx_info.elem_z;
td  = tx_info.delays;
ap  = tx_info.apod;
act = isfinite(td) & ap(:).' ~= 0;       % 励振している素子だけを使う
exA = ex(act); ezA = ez(act); tdA = td(act);

src = tx_info.src_xz;                    % 集束点 / 仮想音源 [x z] （無いと NaN）
useVirtualSource = strcmpi(tx_info.scheme,'focused') && all(isfinite(src)) && ...
                   src(2) > 0 && ~tx_info.single_element;

if useVirtualSource
    Rf   = hypot(exA - src(1), ezA - src(2));
    Tf   = mean(tdA + Rf / c);           % 波面が焦点に収束する時刻 [s]
    D    = hypot(src(1), src(2));
    u    = [src(1) src(2)] / D;
    proj = xg*u(1) + zg*u(2);
    dF   = hypot(xg - src(1), zg - src(2));
    sgn  = ones(Npix,1);  sgn(proj <= D) = -1;
    tau_tx = Tf + sgn .* dF / c;
else
    tau_tx = inf(Npix, 1);
    for e = 1:numel(exA)
        tau_tx = min(tau_tx, tdA(e) + hypot(xg - exA(e), zg - ezA(e)) / c);
    end
end

%% =====================================================================
%  【STEP 2】受信開口の重み（f 値・アポダイゼーション）  <<< 書き換え可
%% =====================================================================
fnum = tx_info.fnumber;
if ~isfinite(fnum) || fnum <= 0
    halfAp = inf(Npix, 1);               % 全開口
else
    halfAp = zg ./ (2*fnum);             % 深さに比例して開口を広げる
end

%% =====================================================================
%  【STEP 3】遅延加算のメインループ  <<< ここを書き換える
%  ---------------------------------------------------------------------
%  初期実装は「最近傍サンプリング」。round() を floor()+線形補間に
%  変えると参照実装と一致するはずである（自作コードの検算に使える）。
%% =====================================================================
bf = zeros(Npix, 1);
if nargout > 1
    delays = nan(Npix, Nel);
end

for e = 1:Nel
    dx  = xg - rxx(e);
    tau = tau_tx + hypot(dx, zg - rxz(e)) / c;    % 往復合計時間 [s]

    w   = double(abs(dx) <= halfAp);              % 矩形アポダイゼーション

    idx = round((tau - t0) * fs + 1);             % ★ 最近傍サンプリング
    ok  = w > 0 & idx >= 1 & idx <= Nt;
    idxc = min(max(idx, 1), Nt);

    val = rf_data(idxc, e);
    val(~ok) = 0;

    bf = bf + w .* val;

    if nargout > 1
        d = tau; d(~ok) = NaN;
        delays(:, e) = d;
    end
end

%% =====================================================================
%  【STEP 4】エンベロープ検波  <<< 書き換え可（例: ヒルベルト以外の手法）
%% =====================================================================
bf2 = reshape(bf, imgSize);
if imgSize(1) >= 4
    bmode_img = abs(hilbert(bf2));
else
    bmode_img = abs(bf2);
end

%% =====================================================================
%  【STEP 5】後処理（TGC・スペックル低減など）を入れるならここ
%% =====================================================================
% bmode_img = ...;

end
