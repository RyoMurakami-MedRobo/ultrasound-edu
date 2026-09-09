function varargout = wave_animator(action, varargin)
%WAVE_ANIMATOR  Drawing logic for wave propagation, delay curves and the
%   before/after alignment comparison.
%
%   A display-only module called by MAIN_GUI. All state is kept in the
%   UserData of the target axes, so the caller only has to hold axes handles.
%
%   ---- Wavefront animation --------------------------------------------
%   tmax = WAVE_ANIMATOR('setup_propagation', ax, S, txIdx, xlim_m, zlim_m)
%       Prepares the transmit-wavefront and echo fields and returns the final
%       animation time [s]. S is the output of SIM_ENGINE.
%   WAVE_ANIMATOR('draw_propagation', ax, t)
%       Renders the snapshot at time t [s] (updates CData only).
%
%   ---- Delay-curve overlay --------------------------------------------
%   WAVE_ANIMATOR('delay_curve', ax, S, txIdx, tau, pt)
%       Displays the RF data as an image and overlays the delay curve, i.e.
%       which sample of every element is summed for the reconstruction point
%       pt. tau is [1 x Nel] of total two-way times [s] (one row of the
%       delays output of a DAS implementation).
%
%   ---- Before / after alignment ---------------------------------------
%   WAVE_ANIMATOR('alignment', axPre, axPost, S, txIdx, tau, pt)
%       Draws the misaligned waveform bundle (left) next to the delay-
%       corrected, in-phase bundle plus its coherent sum (right).
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
        error('wave_animator:action', 'Unknown action "%s".', action);
end
end

%% ======================================================================
function tmax = setupPropagation(ax, S, txIdx, xrange, zrange)
%SETUPPROPAGATION  Precompute the arrival-time maps and build the axes.
tx = S.tx(txIdx);
c  = S.c;

NX = 241; NZ = 241;
xv = linspace(xrange(1), xrange(2), NX);
zv = linspace(max(zrange(1), 0), zrange(2), NZ);
[X, Z] = meshgrid(xv, zv);

% Time at which the transmit wavefront reaches each point
% (same model as DAS_REFERENCE)
TAU = txArrivalTime(X, Z, tx, c);
TAU = reshape(TAU, size(X));

% Echo of every scatterer: arrival time + distance from the scatterer / c.
%   Speckle phantoms hold far too many scatterers to draw, so only the
%   MAXSCAT strongest ones are animated.
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
Rmat = zeros(NZ, NX, Ns);          % distances from each scatterer, precomputed
for k = 1:Ns
    tauS(k) = txArrivalTime(xs(k), zs(k), tx, c);
    Rmat(:,:,k) = hypot(X - xs(k), Z - zs(k));
end

% Pulse length, used as the visual thickness of the wavefront
sigma = sqrt(2*log(2)) / (pi * S.probe.fc * max(S.probe.bandwidth,1)/100);
sigma = max(sigma, 0.3/S.probe.fc);

cla(ax, 'reset');
hold(ax, 'on');
st = struct();
st.X = X; st.Z = Z; st.TAU = TAU; st.c = c; st.sigma = sigma;
st.xs = xs; st.zs = zs; st.rc = rc; st.tauS = tauS;
st.Rmat = Rmat; st.nDropped = nDropped;
st.xv = xv; st.zv = zv;

% Field image (red = transmit wavefront, blue = scattered echo)
st.hImg = imagesc(ax, xv*1e3, zv*1e3, zeros(NZ, NX));
colormap(ax, blueWhiteRed());
set(ax, 'CLim', [-1 1]);

% Element, scatterer and focus markers
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

% End of the animation: when the echo has reached the farthest element
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
%DRAWPROPAGATION  Update the display to the snapshot at time t.
st = ax.UserData;
if isempty(st) || ~isfield(st, 'hImg') || ~isvalid(st.hImg), return; end

sig = st.sigma;

% Transmit wavefront (positive = red)
F = exp(-0.5*((st.TAU - t)/sig).^2);

% Scattered echoes (negative = blue): concentric circles that start
% spreading once the scatterer has been insonified.
%   Overlapping arcs are combined with a maximum rather than a sum, so they
%   do not saturate and hide the transmit wavefront.
Fe   = zeros(size(F));
hitx = [];
for k = 1:numel(st.xs)
    dt = t - st.tauS(k);
    if dt <= 0, continue; end
    amp = min(abs(st.rc(k)), 1) * 0.9;
    Fe  = max(Fe, amp * exp(-0.5*((st.Rmat(:,:,k)/st.c - dt)/sig).^2));

    % Highlight the receive elements the echo has reached
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
    extra = sprintf('  [showing the %d strongest scatterers, %d hidden]', ...
                    numel(st.xs), st.nDropped);
else
    extra = '';
end
st.hTitle.String = sprintf(['t = %.2f \\mus   (red: transmit wavefront / ' ...
                            'blue: scattered echo)%s'], t*1e6, extra);
end

%% ======================================================================
function drawDelayCurve(ax, S, txIdx, tau, pt)
%DRAWDELAYCURVE  Overlay the summed-sample locations on the RF image.
RF = S.RF(:,:,txIdx);
Nel = size(RF, 2);
tus = (S.t0 + (0:size(RF,1)-1)/S.fs) * 1e6;

cla(ax, 'reset');
% Depth-dependent gain for display only, so shallow echoes do not saturate
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
xlabel(ax, 'Receive element index');
ylabel(ax, 'Time [\mus]');
if any(ok)
    span = max(max(tauUs(ok)) - min(tauUs(ok)), 0.5);
    ylim(ax, [min(tauUs(ok)) - 0.8*span - 1, max(tauUs(ok)) + 0.8*span + 1]);
end
xlim(ax, [0.5 Nel+0.5]);
title(ax, sprintf('Delay curve  (%.2f, %.2f) mm  %d/%d el.', ...
      pt(1)*1e3, pt(2)*1e3, sum(ok), Nel));
end

%% ======================================================================
function drawAlignment(axPre, axPost, S, txIdx, tau, pt)
%DRAWALIGNMENT  Misaligned bundle (left) versus aligned bundle plus sum (right).
RF  = S.RF(:,:,txIdx);
Nel = size(RF, 2);
fc  = S.probe.fc;

tau = tau(:).';
ok  = isfinite(tau);
if ~any(ok)
    cla(axPre,'reset');  text(axPre, .5,.5,'No active element','Units','normalized','HorizontalAlignment','center');
    cla(axPost,'reset'); text(axPost,.5,.5,'No active element','Units','normalized','HorizontalAlignment','center');
    return;
end

trel = linspace(-5/fc, 5/fc, 401);          % +/- 5 periods of the centre frequency
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
title(axPre, sprintf('Before alignment (window at t = %.2f \\mus)', tc*1e6));

drawBundle(axPost, trel*1e6, post, ok, scale, gain, [], []);
% Draw the coherent sum underneath the bundle
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
title(axPost, sprintf('After alignment + sum  (%.2f, %.2f) mm', pt(1)*1e3, pt(2)*1e3));
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
    % Before alignment: show how the sample to be summed scatters per element
    d = (tau - tc) * 1e6;
    plot(ax, d(ok), find(ok), '.', 'Color', [1 .3 0], 'MarkerSize', 9);
    plot(ax, d(ok), find(ok), '-', 'Color', [1 .3 0], 'LineWidth', 1.2);
else
    plot(ax, [0 0], [0 Nel+gain+1], '-', 'Color', [1 .3 0], 'LineWidth', 1.2);
end
hold(ax, 'off');
xlabel(ax, 'Relative time [\mus]');
ylabel(ax, 'Receive element index');
xlim(ax, [tus(1) tus(end)]);
ylim(ax, [0, Nel + gain + 1]);
grid(ax, 'on');
end

%% ----------------------------------------------------------------------
function v = sampleRF(col, t0, fs, times)
%SAMPLERF  Linear interpolation of one RF channel at arbitrary times
%  (zero outside the record).
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
%TXARRIVALTIME  Arrival time [s] of the transmit wavefront.
%  Identical model to the local function of the same name in DAS_REFERENCE.
%  The duplication is deliberate: it keeps the display module independent of
%  the beamformer while still drawing exactly the wavefront the beamformer
%  assumes.
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
%BLUEWHITERED  Blue-white-red diverging colormap (no toolbox required).
n = 128;
u = linspace(0, 1, n).';
lower_ = [u, u, ones(n,1)];                 % blue -> white
upper_ = [ones(n,1), flipud(u), flipud(u)]; % white -> red
cmap = [lower_; upper_];
end
