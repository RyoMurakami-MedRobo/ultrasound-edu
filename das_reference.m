function [bmode_img, delays, trace] = das_reference(rf_data, tx_info, rx_pos, grid_x, grid_z, sound_speed, fs)
% Optional third output trace records actual channel increments and coherent RF.
% See docs/HW3_DAS.md; retain the accumulation hook in custom implementations.
%DAS_REFERENCE  Textbook Delay-And-Sum beamformer (reference implementation).
%
%   [bmode_img, delays] = DAS_REFERENCE(rf_data, tx_info, rx_pos, ...
%                                       grid_x, grid_z, sound_speed, fs)
%
%   This is the unified interface required by the custom-DAS benchmarking
%   framework of this tool. Any user-written algorithm (see
%   DAS_CUSTOM_TEMPLATE) must follow the same argument list and return values.
%
%   ---- Inputs ---------------------------------------------------------
%   rf_data     [Nt x Nel] RF data. **Time along the first dimension**,
%               columns are receive elements. If a 3-D array is passed, only
%               the first transmit event is used.
%   tx_info     Transmit settings struct. The following fields are read:
%                 .delays  [1 x Nel] transmit delays [s] (NaN when inactive)
%                 .apod    [1 x Nel] transmit apodisation (0 when inactive)
%                 .elem_x  [1 x Nel] transmit element x coordinates [m]
%                 .elem_z  [1 x Nel] transmit element z coordinates [m]
%                 .t0      time of the first RF sample [s]
%                 .fnumber receive f-number (0 or non-finite = full aperture)
%                 .rx_apod 'rect' or 'hann'
%   rx_pos      [1 x Nel] or [2 x Nel] receive element coordinates [m]. A
%               single row is treated as x coordinates with z = 0. This
%               argument is authoritative for the receive geometry, whereas
%               tx_info.elem_* describes the transmit geometry.
%   grid_x,grid_z Reconstruction grid [m]. Vectors are expanded with meshgrid
%               and the image size becomes [numel(grid_z) x numel(grid_x)].
%               Matrices of equal size are used directly as a point cloud.
%   sound_speed Speed of sound [m/s]
%   fs          Sampling frequency [Hz]
%
%   ---- Outputs --------------------------------------------------------
%   bmode_img   [Nz x Nx] **linear envelope** (no log compression). The caller
%               displays 20*log10(img/max(img)). Metrics such as MSE and FWHM
%               are only meaningful on linear values, which is why this
%               interface returns the linear envelope.
%   delays      [Nz*Nx x Nel] **total two-way time [s]** for every pixel and
%               element. Entries excluded by the f-number or by an inactive
%               transmit element are NaN. Not computed when nargout < 2, so
%               benchmark timings are not inflated.
%
%   ---- Algorithm ------------------------------------------------------
%   1) Transmit arrival time tau_tx(p)
%      * plane / diverging / single element: first-arrival (Huygens) model
%          tau_tx(p) = min_e ( txdelay(e) + |p - e| / c )
%        For these schemes this is exactly equal to the closed forms
%        (x*sin + z*cos, or the virtual-source expression), so no branching
%        is required.
%      * focused: virtual-source model
%          tau_tx(p) = T_F + s * |p - F| / c  (s = -1 before, +1 beyond focus)
%        The first-arrival model would latch onto the weak edge wave from the
%        aperture rim beyond the focus.
%   2) Receive arrival time  tau_rx(p,e) = |p - e_rx| / c
%   3) Convert tau = tau_tx + tau_rx to a sample position and interpolate
%      linearly.
%   4) Apply the f-number aperture limit and the receive apodisation, then
%      sum across elements.
%   5) Envelope detection along depth using the Hilbert transform.
%
%   See also SIM_ENGINE, DAS_CUSTOM_TEMPLATE, WAVE_ANIMATOR.

narginchk(7, 7);

%% ---- Normalise the inputs -------------------------------------------
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
if ~isfinite(fnum) || fnum <= 0, fnum = 0; end   % 0 = full aperture

%% ---- 1) Transmit arrival time ---------------------------------------
tau_tx = txArrivalTime(XI, ZI, tx_info, c);

%% ---- Preparation ----------------------------------------------------
bf = zeros(Npix, 1);
% Optional execution trace: record the exact increment at the summation site.
% Keep this hook when editing DAS. Request only a scan line to bound memory.
trace = [];
if nargout > 2
    trace.contributions = zeros(Npix, Nel);
end
if nargout > 1
    delays = nan(Npix, Nel);
end

% Half-width of the receive aperture implied by the f-number (grows with depth)
if fnum > 0
    halfAp = ZI(:) ./ (2*fnum);
else
    halfAp = inf(Npix, 1);
end

xg = XI(:); zg = ZI(:);

%% ---- 2)-4) Delay-and-sum over the elements --------------------------
for e = 1:Nel
    dx  = xg - rxx(e);
    rrx = hypot(dx, zg - rxz(e));
    tau = tau_tx + rrx / c;                       % total two-way time [s]

    % Receive aperture limit
    inAp = abs(dx) <= halfAp;
    switch rxApod
        case 'hann'
            u = zeros(Npix,1);
            u(inAp) = dx(inAp) ./ max(halfAp(inAp), eps);   % -1..1
            w = inAp .* (0.5 * (1 + cos(pi*u)));
        otherwise
            w = double(inAp);
    end

    % Sample the RF data by linear interpolation
    pos  = (tau - t0) * fs + 1;
    i0   = floor(pos);
    frac = pos - i0;
    ok   = w > 0 & i0 >= 1 & i0 < Nt;
    i0c  = min(max(i0, 1), Nt-1);
    val  = rf_data(i0c, e) .* (1-frac) + rf_data(i0c+1, e) .* frac;
    val(~ok) = 0;

    if nargout > 2, previous = bf; end
    bf = bf + w .* val;
    if nargout > 2
        trace.contributions(:, e) = bf - previous;
    end

    if nargout > 1
        d = tau;
        d(~ok) = NaN;
        delays(:, e) = d;
    end
end

%% ---- 5) Envelope detection ------------------------------------------
bmode_img = envelopeZ(reshape(bf, imgSize));
if nargout > 2, trace.coherent = reshape(bf, imgSize); end
end % das_reference

%% ======================================================================
function env = envelopeZ(x)
%ENVELOPEZ  Envelope detection along the first (depth) dimension.
if size(x,1) < 4
    env = abs(x);                       % degenerate single-point / tiny grid
elseif ~isempty(which('hilbert'))
    env = abs(hilbert(x));              % Signal Processing Toolbox
else
    % Build the analytic signal with an FFT so no toolbox is required
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
%TXARRIVALTIME  Time [s] at which the transmit wavefront reaches each pixel,
%  returned as [Npix x 1].
ex = getdef(tx_info, 'elem_x', []);
ez = getdef(tx_info, 'elem_z', []);
td = getdef(tx_info, 'delays', []);
ap = getdef(tx_info, 'apod', []);

if isempty(ex) || isempty(td)
    error('das_reference:txinfo', 'tx_info must provide elem_x and delays.');
end
if isempty(ez), ez = zeros(size(ex)); end
if isempty(ap), ap = ones(size(ex)); end

act = isfinite(td) & ap(:).' ~= 0;
if ~any(act)
    error('das_reference:noActiveTx', 'No active transmit element.');
end
ex = ex(act); ez = ez(act); td = td(act);
xg = XI(:); zg = ZI(:);

src    = getdef(tx_info, 'src_xz', [NaN NaN]);
scheme = getdef(tx_info, 'scheme', '');
single = getdef(tx_info, 'single_element', false);
isFocused = strcmpi(scheme,'focused') && numel(src)==2 && all(isfinite(src)) && ...
            src(2) > 0 && ~single;

if isFocused
    % ---- Focused transmit: virtual-source model ---------------------
    % The first-arrival (min) model is exact before the focus but beyond it
    % picks up the weak edge wave radiated by the aperture rim, which arrives
    % earlier than the main energy. Focused transmits are therefore treated
    % as a spherical wave passing through the focus:
    %   tau_tx(p) = T_F + s * |p - F| / c   (s = -1 before, +1 beyond focus)
    % T_F, the instant the wavefront converges on the focus, is constant
    % across elements by definition of the delay law, so it can be recovered
    % from the delays themselves regardless of how they were normalised.
    Rf   = hypot(ex - src(1), ez - src(2));
    Tf   = mean(td + Rf / c);
    D    = hypot(src(1), src(2));
    u    = [src(1) src(2)] / D;                 % unit vector along the beam axis
    proj = xg*u(1) + zg*u(2);                   % projection onto the beam axis
    dF   = hypot(xg - src(1), zg - src(2));
    sgn  = ones(numel(xg), 1);
    sgn(proj <= D) = -1;
    tau  = Tf + sgn .* dF / c;
else
    % ---- Plane / diverging / single element: first-arrival model ----
    % For these schemes the first-arrival model matches the closed form
    % exactly.
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
    error('das_reference:rxpos', 'rx_pos must be [1xNel] or [2xNel].');
end
assert(numel(rxx) == Nel, 'das_reference:rxpos', ...
    'rx_pos has %d elements but the RF data has %d columns.', numel(rxx), Nel);
end

%% ======================================================================
function [XI, ZI, imgSize] = normalizeGrid(grid_x, grid_z)
if isvector(grid_x) && isvector(grid_z)
    % Axis vectors -> meshgrid (image size = [numel(grid_z) x numel(grid_x)])
    [XI, ZI] = meshgrid(grid_x(:).', grid_z(:).');
elseif isequal(size(grid_x), size(grid_z))
    XI = grid_x; ZI = grid_z;
else
    error('das_reference:grid', 'grid_x and grid_z have incompatible sizes.');
end
imgSize = size(XI);
end

%% ======================================================================
function v = getdef(s, name, def)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end
