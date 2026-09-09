function out = sim_engine(cfg)
%SIM_ENGINE  Generate RF data from a transducer / transmit / scatterer setup.
%
%   out = SIM_ENGINE(cfg)
%
%   Uses SIMUS from MUST (Matlab UltraSound Toolbox) as the backend. When MUST
%   cannot be found, the function automatically falls back to an analytic
%   cylindrical-wave model (the "mock" backend), so the whole GUI remains
%   usable on machines where MUST is not installed.
%
%   ---- Coordinate system ----------------------------------------------
%   x : parallel to the array, pointing from the first (leftmost) element to
%       the last (rightmost) element. x = 0 at the array centre.
%   z : perpendicular to the array, pointing downward (depth). Elements lie
%       on z = 0.
%   The steering angle theta is measured from the z axis and is positive
%   towards +x.
%
%   ---- Input: cfg -----------------------------------------------------
%   cfg.probe.Nelements   number of elements                (default 64)
%   cfg.probe.pitch       element pitch [m]                 (default 0.30e-3)
%   cfg.probe.fc          centre frequency [Hz]             (default 5e6)
%   cfg.probe.bandwidth   -6 dB fractional bandwidth [%]    (default 75)
%   cfg.probe.kerf        kerf [m]                          (default pitch*0.1)
%   cfg.probe.height      element height [m]                (default 5e-3)
%   cfg.probe.elevfocus   elevation focus [m]               (default 20e-3)
%   cfg.medium.c          speed of sound [m/s]              (default 1540)
%   cfg.acq.fs_factor     fs = fs_factor * fc               (default 4, >= 4)
%   cfg.tx.scheme         'plane' | 'focused' | 'diverging'
%   cfg.tx.angle_deg      steering angle [deg]              (default 0)
%   cfg.tx.focus_mm       focal depth(s) [mm]. Several values request a
%                         multi-focus sequence (one transmit event each).
%   cfg.tx.src_mm         virtual source depth [mm] for diverging waves
%                         (positive value; the source is placed at z < 0)
%   cfg.tx.single_element true for single-element transmit  (default false)
%   cfg.tx.element_index  element fired in single-element mode (default centre)
%   cfg.scat.x/.z/.rc     scatterer x, z [m] and reflection coefficients
%                         (vectors of equal length)
%   cfg.recon.zmax        maximum depth to record [m]       (default 40e-3)
%   cfg.options.force_mock  true to use the mock backend even if MUST exists
%
%   ---- Output: out ----------------------------------------------------
%   out.RF      [Nt x Nelements x Ntx] RF data. **Time runs along the first
%               dimension**, columns are receive elements, the third
%               dimension indexes transmit events.
%   out.t0      time of the first RF sample [s] (always 0 in this engine)
%   out.fs      sampling frequency [Hz]
%   out.c       speed of sound [m/s]
%   out.time    [Nt x 1] time axis = t0 + (0:Nt-1)/fs
%   out.rx_pos  [1 x Nelements] receive element x coordinates [m]
%   out.tx      [1 x Ntx] struct array holding the tx_info of each event
%   out.backend 'MUST' or 'mock'
%   out.elapsed simulation time [s]
%
%   Fields of tx_info (the contract every DAS implementation relies on):
%     .delays   [1 x Nel] transmit delays [s] (NaN on inactive elements,
%               normalised so that min(active) = 0)
%     .apod     [1 x Nel] transmit apodisation (0 on inactive elements)
%     .elem_x   [1 x Nel] transmit element x coordinates [m]
%     .elem_z   [1 x Nel] transmit element z coordinates [m] (0 for a linear
%               array)
%     .t0       time of the first RF sample [s]
%     .c,.fc,.fs  speed of sound, centre frequency, sampling frequency
%     .fnumber  receive f-number (0 means full aperture). The unified DAS
%               interface has no argument for it, so it travels in here.
%     .rx_apod  receive apodisation, 'rect' or 'hann'
%     .scheme,.angle_deg,.focus_mm  transmit settings (for display/annotation)
%     .src_xz   [x z] focal point or virtual source [m] (used by the animator)
%
%   See also DAS_REFERENCE, DAS_CUSTOM_TEMPLATE, WAVE_ANIMATOR, MAIN_GUI.

%% ---------------- Defaults ----------------
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
        ['fs_factor = %.2f is too low (fs = %.2f*fc). The interpolation-based ' ...
         'DAS will alias. A value of 4 or more is recommended.'], fsFactor, fsFactor);
end
fs = fsFactor * fc;

% Element coordinates (array centred on x = 0)
elem_x = ((0:Nel-1) - (Nel-1)/2) * pitch;
elem_z = zeros(1, Nel);

%% ---------------- Transmit sequence ----------------
scheme   = lower(getdef(tx, 'scheme', 'plane'));
angleDeg = getdef(tx, 'angle_deg', 0);
theta    = angleDeg * pi/180;
singleEl = getdef(tx, 'single_element', false);
elIdx    = getdef(tx, 'element_index', round((Nel+1)/2));

switch scheme
    case 'plane'
        srcDepths = NaN;                                   % no focal point
    case 'focused'
        fmm = getdef(tx, 'focus_mm', 20);
        srcDepths = fmm(:).' * 1e-3;                       % positive = focused
        srcDepths(srcDepths <= 0) = 1e-3;
    case 'diverging'
        smm = getdef(tx, 'src_mm', 10);
        srcDepths = -abs(smm(1)) * 1e-3;                   % negative = virtual source
    otherwise
        error('sim_engine:scheme', 'Unknown transmit scheme "%s".', scheme);
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

%% ---------------- Scatterers ----------------
xs = getdef(scat, 'x', 0);   xs = xs(:).';
zs = getdef(scat, 'z', 20e-3); zs = zs(:).';
rc = getdef(scat, 'rc', 1);  rc = rc(:).';
if isscalar(rc) && numel(xs) > 1, rc = repmat(rc, 1, numel(xs)); end
assert(isequal(numel(xs), numel(zs), numel(rc)), ...
    'sim_engine:scat', 'Scatterer x, z and rc must have the same number of elements.');
keep = isfinite(xs) & isfinite(zs) & isfinite(rc) & zs > 0;
xs = xs(keep); zs = zs(keep); rc = rc(keep);
if isempty(xs)
    xs = 0; zs = max(zmax/2, 1e-3); rc = 0;      % keep a valid record length
end

%% ---------------- Record length ----------------
apHalf  = max(abs(elem_x)) + pitch;
zrec    = max([zmax, max(zs)]) * 1.05;
sigma   = pulseSigma(fc, bandwidth);
maxDel  = max([0, max(cellfun(@(d) max(d(isfinite(d))), {txArr.delays}))]);
Tend    = maxDel + 2*hypot(zrec, 2*apHalf)/c + 8*sigma;
Nt      = max(64, ceil(Tend*fs));
tvec    = (0:Nt-1).'/fs;

%% ---------------- Backend selection ----------------
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
            'SIMUS from MUST failed, switching to the mock backend: %s', ME.message);
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

% Normalise the amplitude so both backends share the same display scale
pk = max(abs(RF(:)));
if pk > 0, RF = RF / pk; end

%% ---------------- Output ----------------
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
%  MUST adapter (all MUST-specific code lives inside this function)
%  ----------------------------------------------------------------------
%  Calling conventions verified against the official MUST documentation
%  (biomecardio.com):
%    RF = SIMUS(X,Z,RC,DELAYS,PARAM)   % 2-D syntax; RF has one column per element
%    PARAM.fc, PARAM.pitch and PARAM.width or kerf are required
%    PARAM.fs defaults to 4*fc; PARAM.bandwidth is the pulse-echo -6 dB
%    fractional bandwidth in %; PARAM.c defaults to 1540;
%    PARAM.TXapodization and PARAM.RXdelay default to none / 0.
%    The time origin is t = 0 (the documented examples build their time axis
%    as t = (0:size(RF,1)-1)/param.fs).
%  If a future MUST release changes the argument order, this is the only
%  place that needs editing.
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
% PARAM.t0 is deliberately not set: it is not listed among the SIMUS
% parameters. The RF time origin is t = 0 and sim_engine reports it as out.t0.

opt = struct('WaitBar', false, 'ParPool', false);

Ntx = numel(txArr);
RFc = cell(1, Ntx);
for k = 1:Ntx
    d = txArr(k).delays;
    a = txArr(k).apod;
    d(~isfinite(d)) = 0;             % inactive elements: delay 0 and apodisation 0
    param.TXapodization = a;
    RFc{k} = simus(xs, zs, rc, d, param, opt);
end

% Transmit events may return different row counts, so zero-pad before stacking
Nrow = max([Nt, cellfun(@(r) size(r,1), RFc)]);
RF = zeros(Nrow, Nel, Ntx);
for k = 1:Ntx
    RF(1:size(RFc{k},1), :, k) = RFc{k};
end
end

%% ======================================================================
%  Mock backend: analytic 2-D cylindrical-wave model, no MUST required.
%  For every scatterer it
%     (1) superimposes the spherical waves radiated by all active transmit
%         elements to build the field at the scatterer, and
%     (2) propagates that field back to each receive element and adds it to
%         the RF matrix.
%  Diffraction and element directivity are not modelled rigorously, but the
%  delay structure is physically exact, which is what matters when verifying
%  the alignment performed by a DAS beamformer.
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

        % --- (1) Transmit: field at the scatterer --------------------
        rt   = hypot(xs(is) - elem_x(act), zs(is) - elem_z(act));   % 1 x Nact
        taut = txd(act) + rt/c;
        ampt = apo(act) .* (zs(is)./rt) ./ sqrt(rt);   % obliquity + cylindrical spread
        sfield = gpulse(tvec - taut, fc, sigma) * ampt(:);          % Nt x 1

        % --- (2) Receive: propagate back to every element ------------
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
%  Helpers
%% ======================================================================
function [d, a, srcxz] = txLaw(elem_x, c, theta, srcDepth, scheme, singleEl, elIdx)
%TXLAW  Transmit delay law. The focal point / virtual source sits on the beam
%  axis at D*(sin theta, cos theta).
Nel = numel(elem_x);
a = ones(1, Nel);
if strcmp(scheme, 'plane') || ~isfinite(srcDepth)
    d = elem_x * sin(theta) / c;                 % plane wave
    srcxz = [NaN NaN];
else
    Px = srcDepth * sin(theta);
    Pz = srcDepth * cos(theta);
    R  = hypot(elem_x - Px, Pz);
    if srcDepth > 0
        d = (max(R) - R) / c;                    % focused (outer elements fire first)
    else
        d = (R - min(R)) / c;                    % diverging (equidistant from source)
    end
    srcxz = [Px Pz];
end

if singleEl
    elIdx = min(max(round(elIdx), 1), Nel);
    a(:) = 0; a(elIdx) = 1;
    d(:) = NaN; d(elIdx) = 0;
    srcxz = [elem_x(elIdx) 0];
end

d = d - min(d(isfinite(d)));                     % normalise so min(active) = 0
d(a == 0) = NaN;                                 % NaN marks inactive elements
end

function s = makeEmptyTx()
s = struct('delays',[],'apod',[],'elem_x',[],'elem_z',[],'t0',0,'c',1540, ...
           'fc',5e6,'fs',20e6,'fnumber',1.5,'rx_apod','rect','scheme','plane', ...
           'angle_deg',0,'focus_mm',NaN,'src_xz',[NaN NaN],'single_element',false);
end

function sigma = pulseSigma(fc, bwPercent)
%PULSESIGMA  Gaussian envelope standard deviation [s] for a given fractional
%  bandwidth [%].
bw = max(bwPercent, 1)/100;
sigma = sqrt(2*log(2)) / (pi * fc * bw);
end

function y = gpulse(t, fc, sigma)
%GPULSE  Gaussian-modulated sinusoidal pulse.
y = exp(-0.5*(t/sigma).^2) .* cos(2*pi*fc*t);
end

function s = getsub(c, name)
if isstruct(c) && isfield(c, name) && isstruct(c.(name)), s = c.(name); else, s = struct(); end
end

function v = getdef(s, name, def)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end
