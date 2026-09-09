function [bmode_img, delays] = das_custom_template(rf_data, tx_info, rx_pos, grid_x, grid_z, sound_speed, fs)
%DAS_CUSTOM_TEMPLATE  Skeleton for writing your own DAS beamformer.
%
%   [bmode_img, delays] = DAS_CUSTOM_TEMPLATE(rf_data, tx_info, rx_pos, ...
%                                             grid_x, grid_z, sound_speed, fs)
%
%   HOW TO USE
%     1) Copy this file under a name of your choice (e.g. my_das.m).
%     2) Rename the function so it matches the file name.
%     3) Edit [STEP 1] ... [STEP 5] below.
%     4) Type that name into the "Custom function" field of MAIN_GUI and pick
%        "Custom DAS" as algorithm A or B to compare it with the reference.
%
%   CONTRACT (breaking it breaks the comparison and the metrics)
%     - Keep the argument order and meaning identical to DAS_REFERENCE.
%     - rf_data is [Nt x Nel] (time along the first dimension, columns are
%       receive elements).
%     - bmode_img must be [numel(grid_z) x numel(grid_x)] and must be the
%       **linear envelope** (no log compression). A wrong size makes the
%       difference image and the MSE fail.
%     - delays must be [Npix x Nel] holding the **total two-way time [s]**,
%       with NaN wherever a sample is not summed. Skip it when nargout < 2.
%     - Pixels are ordered column-major (the linear index after reshape).
%
%   INITIAL STATE
%     A deliberately naive DAS (nearest-neighbour sampling plus a rectangular
%     apodisation) so the file runs out of the box. Comparing it with the
%     reference (linear interpolation) shows the interpolation error as a
%     raised sidelobe floor in the difference image and the lateral profile.
%
%   See also DAS_REFERENCE, SIM_ENGINE, MAIN_GUI.

narginchk(7, 7);

%% =====================================================================
%  [STEP 0] Input handling (usually left untouched)
%% =====================================================================
if ndims(rf_data) > 2 %#ok<ISMAT>
    rf_data = rf_data(:,:,1);
end
[Nt, Nel] = size(rf_data);

rxx = rx_pos(:).';                       % receive element x coordinates [m]
rxz = zeros(1, Nel);                     % linear array, so z = 0
assert(numel(rxx) == Nel, 'rx_pos does not match the number of RF columns.');

if isvector(grid_x) && isvector(grid_z)
    [XI, ZI] = meshgrid(grid_x(:).', grid_z(:).');
else
    XI = grid_x; ZI = grid_z;
end
imgSize = size(XI);
xg = XI(:); zg = ZI(:); Npix = numel(xg);

c  = sound_speed;
t0 = tx_info.t0;                         % time of the first RF sample [s]

%% =====================================================================
%  [STEP 1] Transmit arrival time tau_tx(p)   <<< EDIT HERE
%  ---------------------------------------------------------------------
%  Two models are provided, mirroring the reference implementation.
%   (a) First-arrival (Huygens) model - exact for plane, diverging and
%       single-element transmits:
%         tau_tx(p) = min_e ( delays(e) + |p - elem_e| / c )
%   (b) Virtual-source model - used for focused transmits:
%         tau_tx(p) = T_F + s * |p - F| / c   (s = -1 before, +1 beyond focus)
%  Using (a) for a focused transmit makes the beamformer latch onto the weak
%  edge wave beyond the focus and shifts the image by roughly 1 mm in depth.
%  Removing the distinction is an instructive experiment in itself.
%% =====================================================================
ex  = tx_info.elem_x;
ez  = tx_info.elem_z;
td  = tx_info.delays;
ap  = tx_info.apod;
act = isfinite(td) & ap(:).' ~= 0;       % only the elements that actually fire
exA = ex(act); ezA = ez(act); tdA = td(act);

src = tx_info.src_xz;                    % focal point / virtual source [x z], NaN if none
useVirtualSource = strcmpi(tx_info.scheme,'focused') && all(isfinite(src)) && ...
                   src(2) > 0 && ~tx_info.single_element;

if useVirtualSource
    Rf   = hypot(exA - src(1), ezA - src(2));
    Tf   = mean(tdA + Rf / c);           % instant the wavefront converges on the focus [s]
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
%  [STEP 2] Receive aperture weights (f-number, apodisation)  <<< EDITABLE
%% =====================================================================
fnum = tx_info.fnumber;
if ~isfinite(fnum) || fnum <= 0
    halfAp = inf(Npix, 1);               % full aperture
else
    halfAp = zg ./ (2*fnum);             % aperture grows with depth
end

%% =====================================================================
%  [STEP 3] Delay-and-sum main loop   <<< EDIT HERE
%  ---------------------------------------------------------------------
%  The initial version uses nearest-neighbour sampling. Replacing round()
%  with floor() plus linear interpolation should reproduce the reference
%  implementation exactly, which makes a handy self-check.
%% =====================================================================
bf = zeros(Npix, 1);
if nargout > 1
    delays = nan(Npix, Nel);
end

for e = 1:Nel
    dx  = xg - rxx(e);
    tau = tau_tx + hypot(dx, zg - rxz(e)) / c;    % total two-way time [s]

    w   = double(abs(dx) <= halfAp);              % rectangular apodisation

    idx = round((tau - t0) * fs + 1);             % <-- nearest-neighbour sampling
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
%  [STEP 4] Envelope detection   <<< EDITABLE (try something other than Hilbert)
%% =====================================================================
bf2 = reshape(bf, imgSize);
if imgSize(1) >= 4
    bmode_img = abs(hilbert(bf2));
else
    bmode_img = abs(bf2);
end

%% =====================================================================
%  [STEP 5] Post-processing (TGC, speckle reduction, ...) goes here
%% =====================================================================
% bmode_img = ...;

end
