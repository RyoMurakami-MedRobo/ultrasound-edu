function dump_reference(varargin)
%DUMP_REFERENCE  Generate golden fixtures for the Python parity tests.
%
%   dump_reference()
%   dump_reference('OutDir', '<repo>/python/parity/fixtures')
%
%   Runs SIM_ENGINE (mock backend, force_mock = true so the result is
%   deterministic and MUST-independent), DAS_REFERENCE and
%   DAS_CUSTOM_TEMPLATE on a small matrix of configurations that exercises
%   every branch the Python port has to reproduce, and saves one .mat file
%   per case under python/parity/fixtures/.
%
%   The Python side (python/tests/test_parity.py) loads these, rebuilds the
%   same configuration, and asserts agreement to rtol = 1e-9. Regenerate
%   whenever the numerics of sim_engine.m / das_reference.m /
%   das_custom_template.m change (see CLAUDE.md).
%
%   Files are written with save('-v7') on purpose: '-v7.3' is HDF5 and
%   scipy.io.loadmat cannot read it.

p = inputParser();
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(here));
p.addParameter('OutDir', fullfile(here, 'fixtures'));
p.parse(varargin{:});
outDir = p.Results.OutDir;
if ~exist(outDir, 'dir'), mkdir(outDir); end
addpath(repoRoot);

% Common probe / grid (small, to bound the record length and fixture size)
base = struct('n_elements', 32, 'pitch', 0.30e-3, 'fc', 5e6, 'bandwidth', 75, ...
              'c', 1540, 'fs_factor', 4, 'scheme', 'plane', 'angle_deg', 0, ...
              'focus_mm', 18, 'src_mm', 10, 'single_element', false, ...
              'element_index', 16, 'scat_x', 0, 'scat_z', 18e-3, 'scat_rc', 1, ...
              'zmax', 26e-3, 'fnumber', 1.5, 'rx_apod', 'rect', ...
              'nx', 45, 'nz', 61, 'x_half_mm', 10, 'z_min_mm', 3);

C = {};
C{end+1} = merge(base, 'label', 'plane_onaxis');
C{end+1} = merge(base, 'label', 'plane_steered_offaxis', 'angle_deg', 10, ...
                 'scat_x', 3e-3, 'scat_z', 16e-3);
C{end+1} = merge(base, 'label', 'focused_onaxis', 'scheme', 'focused', 'focus_mm', 18);
C{end+1} = merge(base, 'label', 'diverging_offaxis_fullap', 'scheme', 'diverging', ...
                 'src_mm', 10, 'scat_x', 3e-3, 'scat_z', 16e-3, 'fnumber', 0);
C{end+1} = merge(base, 'label', 'single_element', 'single_element', true, ...
                 'element_index', 12, 'scat_z', 15e-3);
C{end+1} = merge(base, 'label', 'hann_apod_multiscat', 'rx_apod', 'hann', ...
                 'scat_x', [2e-3 -3e-3], 'scat_z', [12e-3 20e-3], 'scat_rc', [1 0.6]);
C{end+1} = merge(base, 'label', 'fnum0_multiscat', 'fnumber', 0, ...
                 'scat_x', [0 -4e-3 5e-3], 'scat_z', [10e-3 16e-3 22e-3], 'scat_rc', [1 0.8 0.5]);
C{end+1} = merge(base, 'label', 'multifocus', 'scheme', 'focused', ...
                 'focus_mm', [12 18 24], ...
                 'scat_x', [0 0 0], 'scat_z', [12e-3 18e-3 24e-3], 'scat_rc', [1 1 1]);

for i = 1:numel(C)
    dumpOne(C{i}, outDir);
end
fprintf('Wrote %d fixtures to %s\n', numel(C), outDir);
end

% ------------------------------------------------------------------------
function s = merge(s, varargin)
for k = 1:2:numel(varargin)
    s.(varargin{k}) = varargin{k+1};
end
end

% ------------------------------------------------------------------------
function dumpOne(c, outDir)
cfg = struct();
cfg.probe   = struct('Nelements', c.n_elements, 'pitch', c.pitch, 'fc', c.fc, ...
                     'bandwidth', c.bandwidth);
cfg.medium  = struct('c', c.c);
cfg.acq     = struct('fs_factor', c.fs_factor);
cfg.tx      = struct('scheme', c.scheme, 'angle_deg', c.angle_deg, ...
                     'focus_mm', c.focus_mm, 'src_mm', c.src_mm, ...
                     'single_element', c.single_element, 'element_index', c.element_index);
cfg.scat    = struct('x', c.scat_x, 'z', c.scat_z, 'rc', c.scat_rc);
cfg.recon   = struct('zmax', c.zmax, 'fnumber', c.fnumber, 'rx_apod', c.rx_apod);
cfg.options = struct('force_mock', true);

S = sim_engine(cfg);

% Push the recon f-number / apodisation into every event, like onBeamform()
for k = 1:numel(S.tx)
    S.tx(k).fnumber = c.fnumber;
    S.tx(k).rx_apod = c.rx_apod;
end

xh = c.x_half_mm * 1e-3;
gx = linspace(-xh, xh, c.nx);
z0 = min(c.z_min_mm, c.zmax*1e3 - 1) * 1e-3;
gz = linspace(max(z0, 1e-4), c.zmax, c.nz);

ntx = numel(S.tx);
[Nt, Nel, ~] = size(S.RF);

% Stack tx info
tx_delays  = nan(Nel, ntx);
tx_apod    = zeros(Nel, ntx);
tx_focus_mm = zeros(1, ntx);
tx_src_xz  = zeros(2, ntx);
for k = 1:ntx
    tx_delays(:, k)  = S.tx(k).delays(:);
    tx_apod(:, k)    = S.tx(k).apod(:);
    tx_focus_mm(k)   = S.tx(k).focus_mm;
    tx_src_xz(:, k)  = S.tx(k).src_xz(:);
end

% Beamform: single event directly, multi-focus composited by depth zone
if ntx == 1
    [img_ref, delays_ref]       = das_reference(S.RF(:,:,1), S.tx(1), S.rx_pos, gx, gz, S.c, S.fs);
    [img_custom, delays_custom] = das_custom_template(S.RF(:,:,1), S.tx(1), S.rx_pos, gx, gz, S.c, S.fs);
    has_delays = 1;
else
    img_ref    = compositeAllTx(@das_reference, S, gx, gz);
    img_custom = compositeAllTx(@das_custom_template, S, gx, gz);
    delays_ref = []; delays_custom = [];
    has_delays = 0;
end

% Metrics: A = reference, B = custom (mirrors updateCompare)
met = compareMetrics(img_ref, img_custom, gx, gz);

save(fullfile(outDir, [c.label '.mat']), '-v7', ...
    'c', 'gx', 'gz', 'ntx', 'Nt', 'Nel', ...
    'tx_delays', 'tx_apod', 'tx_focus_mm', 'tx_src_xz', ...
    'has_delays', 'img_ref', 'img_custom', 'delays_ref', 'delays_custom', 'met');

% Save RF / t0 / fs separately-named for clarity
RF = S.RF; t0 = S.t0; fs = S.fs; rx_pos = S.rx_pos; backend = S.backend;
save(fullfile(outDir, [c.label '.mat']), '-append', ...
    'RF', 't0', 'fs', 'rx_pos', 'backend');

fprintf('  %-26s  RF %dx%dx%d  img %dx%d\n', c.label, Nt, Nel, ntx, size(img_ref,1), size(img_ref,2));
end

% ------------------------------------------------------------------------
function img = compositeAllTx(fh, S, gx, gz)
%COMPOSITEALLTX  Copy of runAllTx from main_gui.m (multi-focus depth zones).
Ntx = numel(S.tx);
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

% ------------------------------------------------------------------------
function met = compareMetrics(A, B, gx, gz)
nA = max(A(:)); if nA <= 0, nA = 1; end
An = A / nA; Bn = B / nA;
[fwA, ~, pkA] = latMet(An, gx, gz, []);
[fwB, ~, pkB] = latMet(Bn, gx, gz, pkA(2));
met = struct();
met.fwhm_a_mm = fwA * 1e3;
met.fwhm_b_mm = fwB * 1e3;
met.peak_a_mm = pkA * 1e3;
met.peak_b_mm = pkB * 1e3;
met.mse = mean((An(:) - Bn(:)).^2);
met.max_abs_err = max(abs(An(:) - Bn(:)));
end

function [fw, prof, pk] = latMet(img, gx, gz, forceZ)
if ~isempty(forceZ) && isfinite(forceZ)
    [~, iz] = min(abs(gz - forceZ));
    [~, ix] = max(img(iz, :));
else
    [~, k] = max(img(:));
    [iz, ix] = ind2sub(size(img), k);
end
prof = img(iz, :);
pk = [gx(ix), gz(iz)];
fw = fwhmLin(prof, gx);
end

function w = fwhmLin(prof, gx)
[pkv, ip] = max(prof);
if pkv <= 0, w = NaN; return; end
half = pkv / 2;
il = find(prof(1:ip) <= half, 1, 'last');
ir = find(prof(ip:end) <= half, 1, 'first');
if isempty(il) || isempty(ir), w = NaN; return; end
ir = ir + ip - 1;
xl = cross1(gx(il), prof(il), gx(il+1), prof(il+1), half);
xr = cross1(gx(ir-1), prof(ir-1), gx(ir), prof(ir), half);
w = xr - xl;
end

function x = cross1(x1, y1, x2, y2, yt)
if y2 == y1, x = x1; else, x = x1 + (yt - y1) * (x2 - x1) / (y2 - y1); end
end
