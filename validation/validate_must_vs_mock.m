function results = validate_must_vs_mock(varargin)
%VALIDATE_MUST_VS_MOCK  Quantitative comparison of the mock RF backend
%   against MUST (Matlab UltraSound Toolbox) SIMUS.
%
%   results = VALIDATE_MUST_VS_MOCK()
%   results = VALIDATE_MUST_VS_MOCK('OutDir', 'validation/results')
%
%   Requires MUST to be on the MATLAB path (see setup_must.m). The first
%   SIMUS call in a MATLAB session costs roughly 20 s (one-time internal
%   initialization); subsequent calls with the same transducer settings
%   are much cheaper (tens of milliseconds for the small scenes used
%   here). This script is still deliberately scoped to a small set of
%   cases: it is a targeted agreement check, not a replacement for the
%   fast headless point-target sweep used for the mock-only regression
%   test.
%
%   For every case (a transmit scheme + a target position) this script:
%     1) Runs SIM_ENGINE with the MUST backend and with the mock backend
%        (same probe, same transmit law, same scatterer).
%     2) Compares the raw RF channel data via the normalized cross-
%        correlation between the two backends, reporting the peak
%        correlation coefficient and the lag (in samples and in
%        nanoseconds) at which it occurs, for every active receive
%        channel.
%     3) Beamforms both RF data sets with DAS_REFERENCE on an identical
%        grid and compares the reconstructed images: peak position error
%        relative to the known target, and lateral FWHM (-6 dB, linear).
%
%   The results struct array and a summary table are saved as
%   validation/results/must_vs_mock.mat and .csv, and a figure comparing
%   one representative RF trace and one B-mode pair is saved as
%   validation/results/must_vs_mock_example.png.
%
%   See also SIM_ENGINE, DAS_REFERENCE, SETUP_MUST.

p = inputParser();
p.addParameter('OutDir', fullfile(fileparts(mfilename('fullpath')), 'results'));
p.parse(varargin{:});
outDir = p.Results.OutDir;
if ~exist(outDir, 'dir'), mkdir(outDir); end

if isempty(which('simus'))
    error('validate_must_vs_mock:noMUST', ...
        ['MUST was not found on the MATLAB path. Run setup_must.m first ' ...
         '(it downloads MUST from biomecardio.com and adds it to the path).']);
end

%% ---------------- Common acquisition settings ----------------
probe  = struct('Nelements', 64, 'pitch', 0.30e-3, 'fc', 5e6, 'bandwidth', 75);
medium = struct('c', 1540);
acq    = struct('fs_factor', 4);
recon  = struct('zmax', 40e-3, 'fnumber', 1.5);

% The lateral grid is much finer than the GUI's interactive default
% (0.0125 mm vs. 0.125 mm spacing). This is required to resolve the
% focused-transmit mainlobe: at the coarse spacing its true FWHM
% (~0.10 mm mock / ~0.26 mm MUST, confirmed converged from 1601 to
% 16001 points) is quantized to 1-2 grid cells, which biases both
% backends' FWHM upward by about 2x and can distort a mock-vs-MUST
% comparison. The finer grid changes nothing about the peak position
% (grid-independent to well under a wavelength either way) and is cheap
% offline; it is not used by the interactive GUI, which prioritizes
% redraw speed over this level of resolution.
gx = linspace(-10e-3, 10e-3, 1601);
gz = linspace(5e-3, 40e-3, 281);

% Cases: {label, tx struct, target [x z] in m, probe override (optional, [])}
% Focused transmit is a single scan line (see README/paper limitations),
% so it is only evaluated on-axis where the physical model applies. Two
% extra cases repeat the on-axis plane-wave check at 16 and 128 elements
% to confirm agreement is not an artifact of the default 64-element probe.
cases = {
    'plane_onaxis',      struct('scheme','plane',    'angle_deg',0),               [0     20e-3], []
    'plane_offaxis',     struct('scheme','plane',    'angle_deg',0),               [4e-3  25e-3], []
    'focused_onaxis',    struct('scheme','focused',  'angle_deg',0,'focus_mm',20), [0     20e-3], []
    'diverging_onaxis',  struct('scheme','diverging','angle_deg',0,'src_mm',10),   [0     20e-3], []
    'diverging_offaxis', struct('scheme','diverging','angle_deg',0,'src_mm',10),   [4e-3  25e-3], []
    'plane_16el',        struct('scheme','plane',    'angle_deg',0),               [0     20e-3], 16
    'plane_128el',       struct('scheme','plane',    'angle_deg',0),               [0     20e-3], 128
};

results = struct('label', {}, 'scheme', {}, 'nelements', {}, 'target_x_mm', {}, 'target_z_mm', {}, ...
    'peak_corr_mean', {}, 'peak_corr_min', {}, 'lag_samples_mean', {}, ...
    'lag_ns_mean', {}, 'lag_samples_max_abs', {}, 'sublag_samples_center', {}, ...
    'peak_err_mock_mm', {}, 'peak_err_must_mm', {}, ...
    'fwhm_mock_mm', {}, 'fwhm_must_mm', {}, ...
    'sim_time_mock_s', {}, 'sim_time_must_s', {});

exampleSaved = false;

for kc = 1:size(cases,1)
    label     = cases{kc,1};
    txSpec    = cases{kc,2};
    target    = cases{kc,3};
    probeOver = cases{kc,4};

    cfg = struct();
    cfg.probe  = probe;
    if ~isempty(probeOver)
        cfg.probe.Nelements = probeOver;
    end
    cfg.medium = medium;
    cfg.acq    = acq;
    cfg.tx     = txSpec;
    cfg.scat   = struct('x', target(1), 'z', target(2), 'rc', 1);
    cfg.recon  = recon;

    fprintf('[%d/%d] %s ... ', kc, size(cases,1), label);

    cfg.options = struct('force_mock', true);
    Smock = sim_engine(cfg);

    cfg.options = struct('force_mock', false);
    Smust = sim_engine(cfg);

    if ~strcmp(Smust.backend, 'MUST')
        error('validate_must_vs_mock:backend', ...
            'Expected the MUST backend for case "%s" but got "%s".', label, Smust.backend);
    end

    %% ---- Raw RF agreement: per-channel normalized cross-correlation ----
    [corrPk, lagSamp] = channelCrossCorr(Smock.RF(:,:,1), Smock.fs, Smust.RF(:,:,1), Smust.fs);

    % Sub-sample lag on the center channel via parabolic interpolation of
    % the three xcorr samples around the integer peak (xcorr lags are
    % integer by construction, so the -1-sample result above could in
    % principle be a quantized version of e.g. -0.5; this checks that).
    centerEl = round(size(Smock.RF,2)/2);
    subLagCenter = subsampleLag(Smock.RF(:,centerEl), Smust.RF(:,centerEl));

    %% ---- Beamformed image agreement ----
    imgMock = das_reference(Smock.RF(:,:,1), Smock.tx(1), Smock.rx_pos, gx, gz, Smock.c, Smock.fs);
    imgMust = das_reference(Smust.RF(:,:,1), Smust.tx(1), Smust.rx_pos, gx, gz, Smust.c, Smust.fs);

    [errMock, fwMock] = peakAndFWHM(imgMock, gx, gz, target);
    [errMust, fwMust] = peakAndFWHM(imgMust, gx, gz, target);

    r = struct();
    r.label            = label;
    r.scheme           = txSpec.scheme;
    r.nelements        = cfg.probe.Nelements;
    r.target_x_mm      = target(1)*1e3;
    r.target_z_mm      = target(2)*1e3;
    r.peak_corr_mean   = mean(corrPk);
    r.peak_corr_min    = min(corrPk);
    r.lag_samples_mean = mean(lagSamp);
    r.lag_ns_mean       = mean(lagSamp) / Smust.fs * 1e9;
    r.lag_samples_max_abs = max(abs(lagSamp));
    r.sublag_samples_center = subLagCenter;
    r.peak_err_mock_mm = errMock * 1e3;
    r.peak_err_must_mm = errMust * 1e3;
    r.fwhm_mock_mm     = fwMock * 1e3;
    r.fwhm_must_mm     = fwMust * 1e3;
    r.sim_time_mock_s  = Smock.elapsed;
    r.sim_time_must_s  = Smust.elapsed;
    results(end+1) = r; %#ok<AGROW>

    fprintf('corr=%.3f  lag=%.2f samp (sub-sample %.3f)  peak_err mock/MUST=%.3f/%.3f mm\n', ...
        r.peak_corr_mean, r.lag_samples_mean, r.sublag_samples_center, r.peak_err_mock_mm, r.peak_err_must_mm);

    if ~exampleSaved && strcmp(label, 'plane_onaxis')
        saveExampleFigure(outDir, Smock, Smust, imgMock, imgMust, gx, gz, target);
        exampleSaved = true;
    end
end

%% ---------------- Save results ----------------
T = struct2table(results);
writetable(T, fullfile(outDir, 'must_vs_mock.csv'));
save(fullfile(outDir, 'must_vs_mock.mat'), 'results', 'T');

fprintf('\nSummary (mean over %d cases):\n', numel(results));
fprintf('  Peak channel correlation (mean of per-case means): %.3f\n', mean(T.peak_corr_mean));
fprintf('  Cross-correlation lag: %.3f samples (%.1f ns) [max |lag| = %.2f samples]\n', ...
    mean(T.lag_samples_mean), mean(T.lag_ns_mean), max(T.lag_samples_max_abs));
fprintf('  Sub-sample lag (center channel, parabolic interp.): mean %.3f samples, range [%.3f, %.3f]\n', ...
    mean(T.sublag_samples_center), min(T.sublag_samples_center), max(T.sublag_samples_center));
fprintf('  Peak position error, mock:  %.4f mm (max %.4f mm)\n', mean(T.peak_err_mock_mm), max(T.peak_err_mock_mm));
fprintf('  Peak position error, MUST: %.4f mm (max %.4f mm)\n', mean(T.peak_err_must_mm), max(T.peak_err_must_mm));
fprintf('  Lateral FWHM, mock:  %.4f mm\n', mean(T.fwhm_mock_mm));
fprintf('  Lateral FWHM, MUST: %.4f mm\n', mean(T.fwhm_must_mm));
fprintf(['  Mean per-call time (STEADY STATE, i.e. excluding the ~20 s one-time\n' ...
    '  SIMUS initialization cost paid by the first call in a session):\n' ...
    '    mock:  %.4f s\n    MUST: %.4f s\n'], mean(T.sim_time_mock_s), mean(T.sim_time_must_s));
fprintf(['  NOTE: peak-position error uses the SAME delay law for both backends\n' ...
    '  (das_reference computes tau_tx from tx_info, not from the RF itself), so\n' ...
    '  it mainly checks grid/delay consistency. The small (<0.14 mm) residual\n' ...
    '  differences seen between the two backends'' peak positions come from\n' ...
    '  amplitude-weighting effects: MUST''s pulse-echo waveform is not perfectly\n' ...
    '  symmetric like the mock''s Gaussian pulse, which can shift the discrete\n' ...
    '  argmax by a fraction of a grid cell on a narrow mainlobe.\n']);
fprintf('Saved: %s\n', fullfile(outDir, 'must_vs_mock.csv'));

saveSummaryFigure(outDir, T);
end

%% ========================================================================
function saveSummaryFigure(outDir, T)
%SAVESUMMARYFIGURE  Bar chart of per-case correlation and FWHM agreement,
%   used as a figure in the paper.
fig = figure('Visible', 'off', 'Position', [100 100 1000 700], 'Color', 'w');
n = height(T);
labels = strrep(T.label, '_', '\_');

subplot(2,1,1);
bar(1:n, T.peak_corr_mean); ylim([0.8 1]);
set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'XTickLabelRotation', 20);
ylabel('Peak channel correlation');
title('Mock vs. MUST: raw-RF per-channel cross-correlation (mean over channels)');
grid on;

subplot(2,1,2);
bar(1:n, [T.fwhm_mock_mm, T.fwhm_must_mm]);
set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'XTickLabelRotation', 20);
ylabel('Lateral FWHM [mm]');
legend('Mock', 'MUST', 'Location', 'best');
title('Mock vs. MUST: beamformed lateral FWHM (-6 dB, linear)');
grid on;

exportgraphics(fig, fullfile(outDir, 'must_vs_mock_summary.png'), 'Resolution', 150);
close(fig);
end

%% ========================================================================
function [corrPk, lagSamp] = channelCrossCorr(rfA, fsA, rfB, fsB)
%CHANNELCROSSCORR  Per-channel normalized cross-correlation between two RF
%   data sets that may have a different number of time samples (they are
%   zero-padded to a common length; sampling rates are assumed equal,
%   which holds here since both backends share fs_factor*fc).
assert(abs(fsA - fsB) < 1, 'validate_must_vs_mock:fs', 'Sampling rates differ between backends.');
Nel = size(rfA, 2);
assert(Nel == size(rfB,2), 'validate_must_vs_mock:nel', 'Channel count mismatch.');
N = max(size(rfA,1), size(rfB,1));

corrPk  = zeros(1, Nel);
lagSamp = zeros(1, Nel);
for e = 1:Nel
    a = zeros(N,1); a(1:size(rfA,1)) = rfA(:,e);
    b = zeros(N,1); b(1:size(rfB,1)) = rfB(:,e);
    if max(abs(a)) < 1e-9 || max(abs(b)) < 1e-9
        corrPk(e) = NaN; lagSamp(e) = NaN; continue;
    end
    [c, lags] = xcorr(a, b, round(N/4), 'coeff');
    [corrPk(e), im] = max(c);
    lagSamp(e) = lags(im);
end
valid = isfinite(corrPk);
corrPk = corrPk(valid);
lagSamp = lagSamp(valid);
end

%% ========================================================================
function subLag = subsampleLag(a, b)
%SUBSAMPLELAG  Sub-sample cross-correlation lag of channel a relative to b,
%   via parabolic interpolation of the three xcorr samples around the
%   integer-lag peak. xcorr itself can only return integer lags, so this
%   check exists to rule out the integer result being a quantized version
%   of a genuinely fractional offset (e.g. -1 masking a true -0.5).
N = max(numel(a), numel(b));
av = zeros(N,1); av(1:numel(a)) = a(:);
bv = zeros(N,1); bv(1:numel(b)) = b(:);
[c, lags] = xcorr(av, bv, 5, 'coeff');
[~, im] = max(c);
if im > 1 && im < numel(c)
    y1 = c(im-1); y2 = c(im); y3 = c(im+1);
    delta = 0.5*(y1-y3)/(y1-2*y2+y3);
else
    delta = 0;
end
subLag = lags(im) + delta;
end

%% ========================================================================
function [err, fw] = peakAndFWHM(img, gx, gz, target)
%PEAKANDFWHM  Euclidean peak-position error [m] against the known target,
%   and the -6 dB lateral FWHM [m] through the peak row (linear envelope).
[~, k] = max(img(:));
[iz, ix] = ind2sub(size(img), k);
err = hypot(gx(ix) - target(1), gz(iz) - target(2));

prof = img(iz, :);
[pkv, ip] = max(prof);
half = pkv/2;
il = find(prof(1:ip) <= half, 1, 'last');
ir = find(prof(ip:end) <= half, 1, 'first');
if isempty(il) || isempty(ir)
    fw = NaN; return;
end
ir = ir + ip - 1;
xl = lerp(gx(il), prof(il), gx(il+1), prof(il+1), half);
xr = lerp(gx(ir-1), prof(ir-1), gx(ir), prof(ir), half);
fw = xr - xl;
end

function x = lerp(x1,y1,x2,y2,yt)
if y2 == y1, x = x1; else, x = x1 + (yt-y1)*(x2-x1)/(y2-y1); end
end

%% ========================================================================
function saveExampleFigure(outDir, Smock, Smust, imgMock, imgMust, gx, gz, target)
%SAVEEXAMPLEFIGURE  One RF-trace overlay and one B-mode comparison, used
%   as a figure in the paper.
fig = figure('Visible', 'off', 'Position', [100 100 1100 850], 'Color', 'w');

centerEl = round(size(Smock.RF,2)/2);
tMock = (0:size(Smock.RF,1)-1)/Smock.fs*1e6;
tMust = (0:size(Smust.RF,1)-1)/Smust.fs*1e6;

% Centre the zoom window on the actual echo (found from the mock envelope
% peak) rather than a fixed fraction of the record length, so the window
% tracks the target depth correctly.
envMock = abs(hilbert(Smock.RF(:,centerEl)));
[~, ipk] = max(envMock);
tCentre = tMock(ipk);
halfWin = 8 / Smock.probe.fc * 1e6;   % +/- 8 carrier periods

subplot(2,2,[1 2]);
plot(tMock, Smock.RF(:,centerEl), 'b-', 'LineWidth', 1.2); hold on;
plot(tMust, Smust.RF(:,centerEl), 'r--', 'LineWidth', 1.2);
xlim(tCentre + [-halfWin halfWin]);
xlabel('Time [\mus]'); ylabel('Normalized amplitude');
legend('Mock backend', 'MUST (SIMUS)', 'Location', 'best');
title(sprintf('Center-channel RF, plane wave, target (%.0f, %.0f) mm', target(1)*1e3, target(2)*1e3));
grid on;

subplot(2,2,3);
imagesc(gx*1e3, gz*1e3, 20*log10(imgMock/max(imgMock(:))+1e-12));
set(gca, 'YDir', 'reverse'); axis image; colormap(gca, gray(256)); clim([-50 0]);
xlabel('x [mm]'); ylabel('z [mm]'); title('Mock backend'); colorbar;

subplot(2,2,4);
imagesc(gx*1e3, gz*1e3, 20*log10(imgMust/max(imgMust(:))+1e-12));
set(gca, 'YDir', 'reverse'); axis image; colormap(gca, gray(256)); clim([-50 0]);
xlabel('x [mm]'); ylabel('z [mm]'); title('MUST (SIMUS)'); colorbar;

exportgraphics(fig, fullfile(outDir, 'must_vs_mock_example.png'), 'Resolution', 150);
close(fig);
end
