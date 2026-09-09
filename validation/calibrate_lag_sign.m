function calibrate_lag_sign()
%CALIBRATE_LAG_SIGN  Empirically pins down the sign convention of the
%   cross-correlation lag reported by VALIDATE_MUST_VS_MOCK, using a
%   synthetic pair of signals with a KNOWN, hand-constructed delay.
%
%   VALIDATE_MUST_VS_MOCK calls its internal helper as
%   channelCrossCorr(Smock.RF, Smock.fs, Smust.RF, Smust.fs), i.e. the
%   first signal argument is always the mock backend and the second is
%   always MUST. This script builds a synthetic pulse "a" and a second
%   pulse "b" that is a copy of "a" delayed by a known number of samples
%   (b arrives LATER than a), reproduces the exact same cross-correlation
%   and parabolic sub-sample-refinement logic used there, and reports
%   which sign of the returned lag corresponds to "a is earlier than b".
%
%   Do not rely on reasoning about MATLAB's xcorr lag convention in the
%   abstract: run this once and read the printed conclusion, because that
%   is the convention actually used to interpret the real mock-vs-MUST
%   results in validate_must_vs_mock.m and reported in the paper.

N = 200;
t = (0:N-1)';
a = exp(-0.5*((t-40)/5).^2) .* cos(2*pi*0.15*t);   % pulse "a", peaks near sample 40
knownDelay = 2;                                     % ground truth: b arrives LATER by this many samples
b = [zeros(knownDelay,1); a(1:end-knownDelay)];     % "b" = a, shifted later by knownDelay samples

[corrPk, lagSamp] = channelCrossCorr(a, 1, b, 1);
subLag = subsampleLag(a, b);

fprintf('Ground truth: b (2nd argument) arrives LATER than a (1st argument) by +%d samples.\n', knownDelay);
fprintf('channelCrossCorr(a, fs, b, fs) reports lag = %+.2f samples (corr=%.3f)\n', lagSamp, corrPk);
fprintf('subsampleLag(a, b)             reports lag = %+.4f samples\n', subLag);
fprintf('\n');
if lagSamp < 0
    fprintf(['CONCLUSION: a NEGATIVE reported lag means the FIRST argument (mock, in the\n' ...
        'real comparison) is EARLIER than the SECOND argument (MUST).\n']);
else
    fprintf(['CONCLUSION: a POSITIVE reported lag means the FIRST argument (mock, in the\n' ...
        'real comparison) is EARLIER than the SECOND argument (MUST).\n']);
end
end

%% ---- verbatim copies of the local functions under test, kept in sync
%      with validate_must_vs_mock.m by inspection (both are short). ----
function [corrPk, lagSamp] = channelCrossCorr(rfA, fsA, rfB, fsB)
assert(abs(fsA - fsB) < 1e-9);
Nel = size(rfA, 2);
N = max(size(rfA,1), size(rfB,1));
corrPk  = zeros(1, Nel);
lagSamp = zeros(1, Nel);
for e = 1:Nel
    a = zeros(N,1); a(1:size(rfA,1)) = rfA(:,e);
    b = zeros(N,1); b(1:size(rfB,1)) = rfB(:,e);
    [c, lags] = xcorr(a, b, round(N/4), 'coeff');
    [corrPk(e), im] = max(c);
    lagSamp(e) = lags(im);
end
end

function subLag = subsampleLag(a, b)
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
