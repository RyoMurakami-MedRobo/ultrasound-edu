"""Comparison metrics for two beamformed images.

Python port of the metric helpers in ``main_gui.m`` (``lateralMetrics``,
``fwhm``, ``crossing``, ``todB``, and the ``updateCompare`` math). MATLAB is the
source of truth.

Conventions (see the README "Definition of the metrics"):
  * B-mode images are linear envelopes; ``20*log10`` is display-only.
  * The difference image and the MSE normalise **both** images by ``max(A)``.
  * FWHM is the -6 dB (half-amplitude) full width of the linear envelope, with
    the crossings located by linear interpolation.
"""

from __future__ import annotations

import numpy as np

__all__ = ["fwhm", "lateral_metrics", "compare_images", "to_db"]


def _crossing(x1, y1, x2, y2, yt):
    if y2 == y1:
        return x1
    return x1 + (yt - y1) * (x2 - x1) / (y2 - y1)


def fwhm(prof: np.ndarray, gx: np.ndarray) -> float:
    """Full width at half maximum (= -6 dB) of a linear envelope profile.
    Crossings are located by linear interpolation. Returns NaN when the
    half-maximum level is not crossed on both sides."""
    prof = np.asarray(prof, dtype=float)
    gx = np.asarray(gx, dtype=float)
    ip = int(np.argmax(prof))
    pkv = prof[ip]
    if pkv <= 0:
        return np.nan
    half = pkv / 2.0

    left = np.where(prof[: ip + 1] <= half)[0]      # MATLAB: find(prof(1:ip)<=half,1,'last')
    right = np.where(prof[ip:] <= half)[0]          # MATLAB: find(prof(ip:end)<=half,1,'first')
    if left.size == 0 or right.size == 0:
        return np.nan
    il = int(left[-1])
    ir = int(right[0]) + ip

    xl = _crossing(gx[il], prof[il], gx[il + 1], prof[il + 1], half)
    xr = _crossing(gx[ir - 1], prof[ir - 1], gx[ir], prof[ir], half)
    return xr - xl


def lateral_metrics(img, gx, gz, force_z=None):
    """Lateral profile through the peak row and its -6 dB width, on linear
    values. Returns ``(fw, prof, (peak_x, peak_z))``.

    ``force_z`` pins the analysis row to a depth (used so profile B is taken at
    the same depth as profile A)."""
    img = np.asarray(img, dtype=float)
    gx = np.asarray(gx, dtype=float)
    gz = np.asarray(gz, dtype=float)

    if force_z is not None and np.isfinite(force_z):
        iz = int(np.argmin(np.abs(gz - force_z)))
        ix = int(np.argmax(img[iz, :]))
    else:
        iz, ix = np.unravel_index(np.argmax(img), img.shape)

    prof = img[iz, :]
    return fwhm(prof, gx), prof, (gx[ix], gz[iz])


def to_db(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    m = np.max(p)
    if m <= 0:
        m = 1.0
    return 20.0 * np.log10(p / m + 1e-12)


def compare_images(img_a, img_b, gx, gz):
    """Full comparison of algorithm A vs. B, mirroring ``updateCompare``.

    Returns a dict with ``fwhm_a_mm``, ``fwhm_b_mm``, ``peak_a_mm``,
    ``peak_b_mm``, ``mse``, ``max_abs_err``, ``diff_db`` (image, or None),
    ``prof_a``, ``prof_b``.
    """
    A = np.asarray(img_a, dtype=float)
    gx = np.asarray(gx, dtype=float)
    gz = np.asarray(gz, dtype=float)

    nA = np.max(A)
    if nA <= 0:
        nA = 1.0
    An = A / nA

    fw_a, prof_a, pk_a = lateral_metrics(An, gx, gz)
    out = dict(
        fwhm_a_mm=fw_a * 1e3,
        peak_a_mm=(pk_a[0] * 1e3, pk_a[1] * 1e3),
        prof_a=prof_a,
    )

    if img_b is not None:
        Bn = np.asarray(img_b, dtype=float) / nA        # normalise BOTH by max(A)
        Dm = np.abs(An - Bn)
        out["mse"] = float(np.mean((An.ravel() - Bn.ravel()) ** 2))
        out["max_abs_err"] = float(np.max(Dm))
        out["diff_db"] = 20.0 * np.log10(Dm + 1e-12)
        fw_b, prof_b, pk_b = lateral_metrics(Bn, gx, gz, force_z=pk_a[1])
        out["fwhm_b_mm"] = fw_b * 1e3
        out["peak_b_mm"] = (pk_b[0] * 1e3, pk_b[1] * 1e3)
        out["prof_b"] = prof_b
    else:
        out["mse"] = np.nan
        out["max_abs_err"] = np.nan
        out["diff_db"] = None
        out["fwhm_b_mm"] = np.nan
        out["peak_b_mm"] = (np.nan, np.nan)
        out["prof_b"] = None
    return out
