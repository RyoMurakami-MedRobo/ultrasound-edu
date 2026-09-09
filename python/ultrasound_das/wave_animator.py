"""Wave propagation, delay-curve and before/after alignment computations.

Python port of ``wave_animator.m``. MATLAB is the source of truth.

Deliberate MATLAB<->Python divergences (keep both sides in sync with this in
mind, see ``CLAUDE.md``):

1. ``wave_animator.m`` mixes computation and MATLAB drawing in one file. Here
   the numeric part lives in this module (returning plain arrays) and the
   drawing lives in ``gui.py``. Port numeric changes here; port drawing
   changes to ``gui.py``.
2. ``wave_animator.m`` deliberately *duplicates* its own ``txArrivalTime``
   local function (see the comment at ``wave_animator.m:311``) to stay
   independent of the beamformer. This port instead reuses
   :func:`ultrasound_das.das_reference.tx_arrival_time`. The wavefront model
   is therefore guaranteed identical to the beamformer's by construction; if
   the MATLAB duplicate is ever changed without changing ``das_reference.m``,
   that divergence must be reproduced here on purpose.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .das_reference import tx_arrival_time

__all__ = ["WavePropagation", "delay_curve_image", "alignment_bundles", "blue_white_red"]

_MAXSCAT = 30


def blue_white_red(n: int = 128) -> np.ndarray:
    """Blue-white-red diverging colormap as an (2n, 3) RGB array."""
    u = np.linspace(0.0, 1.0, n).reshape(-1, 1)
    lower = np.hstack([u, u, np.ones((n, 1))])
    upper = np.hstack([np.ones((n, 1)), np.flipud(u), np.flipud(u)])
    return np.vstack([lower, upper])


def _sample_rf(col, t0, fs, times):
    """Linear interpolation of one RF channel at arbitrary times (zero outside
    the record). Mirrors ``sampleRF``."""
    col = np.asarray(col, dtype=float)
    times = np.asarray(times, dtype=float).ravel()
    n = col.size
    pos = (times - t0) * fs                     # 0-based
    i0 = np.floor(pos).astype(np.int64)
    fr = pos - i0
    ok = (i0 >= 0) & (i0 < n - 1)
    i0c = np.clip(i0, 0, n - 2)
    v = col[i0c] * (1.0 - fr) + col[i0c + 1] * fr
    v[~ok] = 0.0
    return v


@dataclass
class WavePropagation:
    """Precomputed arrival-time maps for the wavefront animation, plus a
    :meth:`frame` method returning the field snapshot at time ``t``.

    Mirrors ``setupPropagation`` / ``drawPropagation`` (numeric parts only).
    """

    xv: np.ndarray
    zv: np.ndarray
    X: np.ndarray
    Z: np.ndarray
    TAU: np.ndarray
    c: float
    sigma: float
    xs: np.ndarray
    zs: np.ndarray
    rc: np.ndarray
    tau_s: np.ndarray
    Rmat: np.ndarray
    n_dropped: int
    elem_x: np.ndarray
    elem_z: np.ndarray
    act_elem_x: np.ndarray
    src_xz: tuple
    scheme: str
    tmax: float

    @classmethod
    def setup(cls, result, tx_idx, xrange, zrange, nx=241, nz=241):
        S = result
        tx = S.tx[tx_idx]
        c = S.c

        xv = np.linspace(xrange[0], xrange[1], nx)
        zv = np.linspace(max(zrange[0], 0.0), zrange[1], nz)
        X, Z = np.meshgrid(xv, zv)

        TAU = tx_arrival_time(X, Z, tx, c).reshape(X.shape, order="F")

        xs = np.asarray(S.scat["x"], dtype=float)
        zs = np.asarray(S.scat["z"], dtype=float)
        rc = np.asarray(S.scat["rc"], dtype=float)
        keep = rc != 0
        xs, zs, rc = xs[keep], zs[keep], rc[keep]
        n_dropped = 0
        if xs.size > _MAXSCAT:
            order = np.argsort(np.abs(rc))[::-1][:_MAXSCAT]
            n_dropped = xs.size - _MAXSCAT
            xs, zs, rc = xs[order], zs[order], rc[order]

        ns = xs.size
        tau_s = np.zeros(ns)
        Rmat = np.zeros((nz, nx, ns))
        for k in range(ns):
            tau_s[k] = float(tx_arrival_time(np.array([[xs[k]]]), np.array([[zs[k]]]), tx, c)[0])
            Rmat[:, :, k] = np.hypot(X - xs[k], Z - zs[k])

        fc = S.probe["fc"]
        bw = max(S.probe["bandwidth"], 1) / 100.0
        sigma = np.sqrt(2.0 * np.log(2.0)) / (np.pi * fc * bw)
        sigma = max(sigma, 0.3 / fc)

        elem_x = np.asarray(tx.elem_x, dtype=float)
        elem_z = np.asarray(tx.elem_z, dtype=float)
        act = np.isfinite(tx.delays) & (np.asarray(tx.apod) != 0)

        if ns > 0:
            dmax = np.max(np.hypot(xs[:, None] - elem_x[None, :], zs[:, None] - elem_z[None, :]), axis=1)
            tmax = float(np.max(tau_s + dmax / c) + 6.0 * sigma)
        else:
            tmax = float(np.max(TAU) + 6.0 * sigma)
        tmax = max(tmax, 6.0 * sigma)

        return cls(
            xv=xv, zv=zv, X=X, Z=Z, TAU=TAU, c=c, sigma=sigma,
            xs=xs, zs=zs, rc=rc, tau_s=tau_s, Rmat=Rmat, n_dropped=n_dropped,
            elem_x=elem_x, elem_z=elem_z, act_elem_x=elem_x[act],
            src_xz=tuple(tx.src_xz), scheme=tx.scheme, tmax=tmax,
        )

    def frame(self, t: float):
        """Field snapshot at time ``t`` [s]. Returns ``(field, hit_x)`` where
        ``field`` is in [-1, 1] (red = transmit wavefront, blue = echo) and
        ``hit_x`` are the x coordinates of receive elements the echo just
        reached."""
        sig = self.sigma
        F = np.exp(-0.5 * ((self.TAU - t) / sig) ** 2)

        Fe = np.zeros_like(F)
        hit_x = []
        for k in range(self.xs.size):
            dt = t - self.tau_s[k]
            if dt <= 0:
                continue
            amp = min(abs(self.rc[k]), 1.0) * 0.9
            Fe = np.maximum(Fe, amp * np.exp(-0.5 * ((self.Rmat[:, :, k] / self.c - dt) / sig) ** 2))
            re = np.hypot(self.elem_x - self.xs[k], self.elem_z - self.zs[k])
            hit_x.extend(self.elem_x[np.abs(re / self.c - dt) < 2.0 * sig].tolist())

        field = np.clip(F - Fe, -1.0, 1.0)
        return field, np.asarray(hit_x)


def delay_curve_image(result, tx_idx, tau):
    """Display data for the RF-image + delay-curve overlay.

    Returns ``(disp_rf, t_us, tau_us, ok)`` where ``disp_rf`` is the
    depth-gain-compensated RF image (for display only), ``t_us`` the time axis
    in microseconds, ``tau_us`` the per-element summed-sample time, and ``ok``
    the finite mask. Mirrors ``drawDelayCurve`` (data only)."""
    S = result
    RF = S.RF[:, :, tx_idx]
    t_us = (S.t0 + np.arange(RF.shape[0]) / S.fs) * 1e6
    denom = np.maximum(np.max(np.abs(RF), axis=1, keepdims=True), 1e-12) ** 0.7
    disp_rf = RF / denom
    tau = np.asarray(tau, dtype=float).ravel()
    ok = np.isfinite(tau)
    return disp_rf, t_us, tau * 1e6, ok


def alignment_bundles(result, tx_idx, tau, n_periods=5, n_samples=401, gain=3.0):
    """Before/after alignment waveform bundles.

    Returns a dict with ``trel_us``, ``pre`` [Nsamp x Nel], ``post``
    [Nsamp x Nel], ``ok`` mask, ``scale``, ``gain``, ``tc`` (window centre
    time), ``sum_trace`` (coherent sum of aligned active channels).
    Mirrors ``drawAlignment`` (data only)."""
    S = result
    RF = S.RF[:, :, tx_idx]
    nel = RF.shape[1]
    fc = S.probe["fc"]

    tau = np.asarray(tau, dtype=float).ravel()
    ok = np.isfinite(tau)
    if not np.any(ok):
        return None

    trel = np.linspace(-n_periods / fc, n_periods / fc, n_samples)
    tc = float(np.mean(tau[ok]))

    pre = np.zeros((trel.size, nel))
    post = np.zeros((trel.size, nel))
    for e in range(nel):
        pre[:, e] = _sample_rf(RF[:, e], S.t0, S.fs, tc + trel)
        if ok[e]:
            post[:, e] = _sample_rf(RF[:, e], S.t0, S.fs, tau[e] + trel)

    scale = max(np.max(np.abs(pre)), np.max(np.abs(post)), np.finfo(float).eps)
    sum_trace = np.sum(post[:, ok], axis=1)
    return dict(
        trel_us=trel * 1e6, pre=pre, post=post, ok=ok,
        scale=float(scale), gain=gain, tc=tc, sum_trace=sum_trace, tau=tau,
    )
