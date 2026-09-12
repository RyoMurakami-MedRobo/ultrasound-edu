"""Textbook Delay-And-Sum beamformer (reference implementation).

Python port of ``das_reference.m``. MATLAB is the source of truth.

Unified interface (identical to the MATLAB tool)::

    bmode_img, delays = das_reference(rf_data, tx_info, rx_pos,
                                      grid_x, grid_z, sound_speed, fs,
                                      want_delays=False)

rf_data     [Nt x Nel] RF data, time along the first dimension.
tx_info     a :class:`~ultrasound_das.sim_engine.TxInfo` (or any object /
            mapping exposing the same fields).
rx_pos      [Nel] (x, z = 0) or [2 x Nel] receive element coordinates [m].
grid_x,grid_z  reconstruction axes [m]; expanded with meshgrid. Image size is
            ``(len(grid_z), len(grid_x))``.
bmode_img   [Nz x Nx] **linear envelope** (no log compression).
delays      [Npix x Nel] total two-way time [s] per pixel and element, NaN
            where not summed. Pixels are column-major (MATLAB ``reshape``
            order), i.e. ``bmode_img.ravel(order='F')``. Only returned when
            ``want_delays=True`` (mirrors ``nargout < 2``).
"""

# want_trace=True returns (image, delays, trace) with actual channel increments.
# Keep the accumulation hook in custom implementations; see docs/HW3_DAS.md.

from __future__ import annotations

import numpy as np

__all__ = ["das_reference", "tx_arrival_time", "envelope_z"]


def _get(tx, name, default):
    if isinstance(tx, dict):
        v = tx.get(name, default)
    else:
        v = getattr(tx, name, default)
    if v is None:
        return default
    return v


def _normalize_rx_pos(rx_pos, nel):
    rx_pos = np.asarray(rx_pos, dtype=float)
    if rx_pos.ndim == 1:
        rxx = rx_pos.ravel()
        rxz = np.zeros(rxx.size)
    elif rx_pos.shape[0] == 2:
        rxx, rxz = rx_pos[0, :], rx_pos[1, :]
    elif rx_pos.shape[1] == 2:
        rxx, rxz = rx_pos[:, 0], rx_pos[:, 1]
    else:
        raise ValueError("rx_pos must be [Nel] or [2 x Nel].")
    if rxx.size != nel:
        raise ValueError(
            f"rx_pos has {rxx.size} elements but the RF data has {nel} columns."
        )
    return rxx, rxz


def _normalize_grid(grid_x, grid_z):
    grid_x = np.asarray(grid_x, dtype=float)
    grid_z = np.asarray(grid_z, dtype=float)
    if grid_x.ndim <= 1 and grid_z.ndim <= 1:
        gx = np.atleast_1d(grid_x).ravel()
        gz = np.atleast_1d(grid_z).ravel()
        XI, ZI = np.meshgrid(gx, gz)            # (Nz, Nx)
    elif grid_x.shape == grid_z.shape:
        XI, ZI = grid_x, grid_z
    else:
        raise ValueError("grid_x and grid_z have incompatible sizes.")
    return XI, ZI, XI.shape


def envelope_z(x: np.ndarray) -> np.ndarray:
    """Envelope detection along the first (depth) dimension.

    Uses an FFT-built analytic signal (equivalent to ``scipy.signal.hilbert``
    with ``axis=0``), matching the toolbox-free branch of ``das_reference.m``.
    """
    x = np.asarray(x, dtype=float)
    if x.shape[0] < 4:
        return np.abs(x)
    n = x.shape[0]
    X = np.fft.fft(x, n=n, axis=0)
    h = np.zeros(n)
    h[0] = 1.0
    if n % 2 == 0:
        h[n // 2] = 1.0
        h[1 : n // 2] = 2.0
    else:
        h[1 : (n + 1) // 2] = 2.0
    h = h.reshape([-1] + [1] * (x.ndim - 1))
    return np.abs(np.fft.ifft(X * h, n=n, axis=0))


def tx_arrival_time(XI, ZI, tx_info, c) -> np.ndarray:
    """Time [s] at which the transmit wavefront reaches each pixel, [Npix]
    (column-major flatten of XI / ZI)."""
    ex = np.asarray(_get(tx_info, "elem_x", None), dtype=float)
    td = np.asarray(_get(tx_info, "delays", None), dtype=float)
    ap = _get(tx_info, "apod", None)
    ez = _get(tx_info, "elem_z", None)

    if ex is None or td is None or ex.size == 0 or td.size == 0:
        raise ValueError("tx_info must provide elem_x and delays.")
    ez = np.zeros_like(ex) if ez is None else np.asarray(ez, dtype=float)
    ap = np.ones_like(ex) if ap is None else np.asarray(ap, dtype=float)

    act = np.isfinite(td) & (ap.ravel() != 0)
    if not np.any(act):
        raise ValueError("No active transmit element.")
    ex, ez, td = ex[act], ez[act], td[act]

    xg = XI.ravel(order="F")
    zg = ZI.ravel(order="F")

    src = np.asarray(_get(tx_info, "src_xz", (np.nan, np.nan)), dtype=float)
    scheme = str(_get(tx_info, "scheme", ""))
    single = bool(_get(tx_info, "single_element", False))
    is_focused = (
        scheme.lower() == "focused"
        and src.size == 2
        and np.all(np.isfinite(src))
        and src[1] > 0
        and not single
    )

    if is_focused:
        Rf = np.hypot(ex - src[0], ez - src[1])
        Tf = np.mean(td + Rf / c)
        D = np.hypot(src[0], src[1])
        u = src / D
        proj = xg * u[0] + zg * u[1]
        dF = np.hypot(xg - src[0], zg - src[1])
        sgn = np.ones(xg.size)
        sgn[proj <= D] = -1.0
        tau = Tf + sgn * dF / c
    else:
        tau = np.full(xg.size, np.inf)
        for e in range(ex.size):
            tau = np.minimum(tau, td[e] + np.hypot(xg - ex[e], zg - ez[e]) / c)
    return tau


def das_reference(
    rf_data,
    tx_info,
    rx_pos,
    grid_x,
    grid_z,
    sound_speed,
    fs,
    want_delays: bool = False, want_trace: bool = False,
):
    rf_data = np.asarray(rf_data, dtype=float)
    if rf_data.ndim > 2:
        rf_data = rf_data[:, :, 0]
    nt, nel = rf_data.shape

    rxx, rxz = _normalize_rx_pos(rx_pos, nel)
    XI, ZI, img_size = _normalize_grid(grid_x, grid_z)
    npix = XI.size

    c = float(sound_speed)
    t0 = float(_get(tx_info, "t0", 0.0))
    fnum = float(_get(tx_info, "fnumber", 0.0))
    rx_apod = str(_get(tx_info, "rx_apod", "rect")).lower()
    if not np.isfinite(fnum) or fnum <= 0:
        fnum = 0.0

    tau_tx = tx_arrival_time(XI, ZI, tx_info, c)

    bf = np.zeros(npix)
    delays = np.full((npix, nel), np.nan) if (want_delays or want_trace) else None

    xg = XI.ravel(order="F")
    zg = ZI.ravel(order="F")
    if fnum > 0:
        half_ap = zg / (2.0 * fnum)
    else:
        half_ap = np.full(npix, np.inf)

    contributions = np.zeros((npix, nel)) if want_trace else None

    for e in range(nel):
        dx = xg - rxx[e]
        rrx = np.hypot(dx, zg - rxz[e])
        tau = tau_tx + rrx / c

        in_ap = np.abs(dx) <= half_ap
        if rx_apod == "hann":
            w = np.zeros(npix)
            denom = np.where(in_ap, np.maximum(half_ap, np.finfo(float).eps), 1.0)
            u = np.where(in_ap, dx / denom, 0.0)
            w = in_ap * (0.5 * (1.0 + np.cos(np.pi * u)))
        else:
            w = in_ap.astype(float)

        pos = (tau - t0) * fs                     # 0-based sample position
        i0 = np.floor(pos).astype(np.int64)
        frac = pos - i0
        ok = (w > 0) & (i0 >= 0) & (i0 < nt - 1)
        i0c = np.clip(i0, 0, nt - 2)
        val = rf_data[i0c, e] * (1.0 - frac) + rf_data[i0c + 1, e] * frac
        val = np.where(ok, val, 0.0)

        previous = bf.copy() if want_trace else None
        bf += w * val
        if want_trace:
            contributions[:, e] = bf - previous

        if want_delays or want_trace:
            d = tau.copy()
            d[~ok] = np.nan
            delays[:, e] = d

    bmode_img = envelope_z(bf.reshape(img_size, order="F"))
    if want_trace:
        return bmode_img, delays, {"contributions": contributions,
                                   "coherent": bf.reshape(img_size, order="F")}
    if want_delays:
        return bmode_img, delays
    return bmode_img
