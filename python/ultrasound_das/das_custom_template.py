"""Skeleton for writing your own DAS beamformer (Python port of
``das_custom_template.m``).

HOW TO USE
    1) Copy this file, e.g. to ``my_das.py``.
    2) Rename ``das_custom_template`` to ``my_das``.
    3) In the GUI (Detailed mode) type ``my_das`` into "Custom function name"
       and pick "Custom DAS" as algorithm A or B, or pass the callable to
       :func:`ultrasound_das.pipeline.run_algorithm`.

CONTRACT (identical to :func:`ultrasound_das.das_reference.das_reference`)
    - Keep the argument order and meaning.
    - ``rf_data`` is [Nt x Nel] (time along the first dimension).
    - ``bmode_img`` must be [len(grid_z) x len(grid_x)] and the **linear
      envelope** (no log compression).
    - ``delays`` must be [Npix x Nel] holding the total two-way time [s], NaN
      where a sample is not summed. Pixels are column-major
      (``img.ravel(order='F')``). Only build it when ``want_delays=True``.

INITIAL STATE
    A deliberately naive DAS (nearest-neighbour sampling + rectangular
    apodisation) so the file runs out of the box. Comparing it with the
    reference (linear interpolation) shows the interpolation error as a raised
    sidelobe floor in the difference image and the lateral profile.
"""

# want_trace=True returns (image, delays, trace) with actual channel increments.
# Keep the accumulation hook in custom implementations; see docs/HW3_DAS.md.

from __future__ import annotations

import numpy as np


def _get(tx, name, default):
    v = tx.get(name, default) if isinstance(tx, dict) else getattr(tx, name, default)
    return default if v is None else v


def das_custom_template(
    rf_data, tx_info, rx_pos, grid_x, grid_z, sound_speed, fs, want_delays: bool = False, want_trace: bool = False
):
    # ---- [STEP 0] Input handling (usually left untouched) ----
    rf_data = np.asarray(rf_data, dtype=float)
    if rf_data.ndim > 2:
        rf_data = rf_data[:, :, 0]
    nt, nel = rf_data.shape

    rxx = np.asarray(rx_pos, dtype=float).ravel()
    rxz = np.zeros(nel)
    if rxx.size != nel:
        raise ValueError("rx_pos does not match the number of RF columns.")

    gx = np.atleast_1d(np.asarray(grid_x, dtype=float))
    gz = np.atleast_1d(np.asarray(grid_z, dtype=float))
    if gx.ndim == 1 and gz.ndim == 1:
        XI, ZI = np.meshgrid(gx, gz)
    else:
        XI, ZI = gx, gz
    img_size = XI.shape
    xg = XI.ravel(order="F")
    zg = ZI.ravel(order="F")
    npix = xg.size

    c = float(sound_speed)
    t0 = float(_get(tx_info, "t0", 0.0))

    # ---- [STEP 1] Transmit arrival time tau_tx(p)   <<< EDIT HERE ----
    ex = np.asarray(_get(tx_info, "elem_x", None), dtype=float)
    ez = np.asarray(_get(tx_info, "elem_z", np.zeros_like(ex)), dtype=float)
    td = np.asarray(_get(tx_info, "delays", None), dtype=float)
    ap = np.asarray(_get(tx_info, "apod", np.ones_like(ex)), dtype=float)
    act = np.isfinite(td) & (ap.ravel() != 0)
    exA, ezA, tdA = ex[act], ez[act], td[act]

    src = np.asarray(_get(tx_info, "src_xz", (np.nan, np.nan)), dtype=float)
    use_virtual_source = (
        str(_get(tx_info, "scheme", "")).lower() == "focused"
        and np.all(np.isfinite(src))
        and src[1] > 0
        and not bool(_get(tx_info, "single_element", False))
    )

    if use_virtual_source:
        Rf = np.hypot(exA - src[0], ezA - src[1])
        Tf = np.mean(tdA + Rf / c)
        D = np.hypot(src[0], src[1])
        u = src / D
        proj = xg * u[0] + zg * u[1]
        dF = np.hypot(xg - src[0], zg - src[1])
        sgn = np.ones(npix)
        sgn[proj <= D] = -1.0
        tau_tx = Tf + sgn * dF / c
    else:
        tau_tx = np.full(npix, np.inf)
        for e in range(exA.size):
            tau_tx = np.minimum(tau_tx, tdA[e] + np.hypot(xg - exA[e], zg - ezA[e]) / c)

    # ---- [STEP 2] Receive aperture weights (f-number, apodisation) ----
    fnum = float(_get(tx_info, "fnumber", 0.0))
    if not np.isfinite(fnum) or fnum <= 0:
        half_ap = np.full(npix, np.inf)
    else:
        half_ap = zg / (2.0 * fnum)

    # ---- [STEP 3] Delay-and-sum main loop   <<< EDIT HERE ----
    #  Nearest-neighbour sampling. Replacing round() with floor() + linear
    #  interpolation should reproduce the reference implementation exactly.
    bf = np.zeros(npix)
    delays = np.full((npix, nel), np.nan) if (want_delays or want_trace) else None

    contributions = np.zeros((npix, nel)) if want_trace else None

    for e in range(nel):
        dx = xg - rxx[e]
        tau = tau_tx + np.hypot(dx, zg - rxz[e]) / c
        w = (np.abs(dx) <= half_ap).astype(float)

        # nearest-neighbour, 0-based. floor(x + 0.5) matches MATLAB round()
        # (round-half-away-from-zero) for the non-negative sample positions here.
        idx = np.floor((tau - t0) * fs + 0.5).astype(np.int64)
        ok = (w > 0) & (idx >= 0) & (idx < nt)
        idxc = np.clip(idx, 0, nt - 1)
        val = np.where(ok, rf_data[idxc, e], 0.0)

        previous = bf.copy() if want_trace else None
        bf += w * val
        if want_trace:
            contributions[:, e] = bf - previous
        if want_delays or want_trace:
            d = tau.copy()
            d[~ok] = np.nan
            delays[:, e] = d

    # ---- [STEP 4] Envelope detection   <<< EDITABLE ----
    bf2 = bf.reshape(img_size, order="F")
    if img_size[0] >= 4:
        from .das_reference import envelope_z

        bmode_img = envelope_z(bf2)
    else:
        bmode_img = np.abs(bf2)

    # ---- [STEP 5] Post-processing (TGC, speckle reduction, ...) ----

    if want_trace:
        return bmode_img, delays, {"contributions": contributions,
                                   "coherent": bf.reshape(img_size, order="F")}
    if want_delays:
        return bmode_img, delays
    return bmode_img
