"""Algorithm dispatch and multi-transmit compositing.

Python port of the ``runAlgorithm`` / ``runAllTx`` / ``algoHandle`` logic in
``main_gui.m``. MATLAB is the source of truth.
"""

from __future__ import annotations

import time

import numpy as np

from .das_custom_template import das_custom_template
from .das_reference import das_reference

__all__ = ["run_all_tx", "run_algorithm", "recon_grid", "ALGORITHMS"]

# "MUST das()" from the MATLAB tool has no Python equivalent (MUST is
# MATLAB-only). Only the two pure-Python algorithms are exposed here.
ALGORITHMS = {
    "Reference DAS": das_reference,
    "Custom DAS": das_custom_template,
}


def recon_grid(x_half_mm=12.0, nx=161, z_min_mm=3.0, z_max_mm=40.0, nz=221):
    """Reconstruction axes [m], mirroring ``reconGrid`` in ``main_gui.m``."""
    xh = x_half_mm * 1e-3
    gx = np.linspace(-xh, xh, int(round(nx)))
    z0 = min(z_min_mm, z_max_mm - 1) * 1e-3
    gz = np.linspace(max(z0, 1e-4), z_max_mm * 1e-3, int(round(nz)))
    return gx, gz


def run_all_tx(fh, result, gx, gz):
    """Run beamformer ``fh`` over every transmit event of ``result``.

    A single event is beamformed directly. A multi-focus sequence is composited
    zone by zone in depth (each zone uses the event whose focal depth is
    nearest), exactly like ``runAllTx``.
    """
    S = result
    ntx = len(S.tx)
    if ntx == 1:
        return fh(S.RF[:, :, 0], S.tx[0], S.rx_pos, gx, gz, S.c, S.fs)

    fz = np.array([t.focus_mm for t in S.tx]) * 1e-3
    order = np.argsort(fz)
    fzs = fz[order]
    edges = np.concatenate(([-np.inf], (fzs[:-1] + fzs[1:]) / 2.0, [np.inf]))
    img = np.zeros((gz.size, gx.size))
    for k in range(ntx):
        sub = (gz >= edges[k]) & (gz < edges[k + 1])
        if not np.any(sub):
            continue
        j = int(order[k])
        img[sub, :] = fh(S.RF[:, :, j], S.tx[j], S.rx_pos, gx, gz[sub], S.c, S.fs)
    return img


def run_algorithm(key_or_callable, result, gx, gz):
    """Beamform with a named algorithm or a custom callable.

    Returns ``(img, elapsed_ms, name)``. Mirrors ``runAlgorithm``: a throwaway
    warm-up call on a tiny grid, then the timed call requesting a single output.
    """
    if callable(key_or_callable):
        fh = key_or_callable
        name = getattr(fh, "__name__", "custom")
    else:
        if key_or_callable not in ALGORITHMS:
            raise ValueError(
                f'Unknown algorithm "{key_or_callable}". '
                f"Choices: {sorted(ALGORITHMS)} (MUST das() is MATLAB-only)."
            )
        fh = ALGORITHMS[key_or_callable]
        name = fh.__name__

    S = result
    try:
        fh(S.RF[:, :, 0], S.tx[0], S.rx_pos, gx[: min(4, gx.size)],
           gz[: min(4, gz.size)], S.c, S.fs)
    except Exception:
        pass

    t0 = time.perf_counter()
    img = run_all_tx(fh, S, gx, gz)
    ms = (time.perf_counter() - t0) * 1000.0

    expected = (gz.size, gx.size)
    if img.shape != expected:
        raise ValueError(
            f"{name} returned an image of size {img.shape} but {expected} was "
            "expected. The unified interface requires a linear envelope of size "
            "(len(grid_z), len(grid_x))."
        )
    if np.iscomplexobj(img) or not np.all(np.isfinite(img)):
        img = np.abs(img)
        img[~np.isfinite(img)] = 0.0
    return img, ms, name
