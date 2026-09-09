"""RF data generation from a transducer / transmit / scatterer setup.

Python port of ``sim_engine.m``. MATLAB is the source of truth; keep this file
in lock-step with it (see ``CLAUDE.md``).

Coordinate system
-----------------
x : parallel to the array, from the first (leftmost) to the last (rightmost)
    element. x = 0 at the array centre.
z : perpendicular to the array, pointing downward (depth). Elements lie on z = 0.
The steering angle ``theta`` is measured from the z axis, positive towards +x.

Backends
--------
The MATLAB tool wraps SIMUS from MUST (Matlab UltraSound Toolbox) and falls back
to an analytic "mock" backend when MUST is missing. **MUST has no Python
equivalent, so this port implements the mock backend only.** ``backend`` is
always reported as ``'mock'`` and ``force_mock`` is accepted (and ignored) so
the configuration contract mirrors the MATLAB side. The ``MUST das()`` algorithm
option and ``validation/validate_must_vs_mock.m`` remain MATLAB-only.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

__all__ = ["sim_engine", "SimConfig", "SimResult", "TxInfo", "tx_law", "pulse_sigma"]


# --------------------------------------------------------------------------
#  Configuration / result containers
# --------------------------------------------------------------------------
@dataclass
class TxInfo:
    """Transmit settings for one transmit event (the contract every DAS
    implementation relies on). Mirrors the ``tx_info`` struct of the MATLAB
    tool field for field."""

    delays: np.ndarray            # [Nel] transmit delays [s], NaN on inactive elements
    apod: np.ndarray              # [Nel] transmit apodisation, 0 on inactive elements
    elem_x: np.ndarray            # [Nel] transmit element x coordinates [m]
    elem_z: np.ndarray            # [Nel] transmit element z coordinates [m]
    t0: float = 0.0               # time of the first RF sample [s]
    c: float = 1540.0
    fc: float = 5e6
    fs: float = 20e6
    fnumber: float = 1.5          # receive f-number (0 = full aperture)
    rx_apod: str = "rect"         # 'rect' or 'hann'
    scheme: str = "plane"
    angle_deg: float = 0.0
    focus_mm: float = np.nan
    src_xz: tuple = (np.nan, np.nan)   # focal point / virtual source [x z] [m]
    single_element: bool = False


@dataclass
class SimConfig:
    """Flat configuration. Groups mirror the ``cfg.*`` sub-structs of
    ``sim_engine.m`` (``cfg.probe``, ``cfg.tx`` ...)."""

    # probe
    n_elements: int = 64
    pitch: float = 0.30e-3
    fc: float = 5e6
    bandwidth: float = 75.0
    kerf: float | None = None          # default: pitch * 0.1
    height: float = 5e-3
    elevfocus: float = 20e-3
    # medium
    c: float = 1540.0
    # acquisition
    fs_factor: float = 4.0
    # transmit
    scheme: str = "plane"              # 'plane' | 'focused' | 'diverging'
    angle_deg: float = 0.0
    focus_mm: object = 20.0            # scalar or sequence -> multi-focus
    src_mm: float = 10.0              # virtual-source depth for diverging waves
    single_element: bool = False
    element_index: int | None = None   # default: centre element (1-based like MATLAB)
    # scatterers
    scat_x: object = 0.0
    scat_z: object = 20e-3
    scat_rc: object = 1.0
    # reconstruction
    zmax: float = 40e-3
    fnumber: float = 1.5
    rx_apod: str = "rect"
    # options
    force_mock: bool = False           # accepted for API parity; the port is mock-only


@dataclass
class SimResult:
    RF: np.ndarray                     # [Nt x Nel x Ntx]
    t0: float
    fs: float
    c: float
    time: np.ndarray                   # [Nt]
    rx_pos: np.ndarray                 # [Nel] receive element x coordinates [m]
    rx_pos_z: np.ndarray
    tx: list                           # list[TxInfo]
    backend: str
    elapsed: float
    probe: dict
    scat: dict = field(default_factory=dict)


# --------------------------------------------------------------------------
#  Helpers (mirror the local functions of sim_engine.m)
# --------------------------------------------------------------------------
def pulse_sigma(fc: float, bw_percent: float) -> float:
    """Gaussian envelope standard deviation [s] for a fractional bandwidth [%]."""
    bw = max(bw_percent, 1.0) / 100.0
    return np.sqrt(2.0 * np.log(2.0)) / (np.pi * fc * bw)


def _gpulse(t: np.ndarray, fc: float, sigma: float) -> np.ndarray:
    """Gaussian-modulated sinusoidal pulse."""
    return np.exp(-0.5 * (t / sigma) ** 2) * np.cos(2.0 * np.pi * fc * t)


def tx_law(elem_x, c, theta, src_depth, scheme, single_el, el_idx):
    """Transmit delay law. Mirrors ``txLaw`` in ``sim_engine.m``.

    ``el_idx`` is 1-based, matching MATLAB. Returns ``(delays, apod, src_xz)``.
    """
    elem_x = np.asarray(elem_x, dtype=float)
    nel = elem_x.size
    a = np.ones(nel)

    if scheme == "plane" or not np.isfinite(src_depth):
        d = elem_x * np.sin(theta) / c
        src_xz = np.array([np.nan, np.nan])
    else:
        px = src_depth * np.sin(theta)
        pz = src_depth * np.cos(theta)
        r = np.hypot(elem_x - px, pz)
        if src_depth > 0:
            d = (np.max(r) - r) / c            # focused: outer elements fire first
        else:
            d = (r - np.min(r)) / c            # diverging: equidistant from source
        src_xz = np.array([px, pz])

    if single_el:
        el_idx = int(min(max(round(el_idx), 1), nel))
        a[:] = 0.0
        a[el_idx - 1] = 1.0
        d = np.full(nel, np.nan)
        d[el_idx - 1] = 0.0
        src_xz = np.array([elem_x[el_idx - 1], 0.0])

    d = d - np.nanmin(d[np.isfinite(d)])
    d = np.where(a == 0, np.nan, d)
    return d, a, src_xz


# --------------------------------------------------------------------------
#  Mock backend: analytic 2-D cylindrical-wave model
# --------------------------------------------------------------------------
def _run_mock(xs, zs, rc, tx_arr, elem_x, elem_z, c, fc, sigma, tvec):
    nel = elem_x.size
    nt = tvec.size
    ntx = len(tx_arr)
    fs = 1.0 / (tvec[1] - tvec[0])
    RF = np.zeros((nt, nel, ntx))
    tvec = tvec.reshape(-1, 1)                       # [Nt x 1]

    for k, tx in enumerate(tx_arr):
        txd = tx.delays
        apo = tx.apod
        act = np.isfinite(txd) & (apo != 0)
        if not np.any(act):
            continue
        acc = np.zeros((nt, nel))

        exA = elem_x[act]
        ezA = elem_z[act]
        tdA = txd[act]
        apA = apo[act]

        for is_ in range(len(xs)):
            if rc[is_] == 0:
                continue
            # (1) Transmit: field at the scatterer
            rt = np.hypot(xs[is_] - exA, zs[is_] - ezA)             # [Nact]
            taut = tdA + rt / c                                     # [Nact]
            ampt = apA * (zs[is_] / rt) / np.sqrt(rt)               # [Nact]
            # sum over active tx elements (== MATLAB `gpulse(...) * ampt(:)`);
            # np.einsum avoids a spurious NumPy 2.0 matmul SIMD warning.
            sfield = np.einsum("ij,j->i", _gpulse(tvec - taut, fc, sigma), ampt)  # [Nt]

            # (2) Receive: propagate back to every element
            rr = np.hypot(xs[is_] - elem_x, zs[is_] - elem_z)       # [Nel]
            ampr = (zs[is_] / rr) / np.sqrt(rr)                     # [Nel]
            pos = ((tvec - rr / c) - tvec[0]) * fs                  # [Nt x Nel], 0-based
            i0 = np.floor(pos).astype(np.int64)
            frac = pos - i0
            ok = (i0 >= 0) & (i0 < nt - 1)
            i0c = np.clip(i0, 0, nt - 2)
            v = sfield[i0c] * (1.0 - frac) + sfield[i0c + 1] * frac
            v[~ok] = 0.0

            acc += rc[is_] * (v * ampr)
        RF[:, :, k] = acc
    return RF


# --------------------------------------------------------------------------
#  Public entry point
# --------------------------------------------------------------------------
def sim_engine(cfg: SimConfig | None = None) -> SimResult:
    if cfg is None:
        cfg = SimConfig()

    nel = int(cfg.n_elements)
    pitch = float(cfg.pitch)
    fc = float(cfg.fc)
    bandwidth = float(cfg.bandwidth)
    kerf = pitch * 0.1 if cfg.kerf is None else float(cfg.kerf)  # noqa: F841 (parity only)
    c = float(cfg.c)
    fs_factor = float(cfg.fs_factor)
    zmax = float(cfg.zmax)

    if nel < 2 or nel > 1024:
        raise ValueError("n_elements must be in [2, 1024]")
    if not (pitch > 0 and np.isfinite(pitch)):
        raise ValueError("pitch must be positive and finite")
    if not (fc > 0 and np.isfinite(fc)):
        raise ValueError("fc must be positive and finite")
    if fs_factor < 4:
        import warnings
        warnings.warn(
            f"fs_factor = {fs_factor:.2f} is too low; the interpolation-based "
            "DAS will alias. Use 4 or more.",
            stacklevel=2,
        )
    fs = fs_factor * fc

    elem_x = (np.arange(nel) - (nel - 1) / 2.0) * pitch
    elem_z = np.zeros(nel)

    # ---------------- Transmit sequence ----------------
    scheme = str(cfg.scheme).lower()
    angle_deg = float(cfg.angle_deg)
    theta = angle_deg * np.pi / 180.0
    single_el = bool(cfg.single_element)
    el_idx = int(cfg.element_index) if cfg.element_index is not None else round((nel + 1) / 2)

    if scheme == "plane":
        src_depths = np.array([np.nan])
    elif scheme == "focused":
        fmm = np.atleast_1d(np.asarray(cfg.focus_mm, dtype=float)).ravel()
        src_depths = fmm * 1e-3
        src_depths[src_depths <= 0] = 1e-3
    elif scheme == "diverging":
        smm = np.atleast_1d(np.asarray(cfg.src_mm, dtype=float)).ravel()
        src_depths = np.array([-abs(smm[0]) * 1e-3])
    else:
        raise ValueError(f'Unknown transmit scheme "{scheme}".')

    ntx = src_depths.size
    tx_arr: list[TxInfo] = []
    for k in range(ntx):
        d, a, srcxz = tx_law(elem_x, c, theta, src_depths[k], scheme, single_el, el_idx)
        tx_arr.append(
            TxInfo(
                delays=d,
                apod=a,
                elem_x=elem_x.copy(),
                elem_z=elem_z.copy(),
                t0=0.0,
                c=c,
                fc=fc,
                fs=fs,
                fnumber=float(cfg.fnumber),
                rx_apod=str(cfg.rx_apod),
                scheme=scheme,
                angle_deg=angle_deg,
                focus_mm=float(src_depths[k] * 1e3),
                src_xz=(float(srcxz[0]), float(srcxz[1])),
                single_element=single_el,
            )
        )

    # ---------------- Scatterers ----------------
    xs = np.atleast_1d(np.asarray(cfg.scat_x, dtype=float)).ravel()
    zs = np.atleast_1d(np.asarray(cfg.scat_z, dtype=float)).ravel()
    rc = np.atleast_1d(np.asarray(cfg.scat_rc, dtype=float)).ravel()
    if rc.size == 1 and xs.size > 1:
        rc = np.repeat(rc, xs.size)
    if not (xs.size == zs.size == rc.size):
        raise ValueError("Scatterer x, z and rc must have the same number of elements.")
    keep = np.isfinite(xs) & np.isfinite(zs) & np.isfinite(rc) & (zs > 0)
    xs, zs, rc = xs[keep], zs[keep], rc[keep]
    if xs.size == 0:
        xs = np.array([0.0])
        zs = np.array([max(zmax / 2, 1e-3)])
        rc = np.array([0.0])

    # ---------------- Record length ----------------
    ap_half = np.max(np.abs(elem_x)) + pitch
    zrec = max(zmax, np.max(zs)) * 1.05
    sigma = pulse_sigma(fc, bandwidth)
    max_del = max(
        [0.0] + [np.max(t.delays[np.isfinite(t.delays)]) for t in tx_arr]
    )
    t_end = max_del + 2.0 * np.hypot(zrec, 2.0 * ap_half) / c + 8.0 * sigma
    nt = max(64, int(np.ceil(t_end * fs)))
    tvec = np.arange(nt) / fs

    # ---------------- Backend (mock only) ----------------
    t_start = time.perf_counter()
    RF = _run_mock(xs, zs, rc, tx_arr, elem_x, elem_z, c, fc, sigma, tvec)
    backend = "mock (forced)" if cfg.force_mock else "mock (MUST not found)"
    elapsed = time.perf_counter() - t_start

    pk = np.max(np.abs(RF))
    if pk > 0:
        RF = RF / pk

    return SimResult(
        RF=RF,
        t0=0.0,
        fs=fs,
        c=c,
        time=np.arange(RF.shape[0]) / fs,
        rx_pos=elem_x.copy(),
        rx_pos_z=elem_z.copy(),
        tx=tx_arr,
        backend=backend,
        elapsed=elapsed,
        probe=dict(
            n_elements=nel,
            pitch=pitch,
            fc=fc,
            bandwidth=bandwidth,
            kerf=kerf,
            height=cfg.height,
            elevfocus=cfg.elevfocus,
        ),
        scat=dict(x=xs, z=zs, rc=rc),
    )
