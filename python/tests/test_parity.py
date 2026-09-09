"""Parity tests: the Python port must reproduce the MATLAB reference to
floating-point precision.

Golden fixtures under ``python/parity/fixtures/*.mat`` are produced by
``python/parity/dump_reference.m`` running the real MATLAB
``sim_engine.m`` / ``das_reference.m`` / ``das_custom_template.m``. Regenerate
them whenever the MATLAB numerics change (see ``CLAUDE.md``).

If a test only passes at a loose tolerance, that is a structural port bug
(column-major vs. C order, hilbert axis, 1-based vs. 0-based sample index,
FWHM off-by-one), not float noise -- do not loosen the tolerance.

Run with ``pytest`` or directly: ``python python/tests/test_parity.py``.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import scipy.io as sio

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from ultrasound_das import das_custom_template, das_reference  # noqa: E402
from ultrasound_das import SimConfig, compare_images, run_all_tx, sim_engine  # noqa: E402

FIXTURE_DIR = pathlib.Path(__file__).resolve().parents[1] / "parity" / "fixtures"
RTOL = 1e-9
ATOL = 1e-12


def _fixtures():
    return sorted(FIXTURE_DIR.glob("*.mat"))


def _scalar(v):
    return np.asarray(v).ravel()[0]


def _cfg_from_fixture(mat) -> SimConfig:
    c = mat["c"]
    g = lambda name: c[name][0, 0]  # noqa: E731

    def _vec(name):
        return np.asarray(g(name), dtype=float).ravel()

    return SimConfig(
        n_elements=int(_scalar(g("n_elements"))),
        pitch=float(_scalar(g("pitch"))),
        fc=float(_scalar(g("fc"))),
        bandwidth=float(_scalar(g("bandwidth"))),
        c=float(_scalar(g("c"))),
        fs_factor=float(_scalar(g("fs_factor"))),
        scheme=str(g("scheme")[0]) if np.asarray(g("scheme")).size else "plane",
        angle_deg=float(_scalar(g("angle_deg"))),
        focus_mm=_vec("focus_mm"),
        src_mm=float(_scalar(g("src_mm"))),
        single_element=bool(_scalar(g("single_element"))),
        element_index=int(_scalar(g("element_index"))),
        scat_x=_vec("scat_x"),
        scat_z=_vec("scat_z"),
        scat_rc=_vec("scat_rc"),
        zmax=float(_scalar(g("zmax"))),
        fnumber=float(_scalar(g("fnumber"))),
        rx_apod=str(g("rx_apod")[0]),
        force_mock=True,
    )


def _run_python(mat):
    cfg = _cfg_from_fixture(mat)
    S = sim_engine(cfg)
    for tx in S.tx:
        tx.fnumber = cfg.fnumber
        tx.rx_apod = cfg.rx_apod
    gx = np.asarray(mat["gx"], dtype=float).ravel()
    gz = np.asarray(mat["gz"], dtype=float).ravel()
    return cfg, S, gx, gz


def _check_one(path):
    mat = sio.loadmat(str(path))
    label = path.stem
    cfg, S, gx, gz = _run_python(mat)

    # ---- RF ----
    RF_ref = np.asarray(mat["RF"], dtype=float)
    ntx = int(_scalar(mat["ntx"]))
    if RF_ref.ndim == 2:
        RF_ref = RF_ref[:, :, None]
    assert S.RF.shape == RF_ref.shape, f"{label}: RF shape {S.RF.shape} != {RF_ref.shape}"
    np.testing.assert_allclose(S.RF, RF_ref, rtol=RTOL, atol=ATOL, err_msg=f"{label}: RF")
    np.testing.assert_allclose(S.fs, float(_scalar(mat["fs"])), rtol=RTOL, err_msg=f"{label}: fs")
    np.testing.assert_allclose(S.t0, float(_scalar(mat["t0"])), atol=ATOL, err_msg=f"{label}: t0")

    # ---- transmit delay law ----
    tx_delays = np.asarray(mat["tx_delays"], dtype=float)   # [Nel x Ntx]
    tx_apod = np.asarray(mat["tx_apod"], dtype=float)
    for k in range(ntx):
        np.testing.assert_allclose(
            S.tx[k].delays, tx_delays[:, k], rtol=RTOL, atol=ATOL,
            err_msg=f"{label}: tx[{k}].delays", equal_nan=True,
        )
        np.testing.assert_allclose(
            S.tx[k].apod, tx_apod[:, k], rtol=RTOL, atol=ATOL, err_msg=f"{label}: tx[{k}].apod",
        )

    # ---- beamformed images ----
    has_delays = bool(_scalar(mat["has_delays"]))
    if ntx == 1:
        img_ref_py, delays_ref_py = das_reference(
            S.RF[:, :, 0], S.tx[0], S.rx_pos, gx, gz, S.c, S.fs, want_delays=True
        )
        img_cus_py, delays_cus_py = das_custom_template(
            S.RF[:, :, 0], S.tx[0], S.rx_pos, gx, gz, S.c, S.fs, want_delays=True
        )
    else:
        img_ref_py = run_all_tx(das_reference, S, gx, gz)
        img_cus_py = run_all_tx(das_custom_template, S, gx, gz)
        delays_ref_py = delays_cus_py = None

    np.testing.assert_allclose(
        img_ref_py, np.asarray(mat["img_ref"], dtype=float),
        rtol=1e-8, atol=1e-11, err_msg=f"{label}: img_ref",
    )
    np.testing.assert_allclose(
        img_cus_py, np.asarray(mat["img_custom"], dtype=float),
        rtol=1e-8, atol=1e-11, err_msg=f"{label}: img_custom",
    )

    if has_delays:
        np.testing.assert_allclose(
            delays_ref_py, np.asarray(mat["delays_ref"], dtype=float),
            rtol=RTOL, atol=ATOL, equal_nan=True, err_msg=f"{label}: delays_ref",
        )
        np.testing.assert_allclose(
            delays_cus_py, np.asarray(mat["delays_custom"], dtype=float),
            rtol=RTOL, atol=ATOL, equal_nan=True, err_msg=f"{label}: delays_custom",
        )

    # ---- comparison metrics ----
    met = mat["met"]
    cmp = compare_images(img_ref_py, img_cus_py, gx, gz)
    for key in ("fwhm_a_mm", "fwhm_b_mm", "mse", "max_abs_err"):
        exp = float(_scalar(met[key][0, 0]))
        got = float(cmp[key])
        if np.isnan(exp):
            assert np.isnan(got), f"{label}: {key} expected NaN, got {got}"
        else:
            np.testing.assert_allclose(got, exp, rtol=1e-7, atol=1e-10, err_msg=f"{label}: {key}")
    for key, pk in (("peak_a_mm", cmp["peak_a_mm"]), ("peak_b_mm", cmp["peak_b_mm"])):
        exp = np.asarray(met[key][0, 0], dtype=float).ravel()
        np.testing.assert_allclose(np.asarray(pk), exp, rtol=1e-9, atol=1e-9, err_msg=f"{label}: {key}")

    return label


def test_parity_all_fixtures():
    paths = _fixtures()
    assert paths, f"no fixtures in {FIXTURE_DIR}; run parity/dump_reference.m"
    for path in paths:
        _check_one(path)


def test_sim_engine_amplitude_normalised():
    S = sim_engine(SimConfig())
    assert np.isclose(np.max(np.abs(S.RF)), 1.0)


if __name__ == "__main__":
    ok = 0
    for path in _fixtures():
        label = _check_one(path)
        print(f"  PASS  {label}")
        ok += 1
    test_sim_engine_amplitude_normalised()
    print(f"\n{ok} fixtures OK")
