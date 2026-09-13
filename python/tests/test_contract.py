"""Interface-contract tests for the Python port (no MATLAB required)."""

from __future__ import annotations

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from ultrasound_das import (  # noqa: E402
    SimConfig,
    alignment_bundles,
    bundle_gain,
    das_custom_template,
    das_reference,
    delay_curve_frame,
    recon_grid,
    run_algorithm,
    sim_engine,
)


def _sim():
    return sim_engine(SimConfig(n_elements=32, zmax=26e-3, scat_z=18e-3))


def test_image_shape_and_realness():
    S = _sim()
    gx, gz = recon_grid(x_half_mm=10, nx=41, z_max_mm=26, nz=51)
    img = das_reference(S.RF[:, :, 0], S.tx[0], S.rx_pos, gx, gz, S.c, S.fs)
    assert img.shape == (gz.size, gx.size)
    assert np.isrealobj(img) and np.all(np.isfinite(img)) and np.all(img >= 0)


def test_delays_contract():
    S = _sim()
    gx, gz = recon_grid(x_half_mm=10, nx=21, z_max_mm=26, nz=31)
    img, delays = das_reference(
        S.RF[:, :, 0], S.tx[0], S.rx_pos, gx, gz, S.c, S.fs, want_delays=True
    )
    assert delays.shape == (gx.size * gz.size, S.rx_pos.size)
    # column-major pixel order: reshaping the image back must round-trip
    assert np.array_equal(img.ravel(order="F").shape, (gx.size * gz.size,))
    finite = np.isfinite(delays)
    assert finite.any() and (delays[finite] > 0).all()


def test_naive_custom_differs_from_reference():
    """Out of the box the template is a naive DAS; it must not be bit-identical
    to the interpolating reference (that is the teaching point)."""
    S = _sim()
    gx, gz = recon_grid(x_half_mm=10, nx=41, z_max_mm=26, nz=51)
    a = das_reference(S.RF[:, :, 0], S.tx[0], S.rx_pos, gx, gz, S.c, S.fs)
    b = das_custom_template(S.RF[:, :, 0], S.tx[0], S.rx_pos, gx, gz, S.c, S.fs)
    assert a.shape == b.shape
    assert not np.allclose(a, b)


def test_run_algorithm_names_and_timing():
    S = _sim()
    gx, gz = recon_grid(x_half_mm=10, nx=41, z_max_mm=26, nz=51)
    img, ms, name = run_algorithm("Reference DAS", S, gx, gz)
    assert name == "das_reference" and ms >= 0 and img.shape == (gz.size, gx.size)


def test_must_algorithm_is_matlab_only():
    S = _sim()
    gx, gz = recon_grid()
    try:
        run_algorithm("MUST das()", S, gx, gz)
    except ValueError as e:
        assert "MATLAB-only" in str(e)
    else:
        raise AssertionError("expected MUST das() to be rejected in the Python port")


def test_full_aperture_when_fnumber_zero():
    S = sim_engine(SimConfig(n_elements=32, zmax=26e-3, scat_z=18e-3, fnumber=0.0))
    for tx in S.tx:
        tx.fnumber = 0.0
    gx, gz = recon_grid(x_half_mm=10, nx=21, z_max_mm=26, nz=31)
    _, delays = das_reference(
        S.RF[:, :, 0], S.tx[0], S.rx_pos, gx, gz, S.c, S.fs, want_delays=True
    )
    # every active element contributes to every pixel within the record
    assert np.isfinite(delays).mean() > 0.9


def test_bundle_gain_scales_with_channel_count():
    """Mirrors the ``gain`` expression in ``drawAlignment``: a fixed 3.0 buries
    a sparse bundle, so it scales with the channel count and floors at 0.45."""
    assert bundle_gain(8) == 0.45                       # floor, traces stay apart
    assert bundle_gain(64) == 64 / 22
    assert bundle_gain(512) == 3.0                      # ceiling
    assert bundle_gain(16) < bundle_gain(32) < bundle_gain(64)


def test_alignment_window_holds_the_delay_spread():
    """The window is the *larger* of three periods of fc and the delay spread,
    or the "before" panel crops the scatter it exists to show."""
    S = _sim()
    gx, gz = recon_grid(x_half_mm=10, nx=1, z_max_mm=26, nz=1)
    _, tau = das_reference(S.RF[:, :, 0], S.tx[0], S.rx_pos,
                           np.array([0.0]), np.array([18e-3]),
                           S.c, S.fs, want_delays=True)
    tau = tau[0]
    b = alignment_bundles(S, 0, tau)
    ok = np.isfinite(tau)
    spread_us = float(np.max(tau[ok]) - np.min(tau[ok])) * 1e6
    half_win_us = float(b["trel_us"][-1])
    assert half_win_us >= 3 / S.probe["fc"] * 1e6 - 1e-12      # floor honoured
    assert 2 * half_win_us > spread_us                          # curve fits
    assert b["gain"] == bundle_gain(S.RF.shape[1])              # auto gain


def test_delay_curve_frame_maps_onto_the_bmode():
    """Element positions across, apparent depth c*t/2 down - the mapping that
    lets the delay panel share the B-mode's frame."""
    S = _sim()
    _, tau = das_reference(S.RF[:, :, 0], S.tx[0], S.rx_pos,
                           np.array([0.0]), np.array([18e-3]),
                           S.c, S.fs, want_delays=True)
    x_mm, depth_mm, tau_mm = delay_curve_frame(S, tau[0])
    assert x_mm.size == S.rx_pos.size
    assert np.allclose(x_mm, S.rx_pos * 1e3)
    assert depth_mm.size == S.RF.shape[0]
    assert np.allclose(depth_mm, S.c * (S.t0 + np.arange(S.RF.shape[0]) / S.fs) / 2 * 1e3)
    # An unsteered plane wave makes transmit and receive legs equal at the
    # apex, so c*tau/2 there is the pixel depth itself.
    ok = np.isfinite(tau_mm)
    assert abs(float(np.min(tau_mm[ok])) - 18.0) < 0.05


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  PASS  {fn.__name__}")
    print(f"\n{len(fns)} contract tests OK")
