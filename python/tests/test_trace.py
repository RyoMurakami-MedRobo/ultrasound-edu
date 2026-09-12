"""Execution traces must expose the beamformer's real channel increments."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from ultrasound_das import SimConfig, sim_engine
from ultrasound_das.das_reference import das_reference, envelope_z
from ultrasound_das.das_custom_template import das_custom_template


def test_trace_preserves_output_and_coherent_sum():
    s = sim_engine(SimConfig(n_elements=32, zmax=.028, scat_z=.018, fs_factor=4))
    z = np.linspace(.010, .026, 121)
    for fn in (das_reference, das_custom_template):
        for apod in ('rect', 'hann'):
            s.tx[0].rx_apod = apod
            args = (s.RF[:, :, 0], s.tx[0], s.rx_pos, np.array([0.0]), z, s.c, s.fs)
            baseline = fn(*args)
            image, delays, trace = fn(*args, want_trace=True)
            C = trace['contributions']
            np.testing.assert_array_equal(image, baseline)
            np.testing.assert_allclose(C.sum(axis=1), trace['coherent'][:, 0], atol=1e-12)
            np.testing.assert_allclose(envelope_z(np.cumsum(C, axis=1)[:, -1:]), image, atol=1e-12)
            assert np.all(C[~np.isfinite(delays)] == 0)


def test_trace_reflects_different_sampling():
    s = sim_engine(SimConfig(n_elements=16, zmax=.028, scat_z=.018, fnumber=0))
    z = np.linspace(.016, .020, 121)
    args = (s.RF[:, :, 0], s.tx[0], s.rx_pos, np.array([.001]), z, s.c, s.fs)
    a = das_reference(*args, want_trace=True)[2]['contributions']
    b = das_custom_template(*args, want_trace=True)[2]['contributions']
    assert np.max(np.abs(a-b)) > 1e-5


if __name__ == '__main__':
    test_trace_preserves_output_and_coherent_sum()
    test_trace_reflects_different_sampling()
    print('2 execution trace tests OK')
