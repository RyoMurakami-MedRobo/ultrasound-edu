"""Ultrasound DAS Beamforming Explorer -- Python port.

Python port of the MATLAB educational/verification tool for Delay-And-Sum
beamforming. MATLAB is the source of truth; see ``CLAUDE.md`` at the repo root
for the rule that keeps this port in lock-step.

Scope note: MUST (Matlab UltraSound Toolbox) has no Python equivalent, so this
port implements the analytic **mock** RF backend only. The ``MUST das()``
algorithm and ``validation/validate_must_vs_mock.m`` stay MATLAB-only.
"""

from .das_custom_template import das_custom_template
from .das_reference import das_reference, envelope_z, tx_arrival_time
from .metrics import compare_images, fwhm, lateral_metrics, to_db
from .pipeline import ALGORITHMS, recon_grid, run_algorithm, run_all_tx
from .sim_engine import SimConfig, SimResult, TxInfo, pulse_sigma, sim_engine, tx_law
from .wave_animator import (
    WavePropagation,
    alignment_bundles,
    blue_white_red,
    delay_curve_image,
)

__version__ = "1.0.0"

__all__ = [
    "SimConfig",
    "SimResult",
    "TxInfo",
    "sim_engine",
    "tx_law",
    "pulse_sigma",
    "das_reference",
    "das_custom_template",
    "tx_arrival_time",
    "envelope_z",
    "compare_images",
    "fwhm",
    "lateral_metrics",
    "to_db",
    "run_algorithm",
    "run_all_tx",
    "recon_grid",
    "ALGORITHMS",
    "WavePropagation",
    "delay_curve_image",
    "alignment_bundles",
    "blue_white_red",
]
