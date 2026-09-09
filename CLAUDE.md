# Working in this repository

This project ships the **same tool twice**: the original MATLAB implementation
at the repo root, and a Python port under [`python/`](python/). They must not
drift.

## MATLAB is the source of truth

Every feature and every formula originates in MATLAB. The Python package
mirrors it. When you change behaviour, change MATLAB first, then port the same
change to Python **in the same commit / PR**.

### The mirror map

| MATLAB (root)            | Python (`python/ultrasound_das/`)              |
|--------------------------|-----------------------------------------------|
| `sim_engine.m`           | `sim_engine.py`                               |
| `das_reference.m`        | `das_reference.py`                            |
| `das_custom_template.m`  | `das_custom_template.py`                      |
| `wave_animator.m`        | `wave_animator.py` (numeric) + `gui.py` (drawing) |
| `main_gui.m`             | `gui.py` + `pipeline.py` + `metrics.py`       |
| `setup_must.m`, `validation/`, `paper/` | *no Python equivalent* (MUST is MATLAB-only) |

### Definition of done for any change to the numeric core

A change to `sim_engine.m`, `das_reference.m`, `das_custom_template.m`, or the
metric/pipeline helpers in `main_gui.m` is **not complete** until:

1. The mirrored Python file is edited to match.
2. If any number moved, the golden fixtures are regenerated:
   ```
   matlab -batch "run('python/parity/dump_reference.m')"
   ```
   (or run `python/parity/dump_reference.m` from the MATLAB GUI / the MATLAB
   MCP server). Commit the regenerated `python/parity/fixtures/*.mat`.
3. The parity tests pass:
   ```
   cd python && python tests/test_parity.py && python tests/test_contract.py
   ```
   (or `pytest` if installed). Parity tolerance is `rtol = 1e-9` on RF, delays
   and the linear envelope. **If a test only passes at a looser tolerance, that
   is a structural port bug — fix the port, do not loosen the test.**

### Known deliberate MATLAB↔Python divergences

Keep these in mind so a future edit does not "fix" them back:

- **MUST backend / `MUST das()` / `validate_must_vs_mock.m`** — MATLAB only.
  Python implements the analytic *mock* RF backend exclusively. `SimConfig`
  still carries `force_mock` and reports `backend = 'mock (...)'` for API
  parity.
- **`wave_animator.m` deliberately duplicates its own `txArrivalTime`**
  (`wave_animator.m:311`) to stay independent of the beamformer. The Python
  port instead reuses `das_reference.tx_arrival_time`. If the MATLAB duplicate
  is ever changed *without* changing `das_reference.m`, reproduce that change
  in `wave_animator.py` on purpose.
- **GUI toolkit** — `uifigure`/`uitabgroup`/`uitable` in MATLAB;
  `tk.Tk`/`ttk.Notebook`/`ttk.Treeview` in Python. Match structure and
  behaviour, not widget-for-widget API.
- **`wave_animator.py` splits compute from drawing**; `wave_animator.m` mixes
  them. Numeric changes → `wave_animator.py`; drawing changes → `gui.py`.
- **Anechoic-cyst preset RNG** — `rng(0)` in MATLAB vs.
  `np.random.default_rng(0)` in Python produce *different* scatterer clouds.
  This preset is intentionally excluded from parity tests. Do not try to make
  the streams match.

## Quick commands

```bash
# MATLAB tool
matlab -r main_gui

# Python tool
cd python && python -m ultrasound_das            # GUI
cd python && python tests/test_parity.py         # parity vs MATLAB fixtures
cd python && python tests/test_contract.py       # interface contract (no MATLAB)
```
