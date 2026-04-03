# scripts/README.md

## Purpose

This directory contains the reproducible analysis entry points for the project.

The notebook is kept for figure-driven exploratory work.
The Python scripts are kept for reproducible, comparable runs.

## Files

- `01-check_feasibility.ipynb`
  Interactive notebook with the initial feasibility scan, linearity check, and convergence check.
- `01-check_feasibility.py`
  Script version of the initial mass scan and toy `R_tau(L)` evaluation.
- `02-check_linearity.py`
  Checks whether the response coefficient scales linearly with `phi_ini`.
- `03-check_convergence.py`
  Compares solver and tolerance settings at representative masses.
- `04a-scan_a_unit.py`
  Long-running scan for `A_unit(m_a)` only. Writes partial CSV output as it goes.
- `04b-scan_rtau_from_aunit.py`
  Lightweight post-processing step that reads an `A_unit` CSV and computes `R_{\tau,\max}^{unit}(m_a)` and `phi_needed(m_a)`.
- `04-scan_rtau_needed.py`
  Older combined script kept for reference. Prefer `04a` + `04b` when memory or stability matters.

## Recommended order

1. `02-check_linearity.py`
2. `03-check_convergence.py`
3. `01-check_feasibility.py`
4. `04a-scan_a_unit.py`
5. `04b-scan_rtau_from_aunit.py`

## Typical commands

```bash
./.venv/bin/python scripts/02-check_linearity.py
./.venv/bin/python scripts/03-check_convergence.py --masses 1e-32 1e-30 1e-28 1e-27 1e-26
./.venv/bin/python scripts/04a-scan_a_unit.py --m-max 1e-26 --num-mass 40 --no-show
./.venv/bin/python scripts/04b-scan_rtau_from_aunit.py --a-unit-csv results/04a-a-unit/A_unit_scan.csv --include-cutoff --no-show
```

## Output files

`04a-scan_a_unit.py` writes:

- `results/04a-a-unit/A_unit_scan.csv`
- `results/04a-a-unit/A_unit_vs_m.png`
- `results/04a-a-unit/run_config.json`
- `results/04a-a-unit/summary.md`

`04b-scan_rtau_from_aunit.py` writes:

- `results/04b-rtau-needed/Rtau_phi_needed_scan.csv`
- `results/04b-rtau-needed/Rtau_max_unit_vs_m.png`
- `results/04b-rtau-needed/phi_needed_target_*.png`
- `results/04b-rtau-needed/run_config.json`
- `results/04b-rtau-needed/summary.md`

This split workflow is intended to avoid losing everything when the heavy `A_unit` scan is interrupted.

## Notes

- `04a` is the expensive step.
- `04b` should be comparatively light because it only post-processes saved `A_unit` values.
- Core numerical logic should live in `../patchy_reionization.py`, not be reimplemented separately in each script.
