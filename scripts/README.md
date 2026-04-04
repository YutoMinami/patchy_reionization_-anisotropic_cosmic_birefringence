# scripts/README.md

## Purpose

This directory contains the reproducible analysis entry points for the project.

The notebook is kept for figure-driven exploratory work.
The Python scripts are kept for reproducible, comparable runs.

Important:

- Once a numbered script or notebook has been used for a recorded result, keep it frozen for reproducibility.
- If a correction or better-matched rerun is needed, add a new continued number or explicit versioned suffix instead of rewriting the old numbered artifact.

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
- `07a-compute_phi_amp_max.py`
  Computes the phenomenological upper amplitude `phi_amp_max(m_a)` and writes partial CSV output as it goes. Prefer split mass ranges on WSL2.
- `07b-visualize_phi_amp_max.ipynb`
  Notebook for plotting `phi_amp_max(m_a)` together with `phi_needed(m_a)`. It now defaults to the split `07a` outputs.
- `11-recompute_matched_scan.py`
  Recomputes `A_unit` and `phi_amp_max` from the same high-precision background solution. Prefer split or exact-mass runs on WSL2.
- `12-ratio_with_matched_scan.py`
  Recomputes `phi_needed / phi_amp_max` using the matched `11` outputs rather than the older mixed `04a/07a` combination.
- `12b-visualize_matched_ratio.ipynb`
  Notebook for plotting the matched-rerun ratio from `12`.

## Recommended order

1. `02-check_linearity.py`
2. `03-check_convergence.py`
3. `01-check_feasibility.py`
4. `04a-scan_a_unit.py`
5. `04b-scan_rtau_from_aunit.py`
6. `07a-compute_phi_amp_max.py`
7. `07b-visualize_phi_amp_max.ipynb`
8. `11-recompute_matched_scan.py`
9. `12-ratio_with_matched_scan.py`
10. `12b-visualize_matched_ratio.ipynb`

## Typical commands

```bash
./.venv/bin/python scripts/02-check_linearity.py
./.venv/bin/python scripts/03-check_convergence.py --masses 1e-32 1e-30 1e-28 1e-27 1e-26
./.venv/bin/python scripts/04a-scan_a_unit.py --m-max 1e-26 --num-mass 40 --no-show
./.venv/bin/python scripts/04b-scan_rtau_from_aunit.py --a-unit-csv results/04a-a-unit/A_unit_scan.csv --include-cutoff --no-show
./.venv/bin/python scripts/07a-compute_phi_amp_max.py --mass-max 1e-30 --resume --no-show --output-dir results/07a-phi-amp-max-split
./.venv/bin/python scripts/07a-compute_phi_amp_max.py --mass-min 1e-30 --mass-max 1e-28 --resume --no-show --output-dir results/07a-phi-amp-max-split
./.venv/bin/python scripts/07a-compute_phi_amp_max.py --mass-min 1e-28 --mass-max 1e-27 --resume --no-show --output-dir results/07a-phi-amp-max-split
./.venv/bin/python scripts/07a-compute_phi_amp_max.py --mass-min 1e-27 --resume --no-show --output-dir results/07a-phi-amp-max-split
./.venv/bin/python scripts/11-recompute_matched_scan.py --mass-max 1e-30 --resume --no-show --output-dir results/11-matched-scan-global-split
./.venv/bin/python scripts/11-recompute_matched_scan.py --mass-min 1e-30 --mass-max 1e-28 --resume --no-show --output-dir results/11-matched-scan-global-split
./.venv/bin/python scripts/11-recompute_matched_scan.py --mass-min 1e-28 --resume --no-show --output-dir results/11-matched-scan-global-split
./.venv/bin/python scripts/12-ratio_with_matched_scan.py --no-show
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

`07a-compute_phi_amp_max.py` writes:

- `results/07a-phi-amp-max-split/phi_amp_max_scan.csv`
- `results/07a-phi-amp-max-split/phi_amp_max_vs_m.png`
- `results/07a-phi-amp-max-split/run_config.json`
- `results/07a-phi-amp-max-split/summary.md`

`11-recompute_matched_scan.py` writes:

- `results/11-matched-scan-global-split/matched_scan.csv`
- `results/11-matched-scan-global-split/A_unit_vs_m.png`
- `results/11-matched-scan-global-split/phi_amp_max_vs_m.png`
- `results/11-matched-scan-global-split/run_config.json`
- `results/11-matched-scan-global-split/summary.md`

`12-ratio_with_matched_scan.py` writes:

- `results/12-ratio-with-matched-scan/matched_ratio_scan.csv`
- `results/12-ratio-with-matched-scan/matched_ratio_vs_m.png`
- `results/12-ratio-with-matched-scan/matched_breakdown.png`
- `results/12-ratio-with-matched-scan/run_config.json`
- `results/12-ratio-with-matched-scan/summary.md`

This split workflow is intended to avoid losing everything when the heavy `A_unit` scan is interrupted.

## Notes

- `04a` is the expensive step.
- `07a` can also be expensive because it re-solves the ALP background for each mass.
- `11` is also expensive because it recomputes both `A_unit` and `phi_amp_max` from the same high-precision solution.
- `04b` should be comparatively light because it only post-processes saved `A_unit` values.
- On WSL2, prefer running `07a` in mass chunks with `--mass-min`, `--mass-max`, and `--resume`.
- On WSL2, prefer running `11` in mass chunks or exact-mass jobs using `--mass-min`, `--mass-max`, or `--masses`.
- Core numerical logic should live in `../patchy_reionization.py`, not be reimplemented separately in each script.
