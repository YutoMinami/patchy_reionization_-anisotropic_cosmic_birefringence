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
- `13-acb_reinterpretation_with_matched.py`
  Reinterprets the anisotropic-CB amplitude limit using the old toy Gaussian `C_L^{tau tau}` template and the matched `11` outputs.
- `13b-visualize_acb_reinterpretation.ipynb`
  Notebook for visualizing the `13` toy-template reinterpretation.
- `14-patchy_template_family_bound.py`
  Replaces the old toy Gaussian `C_L^{tau tau}` with a lightweight literature-inspired template family in `D_L^{tau tau}` and recomputes the simplified anisotropic-CB reinterpretation.
- `14b-visualize_patchy_template_family.ipynb`
  Notebook for plotting the template family and the allowed peak `D_L^{tau tau}` normalization from `14`.

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
11. `13-acb_reinterpretation_with_matched.py`
12. `13b-visualize_acb_reinterpretation.ipynb`
13. `14-patchy_template_family_bound.py`
14. `14b-visualize_patchy_template_family.ipynb`

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
./.venv/bin/python scripts/13-acb_reinterpretation_with_matched.py --no-show
./.venv/bin/python scripts/14-patchy_template_family_bound.py --no-show
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

`13-acb_reinterpretation_with_matched.py` writes:

- `results/13-acb-reinterpretation/acb_reinterpretation_scan.csv`
- `results/13-acb-reinterpretation/acb_reinterpretation_summary.csv`
- `results/13-acb-reinterpretation/claa_limit.png`
- `results/13-acb-reinterpretation/template_scale_limit_*.png`
- `results/13-acb-reinterpretation/run_config.json`
- `results/13-acb-reinterpretation/summary.md`

`14-patchy_template_family_bound.py` writes:

- `results/14-patchy-template-family/template_family_shapes.csv`
- `results/14-patchy-template-family/template_family_bound_summary.csv`
- `results/14-patchy-template-family/template_family_shapes.png`
- `results/14-patchy-template-family/dpeak_limit_vs_lpeak_*.png`
- `results/14-patchy-template-family/run_config.json`
- `results/14-patchy-template-family/summary.md`

This split workflow is intended to avoid losing everything when the heavy `A_unit` scan is interrupted.

## Notes

- `04a` is the expensive step.
- `07a` can also be expensive because it re-solves the ALP background for each mass.
- `11` is also expensive because it recomputes both `A_unit` and `phi_amp_max` from the same high-precision solution.
- `04b` should be comparatively light because it only post-processes saved `A_unit` values.
- `13` and `14` are lightweight reinterpretation steps that reuse saved matched outputs.
- On WSL2, prefer running `07a` in mass chunks with `--mass-min`, `--mass-max`, and `--resume`.
- On WSL2, prefer running `11` in mass chunks or exact-mass jobs using `--mass-min`, `--mass-max`, or `--masses`.
- Core numerical logic should live in `../patchy_reionization.py`, not be reimplemented separately in each script.

- `23-check_osc_scale.py`
  Caveat 2 verification Step 1: compute the ALP oscillation comoving wavelength
  `lambda_osc(m_a)` at `z_rei` and compare with the bubble scale `R_eff` and the
  visibility-function width `Delta_chi_vis`. Outputs `N_osc` (number of oscillations
  in the visibility window) to diagnose whether the thin-shell approximation is accurate.
- `23b-visualize_osc_scale.ipynb`
  Notebook for visualizing the `lambda_osc` vs scale comparison from `23`.

## Later scripts

- `15-a_reff_surrogate_bound.py`
  Dvorkin-Smith の `A + R_{\rm eff}` surrogate で `C_L^{\tau\tau}` template を作り、anisotropic-CB reinterpretation を更新する。
- `15b-visualize_a_reff_surrogate.ipynb`
  `15` の surrogate shape と制限図を見る notebook。
- `16-check_ds_figure_reproduction.py`
  Dvorkin-Smith の図を surrogate で目視再現できるかを試す最初のチェック。
- `16b-visualize_ds_figure_reproduction.ipynb`
  `16` の再現図と paper image を並べて見る notebook。
- `17-check_ds_figure_reproduction_corrected.py`
  Fig.16 / Fig.18 のパラメータと軸を補正した図再現チェック。
- `17b-visualize_ds_figure_reproduction_corrected.ipynb`
  `17` の corrected figure check を見る notebook。
- `18-deltatau_variation_budget.py`
  `\Delta\tau` の variation だけで `C_L^{\tau\tau}` と birefringence budget を描く簡約チェック。
- `18b-visualize_deltatau_variation_budget.ipynb`
  `18` の図と表をまとめて見る notebook。
- `19-gbench_budget.py`
  Chandra の `g_{a\gamma\gamma}` benchmark を入れて raw / physical の違いを切り分ける。
- `19b-visualize_gbench_budget.ipynb`
  `19` の raw / physical budget を見る notebook。
- `20-check_atau_normalization.py`
  `A_\tau` が大きすぎる原因を field-unit factor の観点から点検する。
- `20b-visualize_atau_normalization.ipynb`
  `20` の正規化チェックを図で確認する notebook。
- `21-check_natural_unit_normalization.py`
  matched quantity を canonical natural units に写して `A_\tau` の大きさを再評価する。
- `21b-visualize_natural_unit_normalization.ipynb`
  `21` の natural-unit summary を確認する notebook。
- `22-natural_unit_budget.py`
  `21` の natural-unit `A_\tau` を使って `D_L^{\alpha\alpha}` budget を描き直す。
- `22b-visualize_natural_unit_budget.ipynb`
  `22` の自然単位 budget 図を確認する notebook。

## Current takeaway after `21/22`

- huge な `C_L^{\alpha\alpha}` は、主に mixed-unit な `A_\tau` をそのまま observable に入れていたことに由来していた。
- `21` の自然単位化では、matched mass `m_a = 5.878016e-27 eV` で `A_\tau \simeq 0.157` for `g=1.4e-12 GeV^{-1}`、`A_\tau \simeq 0.448` for `g=4.0e-12 GeV^{-1}` になった。
- `22` では patchy contribution の最大値が anisotropic-CB limit に対して、`g=1.4e-12` で約 `1.2%`、`g=4.0e-12` でも約 `10%` に収まり、前の absurdly large な budget は unit-system mismatch が主因だったことが強く示唆される。

