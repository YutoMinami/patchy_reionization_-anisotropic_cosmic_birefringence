# 26-natural_unit_aeff_mres summary

This run combines the fine `25` visibility-averaged scan with the natural-unit reinterpretation logic from `21`.

- source CSV: `results/25-visibility-mres-fine/visibility_aeff_mres_scan.csv`
- `z_eval = 7.7`
- `f_DM = 1.0`

Interpretation:
- `A_tau_unit_physical` is the old thin-shell physical coefficient in natural units.
- `A_tau_eff_physical` replaces it by the visibility-weighted `A_eff` result from `25`.
- This is the quantity that should be used for the next budget re-evaluation near `m_res`.
- `g=1.40e-12 GeV^-1`: max `|A_tau_eff_physical| = 2.316157e-02` at `m/m_res = 0.574349`
- `g=4.00e-12 GeV^-1`: max `|A_tau_eff_physical| = 6.617593e-02` at `m/m_res = 0.574349`
