# 19-gbench_budget summary

This run separates the raw patchy response from the physical birefringence after inserting a benchmark ALP-photon coupling.

Coupling benchmark:
- Chandra / NGC1275 bound: `g_{a\gamma\gamma} \lesssim 1.4 - 4.0 x 10^-12 GeV^-1` for `m_a \lesssim 10^-12 eV`.
- Source: Berg et al., arXiv:1605.01043.

- matched mass: `5.878016e-27 eV`
- raw `A_tau_raw = phi_amp_max * A_unit = 2.020277e+14`
- template anchor: `tau=0.084`, `Delta y=19.0`, `Rbar=5.0`, `siglnR=0.693`
- surrogate: `z_re=11.166`, `R_eff=34.139 Mpc`, `L_peak=288.925`

Interpretation:
- `raw_patchy_alpha_vs_limit.png` is intentionally unphysical as an absolute observable because it omits the factor `g_{a\gamma}/2`.
- The `physical_*` plots insert `A_tau = (g/2) * A_tau_raw`.
- This assumes the current field normalization can be interpreted canonically in GeV units; if not, an additional normalization check is needed.
