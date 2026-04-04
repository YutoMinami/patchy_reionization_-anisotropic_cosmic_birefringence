# 20-check_atau_normalization summary

This run asks whether the current ODE field normalization can be directly multiplied by `g_{a\gamma}`.

Definitions:
- code/raw quantity: `A_tau_raw = phi_amp_max * A_unit`
- naive physical quantity if the code field were already canonical: `A_tau = (g/2) * A_tau_raw`
- more generally, if `phi_phys = U_phi * phi_code`, then `A_tau = (g/2) * U_phi * A_tau_raw`

- matched mass: `5.878016e-27 eV`
- `A_unit = 1.393540e+10`
- `phi_amp_max` (code units) = `1.449744e+04`
- `A_tau_raw` (code units) = `2.020277e+14`
- template anchor: `tau=0.084`, `Delta y=19.0`, `delta_tau=0.084`, `R_eff=34.139 Mpc`
- peak at `L=289` with `C_tau_tau=2.174112e-09` and `C_alphaalpha_limit=7.496940e-09`

Interpretation:
- If `U_phi_limit_peak` is extremely small or extremely large, then the present code-field normalization is unlikely to be the canonical field normalization that should be used with `g` directly.
- This does not by itself prove the physics is wrong; it isolates a missing unit-conversion / canonical-normalization step.

- `g=1.40e-12 GeV^-1`: `A_tau_naive=1.414194e+02`, `C_peak_naive=4.348105e-05`, `U_phi_limit_peak=1.313082e-02`
- `g=4.00e-12 GeV^-1`: `A_tau_naive=4.040555e+02`, `C_peak_naive=3.549474e-04`, `U_phi_limit_peak=4.595788e-03`
