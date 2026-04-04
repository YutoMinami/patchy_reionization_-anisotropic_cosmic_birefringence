# 18-deltatau_variation_budget summary

This run fixes the patchy template shape and varies only the optical-depth fluctuation amplitude proxy.

- matched mass: `5.878016e-27 eV`
- `A_tau_eff = A_unit * phi_amp_max = 2.020277e+14`
- template anchor: `tau=0.084`, `Delta y=19.0`, `Rbar=5.0`, `siglnR=0.693`
- surrogate: `z_re=11.166`, `R_eff=34.139 Mpc`, `L_peak=288.925`
- genuine benchmark: `0.30` times the anisotropic-CB limit

Interpretation:
- `Delta y` is fixed to `19.0` in this scan.
- `delta_tau` is interpreted as the optical-depth fluctuation amplitude proxy, normalized to the fiducial `delta_tau_fid = 0.084`.
- Therefore `C_L^{tau tau}` scales as `(delta_tau / delta_tau_fid)^2`.
- The patchy contribution to birefringence is `A_tau_eff^2 C_L^{tau tau}`.
- The total benchmark curve is `C_L^{alpha alpha,phi,bench} + A_tau_eff^2 C_L^{tau tau}`.
