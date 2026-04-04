# 07a-compute_phi_amp_max summary

- input A_unit grid: `results/04a-a-unit/A_unit_scan.csv`
- mass_min: `1e-27`
- mass_max: `None`
- evaluation redshift: `7.7`
- assumption: `rho_phi(z_eval) <= 1 * rho_DM(z_eval)`
- solver: `DOP853`, `rtol=1e-07`, `atol_rel=1e-10`, `dense_output=False`
- output csv: `phi_amp_max_scan.csv`

This is a phenomenological amplitude bound in the same field-amplitude units as the current ODE variable.
For WSL2 stability, prefer splitting heavy runs by mass range and using --resume.
