# 23-check_osc_scale summary

Step 1 of Caveat 2 verification: compare lambda_osc(m_a) to R_eff and Delta_chi_vis.

## Visibility function (tanh model)
- `z_rei = 7.7`, `delta_y = 19.0`
- `chi_rei = 9036.7 Mpc`
- `sigma_chi = 2507.3 Mpc`  (1-sigma comoving width of g(chi))
- `fwhm_chi = 5904.3 Mpc`
- `tau_total = 0.0508`
- `delta_chi_vis (2 sigma, full tanh) = 5014.7 Mpc`

## Patchy reionization epoch width
- `delta_z_patchy = 3.0`
- `dchi/dz at z_rei = 307.9 Mpc/dz`
- `delta_chi_patchy = 923.8 Mpc`

## Bubble scale reference
- `R_eff (fiducial, from run 22) = 34.1 Mpc`
- `R_eff range = 5.0--50.0 Mpc`

## At the reference mass
- `m_best = 5.878016e-27 eV`
- `lambda_osc = 0.0595 Mpc = 59.5 kpc`
- `N_osc (full vis)   = 84321.1`
- `N_osc (patchy)     = 15533.8`

## Crossover mass (N_osc = 1, patchy epoch)
- `m_crossover = 3.784e-31 eV`

## Spatial resonance mass (lambda_osc = R_eff)
- `m_res (fiducial, R_eff=34.139 Mpc) = 1.024e-29 eV`
- `m_res range (5.0--50.0 Mpc) = 6.99e-30--6.99e-29 eV`
  m_res is ~3 orders of magnitude below m_best.

## Interpretation
- At `m_best`, `lambda_osc << R_eff`: no spatial resonance at the reference mass.
- Spatial resonance occurs at `m_res ~ 1e-29 eV`; A_unit there is comparable to m_best.
- `N_osc >> 1` at `m_best`: dominant concern is temporal phase averaging.
  The full visibility-weighted A_eff is quantified in script 24.
