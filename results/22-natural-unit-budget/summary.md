# 22-natural_unit_budget summary

This run rebuilds the birefringence budget using the natural-unit quantities from `21`.

- matched mass: `5.878016e-27 eV`
- template anchor: `tau=0.084`, `Delta y=19.0`, `Rbar=5.0`, `siglnR=0.693`
- surrogate: `z_re=11.166`, `R_eff=34.139 Mpc`, `L_peak=288.925`
- genuine benchmark: `0.30` times the anisotropic-CB limit

Natural-unit inputs:
- `g=1.40e-12 GeV^-1`: `phi_amp_max_nat=4.823366e+09 GeV`, `A_unit_nat=4.648350e+01 GeV`, `A_tau_physical=1.569449e-01`
- `g=4.00e-12 GeV^-1`: `phi_amp_max_nat=4.823366e+09 GeV`, `A_unit_nat=4.648350e+01 GeV`, `A_tau_physical=4.484139e-01`

Interpretation:
- Unlike `18` and the raw part of `19`, these curves use the natural-unit `A_tau` from `21`.
- The main comparison is done in `D_L^{alpha alpha}` so the plotted amplitudes can be compared to the paper-style `D_L` language directly.
- If these curves now sit near or below the anisotropic-CB limit, the previous huge amplitudes were mainly a unit mismatch.
