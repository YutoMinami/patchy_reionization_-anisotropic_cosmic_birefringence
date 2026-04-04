# 16-check_ds_figure_reproduction summary

This run is a human-facing reproduction check for the Dvorkin-Smith-inspired surrogate.

What to inspect by eye:
- In the Fig.16-like panel, changing `(tau, Delta y)` should mainly move the overall amplitude.
- In the Fig.18-like panel, changing `(Rbar, sigma_lnR)` at roughly fixed `R_eff` should keep the peak position approximately fixed.

- `fig16_like`: tau=0.068, Delta y=19.000, A=0.0700, R_eff=34.139 Mpc, L_peak=279.41, D_peak=2.500e-05
- `fig16_like`: tau=0.110, Delta y=14.700, A=0.0863, R_eff=34.139 Mpc, L_peak=300.08, D_peak=3.081e-05
- `fig16_like`: tau=0.110, Delta y=19.000, A=0.1119, R_eff=34.139 Mpc, L_peak=300.06, D_peak=3.996e-05
- `fig18_like`: tau=0.068, Delta y=19.000, A=0.0700, R_eff=34.139 Mpc, L_peak=279.41, D_peak=2.500e-05
- `fig18_like`: tau=0.068, Delta y=19.000, A=0.0700, R_eff=19.346 Mpc, L_peak=493.06, D_peak=2.500e-05
- `fig18_like`: tau=0.068, Delta y=19.000, A=0.0700, R_eff=15.620 Mpc, L_peak=610.66, D_peak=2.500e-05
