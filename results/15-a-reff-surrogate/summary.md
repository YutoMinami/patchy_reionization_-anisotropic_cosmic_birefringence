# 15-a_reff_surrogate_bound summary

This run replaces the free unit-peak normalization of `14` with a Dvorkin-Smith-inspired surrogate.

Approximation choices:
- The observable reionization-amplitude combination is modeled by Eq. (77)-inspired `A(tau, Delta y)` using a tanh-like ionization history.
- The characteristic angular scale is modeled through `R_eff = Rbar * exp(4 sigma_lnR^2)` from Eq. (78).
- The shape remains a lognormal bump in `D_L^{tau tau}`, centered at `L_peak ~= chi(z_re)/R_eff`.
- The surrogate amplitude is calibrated to `A_fid=0.070` and `D_peak_fid=2.500e-05` at `(tau, Delta y)=(0.068, 19.0)`.

- `m=2.030918e-27 eV`: strongest tension `scale<4.707526e-29` for `(tau=0.110, Delta y=19, b=2, Rbar=8, siglnR=0.693)`; weakest tension `scale<1.206087e-21` for `(tau=0.110, Delta y=14.7, b=1, Rbar=5, siglnR=0.409)`.
- `m=5.878016e-27 eV`: strongest tension `scale<4.626459e-29` for `(tau=0.110, Delta y=19, b=2, Rbar=8, siglnR=0.693)`; weakest tension `scale<1.185318e-21` for `(tau=0.110, Delta y=14.7, b=1, Rbar=5, siglnR=0.409)`.
- `m=1.000000e-26 eV`: strongest tension `scale<1.642466e-27` for `(tau=0.110, Delta y=19, b=2, Rbar=8, siglnR=0.693)`; weakest tension `scale<4.208065e-20` for `(tau=0.110, Delta y=14.7, b=1, Rbar=5, siglnR=0.409)`.
