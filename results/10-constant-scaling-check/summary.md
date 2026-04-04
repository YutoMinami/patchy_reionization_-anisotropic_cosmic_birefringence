# 10-check_constant_scaling summary

- input normalization csv: `results/09-normalization-check/normalization_check.csv`
- high-mass threshold: `1.000000e-31 eV`

- all-mass mean scaling ratio: `1.065644e+00` with rms fractional scatter `2.510246e-01`
- low-mass mean scaling ratio: `1.171085e+00` with rms fractional scatter `3.107253e-01`
- high-mass mean scaling ratio: `9.793743e-01` with rms fractional scatter `7.383075e-02`

Interpretation:
- If the rms scatter is small, a constant rescaling is a good approximation in that mass range.
- If the rms scatter is large, the mismatch is mass-dependent and should not be fixed by one global factor.
