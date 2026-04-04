# 21-check_natural_unit_normalization summary

This run reinterprets the matched background quantities in canonical natural units.

Conversions used:
- `H_GeV = ħ H_SI / (1 GeV in J)`
- `m_GeV = m_eV x 1e-9`
- `dη/dτ` from the code is converted from meters to GeV^-1
- `rho_DM` is converted from `J/m^3` to `GeV^4`

- matched mass: `5.878016e-27 eV`
- `phi_amp_max_old` = `1.449744e+04` (code units)
- `phi_amp_max_nat` = `4.823366e+09 GeV`
- `A_unit_old` = `1.393540e+10` (mixed units)
- `A_unit_nat` = `4.648350e+01 GeV`
- `A_tau_raw_nat = phi_amp_max_nat * A_unit_nat = 2.242069e+11 GeV`

- `g=1.40e-12 GeV^-1` -> `A_tau_physical = 1.569449e-01`
- `g=4.00e-12 GeV^-1` -> `A_tau_physical = 4.484139e-01`

If these natural-unit quantities look much more reasonable than the old mixed-unit quantities,
then the large `C_L^{alpha alpha}` problem was mainly a unit-system mismatch.
