from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patchy_reionization import (
    EV_TO_J,
    HBAR,
    H_of_a,
    OMEGA_B,
    OMEGA_M,
    RHO_CRIT0,
)

C_SI = 299792458.0
GEV_TO_J = 1.0e9 * EV_TO_J
GEV_INV_TO_M = HBAR * C_SI / GEV_TO_J
GEV4_TO_J_M3 = GEV_TO_J / (GEV_INV_TO_M**3)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--visibility-csv",
        type=Path,
        default=Path("results/25-visibility-mres-fine/visibility_aeff_mres_scan.csv"),
    )
    p.add_argument("--z-eval", type=float, default=7.7)
    p.add_argument("--f-dm", type=float, default=1.0)
    p.add_argument("--g-values", type=float, nargs="+", default=[1.4e-12, 4.0e-12])
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/26-natural-unit-aeff-mres"),
    )
    return p


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, float | str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def rho_dm_of_a(a: float) -> float:
    omega_dm = OMEGA_M - OMEGA_B
    return omega_dm * RHO_CRIT0 * a**-3


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    vis_rows = read_csv(args.visibility_csv)
    a_eval = 1.0 / (1.0 + args.z_eval)
    H_SI = float(H_of_a(a_eval))
    H_GeV = HBAR * H_SI / GEV_TO_J
    rho_dm_GeV4 = args.f_dm * rho_dm_of_a(a_eval) / GEV4_TO_J_M3

    rows: list[dict[str, float | str]] = []
    for row in vis_rows:
        mass = float(row["mass_eV"])
        phi_unit = float(row["phi_unit_at_zrei"])
        dphi_dx = float(row["dphi_dx_unit_at_zrei"])
        aeff_over_aunit = float(row["Aeff_over_Aunit"])
        a_unit_old = float(row["A_unit_thinshell"])
        a_eff_old = float(row["A_eff_visibility"])

        m_GeV = mass * 1.0e-9
        deta_dtau_m = float((a_unit_old / (a_eval * H_SI * dphi_dx)) if dphi_dx != 0.0 else np.nan)
        deta_dtau_GeVinv = deta_dtau_m / GEV_INV_TO_M

        dotphi_nat = H_GeV * dphi_dx
        rho_unit_GeV4 = 0.5 * dotphi_nat**2 + 0.5 * m_GeV**2 * phi_unit**2
        phi_amp_max_nat_GeV = np.sqrt(rho_dm_GeV4 / rho_unit_GeV4)

        a_unit_nat = a_eval * H_GeV * dphi_dx * deta_dtau_GeVinv
        a_eff_nat = aeff_over_aunit * a_unit_nat

        for g in args.g_values:
            a_tau_unit_phys = 0.5 * g * phi_amp_max_nat_GeV * a_unit_nat
            a_tau_eff_phys = 0.5 * g * phi_amp_max_nat_GeV * a_eff_nat
            rows.append(
                {
                    "mass_eV": mass,
                    "m_res_fid_eV": float(row["m_res_fid_eV"]),
                    "mass_over_mres": float(row["mass_over_mres"]),
                    "g_GeV_inv": g,
                    "phi_unit": phi_unit,
                    "dphi_dx_unit": dphi_dx,
                    "rho_unit_GeV4": rho_unit_GeV4,
                    "rho_dm_GeV4": rho_dm_GeV4,
                    "phi_amp_max_nat_GeV": phi_amp_max_nat_GeV,
                    "A_unit_old_mixed": a_unit_old,
                    "A_eff_old_mixed": a_eff_old,
                    "Aeff_over_Aunit": aeff_over_aunit,
                    "A_unit_nat": a_unit_nat,
                    "A_eff_nat": a_eff_nat,
                    "A_tau_unit_physical": a_tau_unit_phys,
                    "A_tau_eff_physical": a_tau_eff_phys,
                    "lambda_osc_Mpc": float(row["lambda_osc_Mpc"]),
                    "N_osc_patchy": float(row["N_osc_patchy"]),
                }
            )

    write_csv(args.output_dir / "natural_unit_aeff_mres_scan.csv", list(rows[0].keys()), rows)

    summary_lines = [
        "# 26-natural_unit_aeff_mres summary",
        "",
        "This run combines the fine `25` visibility-averaged scan with the natural-unit reinterpretation logic from `21`.",
        "",
        f"- source CSV: `{args.visibility_csv}`",
        f"- `z_eval = {args.z_eval}`",
        f"- `f_DM = {args.f_dm}`",
        "",
        "Interpretation:",
        "- `A_tau_unit_physical` is the old thin-shell physical coefficient in natural units.",
        "- `A_tau_eff_physical` replaces it by the visibility-weighted `A_eff` result from `25`.",
        "- This is the quantity that should be used for the next budget re-evaluation near `m_res`.",
    ]
    for g in args.g_values:
        subset = [r for r in rows if abs(float(r["g_GeV_inv"]) - g) < 1e-30]
        peak = max(subset, key=lambda r: abs(float(r["A_tau_eff_physical"])))
        summary_lines.append(
            f"- `g={g:.2e} GeV^-1`: max `|A_tau_eff_physical| = {abs(float(peak['A_tau_eff_physical'])):.6e}` at "
            f"`m/m_res = {float(peak['mass_over_mres']):.6f}`"
        )

    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps({k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}, indent=2)
    )
    print(f"Saved outputs to: {args.output_dir}")


if __name__ == "__main__":
    main()
