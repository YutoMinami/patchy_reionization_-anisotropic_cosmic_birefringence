from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib

if "--no-show" in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patchy_reionization import (
    EV_TO_J,
    HBAR,
    H_of_a,
    MPC_TO_M,
    OMEGA_B,
    OMEGA_M,
    RHO_CRIT0,
    d_eta_d_tau,
    m_eV_to_omega_SI,
)


C_SI = 299792458.0
GEV_TO_J = 1.0e9 * EV_TO_J
GEV_INV_TO_M = HBAR * C_SI / GEV_TO_J
GEV4_TO_J_M3 = GEV_TO_J / (GEV_INV_TO_M**3)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matched-csv",
        type=Path,
        default=Path("results/11-matched-scan-global-split/matched_scan.csv"),
    )
    parser.add_argument("--mass", type=float, default=5.878016072274924e-27)
    parser.add_argument("--z-eval", type=float, default=7.7)
    parser.add_argument("--f-dm", type=float, default=1.0)
    parser.add_argument("--g-values", type=float, nargs="+", default=[1.4e-12, 4.0e-12])
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/21-natural-unit-normalization"),
    )
    parser.add_argument("--no-show", action="store_true")
    return parser


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def nearest_row(rows: list[dict[str, str]], target: float) -> dict[str, str]:
    masses = np.array([float(row["mass_eV"]) for row in rows])
    idx = int(np.argmin(np.abs(masses - target)))
    return rows[idx]


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

    matched = nearest_row(read_csv(args.matched_csv), args.mass)
    mass = float(matched["mass_eV"])
    phi_unit = float(matched["phi_unit_at_eval"])
    dphi_dx = float(matched["dphi_dx_unit_at_eval"])
    A_unit_old = float(matched["A_unit"])
    phi_amp_max_old = float(matched["phi_amp_max"])

    a_eval = 1.0 / (1.0 + args.z_eval)
    H_SI = float(H_of_a(a_eval))
    H_GeV = HBAR * H_SI / GEV_TO_J
    m_GeV = mass * 1.0e-9
    deta_dtau_m = d_eta_d_tau(a_eval)
    deta_dtau_GeVinv = deta_dtau_m / GEV_INV_TO_M

    # Canonical natural-unit reinterpretation:
    # phi in GeV, H in GeV, eta in GeV^-1.
    dotphi_nat = H_GeV * dphi_dx
    rho_unit_GeV4 = 0.5 * dotphi_nat**2 + 0.5 * m_GeV**2 * phi_unit**2

    rho_dm_SI = args.f_dm * rho_dm_of_a(a_eval)
    rho_dm_GeV4 = rho_dm_SI / GEV4_TO_J_M3
    phi_amp_max_nat_GeV = np.sqrt(rho_dm_GeV4 / rho_unit_GeV4)

    A_unit_nat_GeV = a_eval * H_GeV * dphi_dx * deta_dtau_GeVinv
    A_tau_raw_nat_GeV = phi_amp_max_nat_GeV * A_unit_nat_GeV

    rows: list[dict[str, float | str]] = []
    for g in args.g_values:
        a_tau_phys = 0.5 * g * A_tau_raw_nat_GeV
        rows.append(
            {
                "mass_eV": mass,
                "z_eval": args.z_eval,
                "a_eval": a_eval,
                "H_SI": H_SI,
                "H_GeV": H_GeV,
                "m_GeV": m_GeV,
                "phi_unit": phi_unit,
                "dphi_dx_unit": dphi_dx,
                "rho_unit_GeV4": rho_unit_GeV4,
                "rho_dm_GeV4": rho_dm_GeV4,
                "phi_amp_max_old_code_units": phi_amp_max_old,
                "phi_amp_max_nat_GeV": phi_amp_max_nat_GeV,
                "A_unit_old_mixed": A_unit_old,
                "A_unit_nat_GeV": A_unit_nat_GeV,
                "A_tau_raw_nat_GeV": A_tau_raw_nat_GeV,
                "g_GeV_inv": g,
                "A_tau_physical_dimless": a_tau_phys,
            }
        )

    write_csv(args.output_dir / "natural_unit_summary.csv", list(rows[0].keys()), rows)

    plt.figure(figsize=(7.8, 5.2))
    labels = [
        "phi_amp_max old",
        "phi_amp_max nat [GeV]",
        "A_unit old",
        "A_unit nat [GeV]",
        "A_tau raw nat [GeV]",
    ]
    values = [
        phi_amp_max_old,
        phi_amp_max_nat_GeV,
        abs(A_unit_old),
        abs(A_unit_nat_GeV),
        abs(A_tau_raw_nat_GeV),
    ]
    plt.bar(range(len(values)), values)
    plt.yscale("log")
    plt.xticks(range(len(values)), labels, rotation=20, ha="right")
    plt.title("Old mixed-unit quantities vs natural-unit reinterpretation")
    plt.tight_layout()
    plt.savefig(args.output_dir / "natural_unit_bar_summary.png", dpi=160, bbox_inches="tight")
    plt.close()

    summary_lines = [
        "# 21-check_natural_unit_normalization summary",
        "",
        "This run reinterprets the matched background quantities in canonical natural units.",
        "",
        "Conversions used:",
        "- `H_GeV = ħ H_SI / (1 GeV in J)`",
        "- `m_GeV = m_eV x 1e-9`",
        "- `dη/dτ` from the code is converted from meters to GeV^-1",
        "- `rho_DM` is converted from `J/m^3` to `GeV^4`",
        "",
        f"- matched mass: `{mass:.6e} eV`",
        f"- `phi_amp_max_old` = `{phi_amp_max_old:.6e}` (code units)",
        f"- `phi_amp_max_nat` = `{phi_amp_max_nat_GeV:.6e} GeV`",
        f"- `A_unit_old` = `{A_unit_old:.6e}` (mixed units)",
        f"- `A_unit_nat` = `{A_unit_nat_GeV:.6e} GeV`",
        f"- `A_tau_raw_nat = phi_amp_max_nat * A_unit_nat = {A_tau_raw_nat_GeV:.6e} GeV`",
        "",
    ]
    for row in rows:
        summary_lines.append(
            f"- `g={row['g_GeV_inv']:.2e} GeV^-1` -> "
            f"`A_tau_physical = {row['A_tau_physical_dimless']:.6e}`"
        )
    summary_lines += [
        "",
        "If these natural-unit quantities look much more reasonable than the old mixed-unit quantities,",
        "then the large `C_L^{alpha alpha}` problem was mainly a unit-system mismatch.",
    ]
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps({k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}, indent=2)
    )
    print(f"Saved outputs to: {args.output_dir}")


if __name__ == "__main__":
    main()
