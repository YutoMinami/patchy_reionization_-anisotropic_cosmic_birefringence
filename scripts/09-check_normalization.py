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

from patchy_reionization import H_of_a, X_REI, d_eta_d_tau, m_eV_to_omega_SI


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-unit-csv", type=Path, default=Path("results/04a-a-unit/A_unit_scan.csv"))
    parser.add_argument(
        "--phi-max-csv",
        type=Path,
        default=Path("results/07a-phi-amp-max-split/phi_amp_max_scan.csv"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/09-normalization-check"))
    parser.add_argument("--no-show", action="store_true")
    return parser


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    a_rows = read_csv(args.a_unit_csv)
    p_rows = read_csv(args.phi_max_csv)

    a_by_mass = {float(row["mass_eV"]): float(row["A_unit"]) for row in a_rows}
    a_eval = float(np.exp(X_REI))
    H_eval = float(H_of_a(a_eval))
    deta_dtau_eval = float(d_eta_d_tau(a_eval))

    rows: list[dict[str, float]] = []
    for row in p_rows:
        mass = float(row["mass_eV"])
        if mass not in a_by_mass:
            continue

        phi_eval = float(row["phi_unit_at_eval"])
        dphi_dx = float(row["dphi_dx_unit_at_eval"])
        dotphi_phys = float(row["dotphi_phys_unit_at_eval"])
        rho_unit = float(row["rho_unit_at_eval"])
        phi_amp_max = float(row["phi_amp_max"])
        A_unit_saved = float(a_by_mass[mass])

        dotphi_phys_from_dx = H_eval * dphi_dx
        dotphi_conf = a_eval * dotphi_phys
        A_unit_reconstructed = dotphi_conf * deta_dtau_eval
        omega = m_eV_to_omega_SI(mass)
        kinetic = 0.5 * dotphi_phys**2
        potential = 0.5 * (omega * phi_eval) ** 2
        kinetic_fraction = kinetic / rho_unit if rho_unit > 0.0 else float("nan")
        potential_fraction = potential / rho_unit if rho_unit > 0.0 else float("nan")
        rel_diff_A = (
            (A_unit_reconstructed - A_unit_saved) / A_unit_saved if A_unit_saved != 0.0 else float("nan")
        )
        rel_diff_dotphi = (
            (dotphi_phys_from_dx - dotphi_phys) / dotphi_phys if dotphi_phys != 0.0 else float("nan")
        )

        rows.append(
            {
                "mass_eV": mass,
                "A_unit_saved": A_unit_saved,
                "A_unit_reconstructed": A_unit_reconstructed,
                "A_unit_abs": abs(A_unit_saved),
                "rel_diff_A_unit": rel_diff_A,
                "phi_unit_at_eval": phi_eval,
                "dphi_dx_unit_at_eval": dphi_dx,
                "dotphi_phys_unit_at_eval": dotphi_phys,
                "dotphi_phys_from_dx": dotphi_phys_from_dx,
                "rel_diff_dotphi": rel_diff_dotphi,
                "omega_phi_abs": abs(omega * phi_eval),
                "kinetic_density": kinetic,
                "potential_density": potential,
                "rho_unit_at_eval": rho_unit,
                "kinetic_fraction": kinetic_fraction,
                "potential_fraction": potential_fraction,
                "phi_amp_max": phi_amp_max,
                "A_times_phi_amp_max": abs(A_unit_saved) * phi_amp_max,
            }
        )

    if not rows:
        raise RuntimeError("No overlapping masses found between the 04a and 07a CSV files.")

    rows.sort(key=lambda item: float(item["mass_eV"]))
    write_csv(args.output_dir / "normalization_check.csv", list(rows[0].keys()), rows)

    masses = np.array([row["mass_eV"] for row in rows])
    rel_diff_A = np.abs(np.array([row["rel_diff_A_unit"] for row in rows]))
    kinetic_fraction = np.array([row["kinetic_fraction"] for row in rows])
    potential_fraction = np.array([row["potential_fraction"] for row in rows])
    A_phi_max = np.array([row["A_times_phi_amp_max"] for row in rows])

    plt.figure(figsize=(7.5, 5.2))
    plt.loglog(masses, rel_diff_A)
    plt.xlabel(r"$m_a\ [{\rm eV}]$")
    plt.ylabel(r"$|A_{\rm recon}/A_{\rm saved} - 1|$")
    plt.title(r"$A_{\rm unit}$ reconstruction consistency")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output_dir / "aunit_reconstruction_diff.png", dpi=160, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(7.5, 5.2))
    plt.semilogx(masses, kinetic_fraction, label="kinetic fraction")
    plt.semilogx(masses, potential_fraction, label="potential fraction")
    plt.xlabel(r"$m_a\ [{\rm eV}]$")
    plt.ylabel(r"energy fraction")
    plt.title(r"Kinetic vs potential contribution at $z_{\rm eval}$")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "energy_fraction_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(7.5, 5.2))
    plt.loglog(masses, A_phi_max)
    plt.xlabel(r"$m_a\ [{\rm eV}]$")
    plt.ylabel(r"$|A_{\rm unit}|\,\phi_{\rm amp,max}$")
    plt.title(r"Maximum effective response scale")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output_dir / "aunit_times_phiampmax_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    max_rel_A = float(np.nanmax(rel_diff_A))
    min_kin = rows[int(np.nanargmin(kinetic_fraction))]
    max_kin = rows[int(np.nanargmax(kinetic_fraction))]
    max_Aphi = rows[int(np.nanargmax(A_phi_max))]

    summary_lines = [
        "# 09-check_normalization summary",
        "",
        f"- input A_unit csv: `{args.a_unit_csv}`",
        f"- input phi_amp_max csv: `{args.phi_max_csv}`",
        f"- max relative A_unit reconstruction error: `{max_rel_A:.6e}`",
        f"- kinetic-fraction minimum: `{min_kin['kinetic_fraction']:.6e}` at `m={min_kin['mass_eV']:.6e} eV`",
        f"- kinetic-fraction maximum: `{max_kin['kinetic_fraction']:.6e}` at `m={max_kin['mass_eV']:.6e} eV`",
        f"- max |A_unit| * phi_amp_max: `{max_Aphi['A_times_phi_amp_max']:.6e}` at `m={max_Aphi['mass_eV']:.6e} eV`",
        "",
        "This is a consistency check inside the current field normalization and phenomenological rho_phi proxy.",
    ]
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps(
            {
                "a_unit_csv": str(args.a_unit_csv),
                "phi_max_csv": str(args.phi_max_csv),
                "a_eval": a_eval,
                "H_eval": H_eval,
                "deta_dtau_eval": deta_dtau_eval,
                "max_relative_A_unit_error": max_rel_A,
            },
            indent=2,
        )
    )

    print(f"Saved outputs to: {args.output_dir}")
    print(f"max relative A_unit reconstruction error = {max_rel_A:.6e}")
    print(
        f"kinetic fraction range = [{min_kin['kinetic_fraction']:.6e}, {max_kin['kinetic_fraction']:.6e}]"
    )
    print(
        f"max |A_unit| * phi_amp_max = {max_Aphi['A_times_phi_amp_max']:.6e} "
        f"at m={max_Aphi['mass_eV']:.6e} eV"
    )

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
