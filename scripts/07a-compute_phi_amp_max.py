from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
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
    OMEGA_B,
    OMEGA_M,
    RHO_CRIT0,
    H_of_a,
    SolveConfig,
    evaluate_state_at_x,
    m_eV_to_omega_SI,
    solve_phi_background,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-unit-csv", type=Path, default=Path("results/04a-a-unit/A_unit_scan.csv"))
    parser.add_argument("--mass-min", type=float, default=None)
    parser.add_argument("--mass-max", type=float, default=None)
    parser.add_argument("--z-eval", type=float, default=7.7)
    parser.add_argument(
        "--f-dm",
        type=float,
        default=1.0,
        help="Assume rho_phi(z_eval) <= f_dm * rho_DM(z_eval).",
    )
    parser.add_argument("--z-ini", type=float, default=1.0e7)
    parser.add_argument("--method", type=str, default="DOP853")
    parser.add_argument("--rtol", type=float, default=1.0e-7)
    parser.add_argument("--atol-rel", type=float, default=1.0e-10)
    parser.add_argument(
        "--dense-output",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Keep solve_ivp dense_output enabled. Off by default to reduce memory pressure.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/07a-phi-amp-max"))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    return parser


def read_a_unit_csv(path: Path, mass_min: float | None = None, mass_max: float | None = None) -> list[float]:
    with path.open() as f:
        reader = csv.DictReader(f)
        masses = [float(row["mass_eV"]) for row in reader]
    if mass_min is not None:
        masses = [mass for mass in masses if mass >= mass_min]
    if mass_max is not None:
        masses = [mass for mass in masses if mass <= mass_max]
    return masses


def rho_dm_of_a(a: float) -> float:
    omega_dm = OMEGA_M - OMEGA_B
    return omega_dm * RHO_CRIT0 * a**-3


FIELDNAMES = [
    "mass_eV",
    "phi_unit_at_eval",
    "dphi_dx_unit_at_eval",
    "dotphi_phys_unit_at_eval",
    "rho_unit_at_eval",
    "rho_phi_max_at_eval",
    "phi_amp_max",
]


def read_existing_rows(path: Path) -> list[dict[str, float]]:
    if not path.exists():
        return []
    with path.open() as f:
        reader = csv.DictReader(f)
        return [{key: float(value) for key, value in row.items()} for row in reader]


def write_rows(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    masses = read_a_unit_csv(args.a_unit_csv, mass_min=args.mass_min, mass_max=args.mass_max)
    if not masses:
        raise RuntimeError("No masses selected. Check --mass-min/--mass-max against the input CSV.")
    x_eval = np.log(1.0 / (1.0 + args.z_eval))
    a_eval = np.exp(x_eval)
    rho_dm_eval = rho_dm_of_a(a_eval)
    rho_phi_max = args.f_dm * rho_dm_eval

    config = SolveConfig(
        z_ini=args.z_ini,
        method=args.method,
        rtol=args.rtol,
        atol_rel=args.atol_rel,
        dense_output=args.dense_output,
    )

    csv_path = args.output_dir / "phi_amp_max_scan.csv"
    rows = read_existing_rows(csv_path) if args.resume else []
    completed = {row["mass_eV"] for row in rows}

    pending = [mass for mass in masses if mass not in completed]
    if rows:
        print(f"[07a] resume: found {len(rows)} saved rows, {len(pending)} masses left")

    for idx, mass in enumerate(pending, start=1):
        sol = solve_phi_background(m_eV=mass, phi_ini=1.0, config=config)
        phi_eval, dphi_dx_eval = evaluate_state_at_x(sol, x_target=x_eval)
        H_eval = H_of_a(a_eval)
        dotphi_phys = H_eval * dphi_dx_eval
        omega = m_eV_to_omega_SI(mass)
        rho_unit = 0.5 * dotphi_phys**2 + 0.5 * omega**2 * phi_eval**2
        phi_amp_max = np.sqrt(rho_phi_max / rho_unit) if rho_unit > 0.0 else np.nan
        rows.append(
            {
                "mass_eV": float(mass),
                "phi_unit_at_eval": float(phi_eval),
                "dphi_dx_unit_at_eval": float(dphi_dx_eval),
                "dotphi_phys_unit_at_eval": float(dotphi_phys),
                "rho_unit_at_eval": float(rho_unit),
                "rho_phi_max_at_eval": float(rho_phi_max),
                "phi_amp_max": float(phi_amp_max),
            }
        )
        rows.sort(key=lambda row: row["mass_eV"])
        write_rows(csv_path, rows)
        print(
            f"[07a] {idx}/{len(pending)} pending: "
            f"m={mass:.6e} eV, rho_unit={rho_unit:.6e}, phi_amp_max={phi_amp_max:.6e}"
        )

    if not rows:
        raise RuntimeError("No rows available to save or plot.")

    mass_arr = np.array([row["mass_eV"] for row in rows])
    phi_max_arr = np.array([row["phi_amp_max"] for row in rows])

    plt.figure(figsize=(7.5, 5.2))
    plt.loglog(mass_arr, phi_max_arr)
    plt.xlabel(r"$m_a\ [{\rm eV}]$")
    plt.ylabel(r"$\phi_{\rm amp,max}(m_a)$")
    plt.title(rf"Allowed amplitude assuming $\rho_\phi(z_{{\rm eval}})\leq {args.f_dm:g}\,\rho_{{DM}}$")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output_dir / "phi_amp_max_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    payload = {
        "solve_config": asdict(config),
        "z_eval": args.z_eval,
        "a_eval": a_eval,
        "f_dm": args.f_dm,
        "rho_dm_eval": rho_dm_eval,
        "rho_phi_max": rho_phi_max,
        "assumption": "Phenomenological bound assuming rho_phi(z_eval) <= f_dm * rho_DM(z_eval) and canonical scalar energy density proxy rho_phi = 0.5 dotphi_phys^2 + 0.5 omega^2 phi^2 in the same field units as the ODE variable.",
    }
    (args.output_dir / "run_config.json").write_text(json.dumps(payload, indent=2))
    (args.output_dir / "summary.md").write_text(
        "\n".join(
            [
                "# 07a-compute_phi_amp_max summary",
                "",
                f"- input A_unit grid: `{args.a_unit_csv}`",
                f"- mass_min: `{args.mass_min}`",
                f"- mass_max: `{args.mass_max}`",
                f"- evaluation redshift: `{args.z_eval}`",
                f"- assumption: `rho_phi(z_eval) <= {args.f_dm:g} * rho_DM(z_eval)`",
                f"- solver: `{args.method}`, `rtol={args.rtol:g}`, `atol_rel={args.atol_rel:g}`, `dense_output={args.dense_output}`",
                f"- output csv: `{csv_path.name}`",
                "",
                "This is a phenomenological amplitude bound in the same field-amplitude units as the current ODE variable.",
                "For WSL2 stability, prefer splitting heavy runs by mass range and using --resume.",
            ]
        )
        + "\n"
    )

    print(f"Saved outputs to: {args.output_dir}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
