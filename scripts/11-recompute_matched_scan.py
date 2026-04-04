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
    d_eta_d_tau,
    evaluate_state_at_x,
    m_eV_to_omega_SI,
    solve_phi_background,
)


FIELDNAMES = [
    "mass_eV",
    "A_unit",
    "phi_unit_at_eval",
    "dphi_dx_unit_at_eval",
    "dotphi_phys_unit_at_eval",
    "rho_unit_at_eval",
    "rho_phi_max_at_eval",
    "phi_amp_max",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-m-min", type=float, default=1.0e-35)
    parser.add_argument("--grid-m-max", type=float, default=1.0e-26)
    parser.add_argument("--num-mass", type=int, default=40)
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    parser.add_argument("--mass-min", type=float, default=None)
    parser.add_argument("--mass-max", type=float, default=None)
    parser.add_argument("--z-ini", type=float, default=1.0e7)
    parser.add_argument("--z-eval", type=float, default=7.7)
    parser.add_argument("--f-dm", type=float, default=1.0)
    parser.add_argument("--output-dir", type=Path, default=Path("results/11-matched-scan"))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    return parser


def rho_dm_of_a(a: float) -> float:
    omega_dm = OMEGA_M - OMEGA_B
    return omega_dm * RHO_CRIT0 * a**-3


def append_csv_row(path: Path, row: dict[str, float]) -> None:
    needs_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        if needs_header:
            writer.writeheader()
        writer.writerow(row)


def read_existing_csv(path: Path) -> dict[float, dict[str, float]]:
    if not path.exists():
        return {}
    out: dict[float, dict[str, float]] = {}
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            mass = float(row["mass_eV"])
            out[mass] = {key: float(value) for key, value in row.items()}
    return out


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    config = SolveConfig(
        z_ini=args.z_ini,
        method="DOP853",
        rtol=1.0e-9,
        atol_rel=1.0e-12,
        dense_output=True,
    )
    if args.masses:
        m_grid = np.array(sorted(float(mass) for mass in args.masses))
    else:
        full_m_grid = np.logspace(np.log10(args.grid_m_min), np.log10(args.grid_m_max), args.num_mass)
        m_grid = full_m_grid
        if args.mass_min is not None:
            m_grid = m_grid[m_grid >= args.mass_min]
        if args.mass_max is not None:
            m_grid = m_grid[m_grid <= args.mass_max]
    if len(m_grid) == 0:
        raise RuntimeError("No masses selected. Check --masses or --mass-min/--mass-max against the global grid.")
    x_eval = np.log(1.0 / (1.0 + args.z_eval))
    a_eval = float(np.exp(x_eval))
    rho_dm_eval = rho_dm_of_a(a_eval)
    rho_phi_max = args.f_dm * rho_dm_eval

    csv_path = args.output_dir / "matched_scan.csv"
    if csv_path.exists() and not args.resume:
        csv_path.unlink()

    existing = read_existing_csv(csv_path) if args.resume else {}
    rows: list[dict[str, float]] = []

    for idx, mass in enumerate(m_grid, start=1):
        mass_key = float(mass)
        if mass_key in existing:
            rows.append(existing[mass_key])
            print(f"[11] {idx}/{len(m_grid)}: m={mass:.6e} eV [resume]")
            continue

        sol = solve_phi_background(m_eV=mass, phi_ini=1.0, config=config)
        phi_eval, dphi_dx_eval = evaluate_state_at_x(sol, x_target=x_eval)
        H_eval = float(H_of_a(a_eval))
        dotphi_phys = H_eval * dphi_dx_eval
        dotphi_conf = a_eval * dotphi_phys
        A_unit = dotphi_conf * d_eta_d_tau(a_eval)
        omega = m_eV_to_omega_SI(mass)
        rho_unit = 0.5 * dotphi_phys**2 + 0.5 * omega**2 * phi_eval**2
        phi_amp_max = np.sqrt(rho_phi_max / rho_unit) if rho_unit > 0.0 else np.nan

        row = {
            "mass_eV": mass_key,
            "A_unit": float(A_unit),
            "phi_unit_at_eval": float(phi_eval),
            "dphi_dx_unit_at_eval": float(dphi_dx_eval),
            "dotphi_phys_unit_at_eval": float(dotphi_phys),
            "rho_unit_at_eval": float(rho_unit),
            "rho_phi_max_at_eval": float(rho_phi_max),
            "phi_amp_max": float(phi_amp_max),
        }
        rows.append(row)
        append_csv_row(csv_path, row)
        print(
            f"[11] {idx}/{len(m_grid)}: m={mass:.6e} eV, "
            f"A_unit={A_unit:.6e}, phi_amp_max={phi_amp_max:.6e}"
        )

    rows.sort(key=lambda row: row["mass_eV"])
    mass_arr = np.array([row["mass_eV"] for row in rows])
    A_arr = np.array([row["A_unit"] for row in rows])
    phi_arr = np.array([row["phi_amp_max"] for row in rows])

    plt.figure(figsize=(7.2, 5.0))
    plt.loglog(mass_arr, np.abs(A_arr))
    plt.xlabel(r"$m_a\ [{\rm eV}]$")
    plt.ylabel(r"$|A_{\rm unit}(m_a)|$")
    plt.title("Matched high-precision A_unit")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output_dir / "A_unit_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(7.2, 5.0))
    plt.loglog(mass_arr, phi_arr)
    plt.xlabel(r"$m_a\ [{\rm eV}]$")
    plt.ylabel(r"$\phi_{\rm amp,max}(m_a)$")
    plt.title(rf"Matched high-precision $\phi_{{\rm amp,max}}$ with $f_{{\rm DM}}={args.f_dm:g}$")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output_dir / "phi_amp_max_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    payload = {
        "solve_config": asdict(config),
        "run_args": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
        "a_eval": a_eval,
        "rho_dm_eval": rho_dm_eval,
        "rho_phi_max": rho_phi_max,
    }
    (args.output_dir / "run_config.json").write_text(json.dumps(payload, indent=2))
    (args.output_dir / "summary.md").write_text(
        "\n".join(
            [
                "# 11-recompute_matched_scan summary",
                "",
                f"- global mass grid: `{args.grid_m_min:.3e}` to `{args.grid_m_max:.3e}` eV",
                f"- explicit masses: `{args.masses}`",
                f"- selected mass_min: `{args.mass_min}`",
                f"- selected mass_max: `{args.mass_max}`",
                f"- number of grid masses: `{args.num_mass}`",
                f"- number of selected masses: `{len(m_grid)}`",
                f"- z_ini: `{args.z_ini:.3e}`",
                f"- z_eval: `{args.z_eval:.3e}`",
                f"- f_dm: `{args.f_dm:g}`",
                f"- output csv: `{csv_path.name}`",
                "",
                "This scan recomputes A_unit and phi_amp_max from the same high-precision background solution.",
            ]
        )
        + "\n"
    )
    print(f"Saved outputs to: {args.output_dir}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
