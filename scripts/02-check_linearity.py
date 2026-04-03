from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patchy_reionization import SolveConfig, X_REI
from patchy_reionization import compute_A_from_solution, compute_A_unit, solve_phi_background


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mass", type=float, default=3.919406774847229e-28)
    parser.add_argument("--z-ini", type=float, default=1.0e5)
    parser.add_argument("--phi-list", type=float, nargs="+", default=[1.0, 1.0e-3, 1.0e-6, 1.0e-9])
    return parser


def compute_legacy_A(m_eV: float, phi_ini: float, z_ini: float) -> float:
    config = SolveConfig(
        z_ini=z_ini,
        method="RK45",
        rtol=1.0e-8,
        atol_rel=1.0e-10,
        dense_output=False,
        n_eval=4000,
    )
    sol = solve_phi_background(m_eV=m_eV, phi_ini=phi_ini, config=config)
    return compute_A_from_solution(sol, x_target=X_REI)


def compute_robust_A(m_eV: float, phi_ini: float, z_ini: float) -> float:
    config = SolveConfig(
        z_ini=z_ini,
        method="DOP853",
        rtol=1.0e-9,
        atol_rel=1.0e-12,
        dense_output=True,
        n_eval=None,
    )
    sol = solve_phi_background(m_eV=m_eV, phi_ini=phi_ini, config=config)
    return compute_A_from_solution(sol, x_target=X_REI)


def main() -> None:
    args = build_parser().parse_args()
    robust_config = SolveConfig(
        z_ini=args.z_ini,
        method="DOP853",
        rtol=1.0e-9,
        atol_rel=1.0e-12,
        dense_output=True,
    )
    A_unit = compute_A_unit(args.mass, config=robust_config)

    header = (
        "phi_ini           legacy_A          robust_A          "
        "unit_rescaled      robust/unit-1     legacy/unit-1"
    )
    print(header)
    print("-" * len(header))

    for phi_ini in args.phi_list:
        A_legacy = compute_legacy_A(args.mass, phi_ini=phi_ini, z_ini=args.z_ini)
        A_robust = compute_robust_A(args.mass, phi_ini=phi_ini, z_ini=args.z_ini)
        A_unit_rescaled = A_unit * phi_ini
        robust_rel = (A_robust / A_unit_rescaled) - 1.0
        legacy_rel = (A_legacy / A_unit_rescaled) - 1.0
        print(
            f"{phi_ini:10.3e}  "
            f"{A_legacy:16.6e}  "
            f"{A_robust:16.6e}  "
            f"{A_unit_rescaled:16.6e}  "
            f"{robust_rel:16.6e}  "
            f"{legacy_rel:16.6e}"
        )

    print()
    print("Interpretation:")
    print("- robust_A should track unit_rescaled if the linear ODE is being solved stably.")
    print("- legacy_A reproduces the original notebook-style setup for comparison.")
    print("- If legacy_A drifts or flips sign while robust_A stays linear, the issue is numerical, not physical.")


if __name__ == "__main__":
    main()
