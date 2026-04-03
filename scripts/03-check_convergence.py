from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patchy_reionization import SolveConfig, compute_A_unit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--masses",
        type=float,
        nargs="+",
        default=[1.0e-32, 1.0e-30, 3.919406774847229e-28],
    )
    parser.add_argument("--z-ini", type=float, default=1.0e5)
    return parser


def main() -> None:
    args = build_parser().parse_args()

    reference = SolveConfig(
        z_ini=args.z_ini,
        method="DOP853",
        rtol=1.0e-10,
        atol_rel=1.0e-13,
        dense_output=True,
    )
    test_configs = [
        ("dop853_medium", SolveConfig(z_ini=args.z_ini, method="DOP853", rtol=1.0e-9, atol_rel=1.0e-12, dense_output=True)),
        ("rk45_medium", SolveConfig(z_ini=args.z_ini, method="RK45", rtol=1.0e-8, atol_rel=1.0e-11, dense_output=True)),
        ("rk45_loose", SolveConfig(z_ini=args.z_ini, method="RK45", rtol=1.0e-7, atol_rel=1.0e-10, dense_output=True)),
    ]

    print("mass_eV           config           A_unit            rel_diff_vs_ref")
    print("-------------------------------------------------------------------")
    for mass in args.masses:
        A_ref = compute_A_unit(mass, config=reference)
        print(f"{mass:12.4e}  {'reference':14s}  {A_ref:16.6e}  {0.0:16.6e}")
        for name, config in test_configs:
            A_val = compute_A_unit(mass, config=config)
            rel_diff = (A_val / A_ref) - 1.0
            print(f"{mass:12.4e}  {name:14s}  {A_val:16.6e}  {rel_diff:16.6e}")
        print()


if __name__ == "__main__":
    main()
