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

from patchy_reionization import Cl_phiphi_toy, Cl_tautau_toy


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-unit-csv", type=Path, default=Path("results/04a-a-unit/A_unit_scan.csv"))
    parser.add_argument("--L-max", type=int, default=1200)
    parser.add_argument("--ctau0", type=float, default=1.0e-10)
    parser.add_argument("--cphi0", type=float, default=1.0e-10)
    parser.add_argument("--Lc", type=float, default=300.0)
    parser.add_argument("--sigmaL", type=float, default=100.0)
    parser.add_argument("--Lstar", type=float, default=300.0)
    parser.add_argument("--nphi-values", type=float, nargs="+", default=[1.0, 2.0, 3.0])
    parser.add_argument("--include-cutoff", action="store_true")
    parser.add_argument("--Lcut", type=float, default=800.0)
    parser.add_argument("--targets", type=float, nargs="+", default=[0.1, 1.0, 10.0])
    parser.add_argument("--output-dir", type=Path, default=Path("results/04b-rtau-needed"))
    parser.add_argument("--no-show", action="store_true")
    return parser


def make_case_specs(args: argparse.Namespace) -> list[dict]:
    cases = []
    for nphi in args.nphi_values:
        cases.append({"name": f"nphi_{nphi:g}_no_cut", "nphi": nphi, "Lcut": None})
        if args.include_cutoff:
            cases.append({"name": f"nphi_{nphi:g}_Lcut_{args.Lcut:g}", "nphi": nphi, "Lcut": args.Lcut})
    return cases


def safe_phi_needed(target: float, rtau_max_unit: float) -> float:
    if not np.isfinite(rtau_max_unit) or rtau_max_unit <= 0.0:
        return float("nan")
    return float(np.sqrt(target / rtau_max_unit))


def read_a_unit_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    masses = []
    avals = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            masses.append(float(row["mass_eV"]))
            avals.append(float(row["A_unit"]))
    return np.array(masses), np.array(avals)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    masses, A_unit_vals = read_a_unit_csv(args.a_unit_csv)
    cases = make_case_specs(args)
    L = np.arange(2, args.L_max + 1)
    Cl_tau = Cl_tautau_toy(L, Ctau0=args.ctau0, Lc=args.Lc, sigmaL=args.sigmaL)

    rows: list[dict] = []
    best_by_case: list[dict] = []

    plt.figure(figsize=(7, 5))
    for case in cases:
        rtau_max_vals = []
        peaks = []
        for mass, A_unit in zip(masses, A_unit_vals, strict=True):
            Cl_phi = Cl_phiphi_toy(L, Cphi0=args.cphi0, Lstar=args.Lstar, nphi=case["nphi"], Lcut=case["Lcut"])
            curve = (A_unit**2) * (Cl_tau / Cl_phi)
            idx = int(np.nanargmax(np.abs(curve)))
            rtau_max = float(np.nanmax(np.abs(curve)))
            peak = int(L[idx])
            row = {
                "case_name": case["name"],
                "mass_eV": float(mass),
                "A_unit": float(A_unit),
                "Rtau_max_unit": rtau_max,
                "L_peak": peak,
                "nphi": float(case["nphi"]),
                "Lcut": "" if case["Lcut"] is None else float(case["Lcut"]),
            }
            for target in args.targets:
                row[f"phi_needed_target_{target:g}"] = safe_phi_needed(target, rtau_max)
            rows.append(row)
            rtau_max_vals.append(rtau_max)
            peaks.append(peak)

        rtau_max_vals = np.array(rtau_max_vals)
        best_idx = int(np.nanargmax(rtau_max_vals))
        best_by_case.append(
            {
                "case_name": case["name"],
                "m_best_eV": float(masses[best_idx]),
                "A_unit_at_best": float(A_unit_vals[best_idx]),
                "Rtau_max_unit": float(rtau_max_vals[best_idx]),
                "L_peak": int(peaks[best_idx]),
                "phi_needed": {f"{t:g}": safe_phi_needed(t, float(rtau_max_vals[best_idx])) for t in args.targets},
            }
        )
        plt.loglog(masses, rtau_max_vals, label=case["name"])

    plt.xlabel(r"$m_a \ [{\rm eV}]$")
    plt.ylabel(r"$R_{\tau,\max}^{\rm unit}(m_a)$")
    plt.title("Toy maximum patchy-to-genuine ratio")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "Rtau_max_unit_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    for target in args.targets:
        plt.figure(figsize=(7, 5))
        for case in cases:
            case_rows = [row for row in rows if row["case_name"] == case["name"]]
            case_m = np.array([row["mass_eV"] for row in case_rows])
            case_phi = np.array([row[f"phi_needed_target_{target:g}"] for row in case_rows])
            plt.loglog(case_m, case_phi, label=case["name"])
        plt.xlabel(r"$m_a \ [{\rm eV}]$")
        plt.ylabel(r"$\phi_{\rm needed}(m_a)$")
        plt.title(rf"Required amplitude for $R_{{\tau,\max}}={target:g}$")
        plt.grid(True, which="both", alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(args.output_dir / f"phi_needed_target_{target:g}.png", dpi=160, bbox_inches="tight")
        plt.close()

    fieldnames = [
        "case_name", "mass_eV", "A_unit", "Rtau_max_unit", "L_peak", "nphi", "Lcut",
        *[f"phi_needed_target_{t:g}" for t in args.targets],
    ]
    write_csv(args.output_dir / "Rtau_phi_needed_scan.csv", fieldnames, rows)
    (args.output_dir / "summary.md").write_text(
        "\n".join(
            ["# 04b-scan_rtau_from_aunit summary", "", "## Best mass by case", ""]
            + [
                f"- `{item['case_name']}`: `m_best={item['m_best_eV']:.6e} eV`, "
                f"`Rtau_max_unit={item['Rtau_max_unit']:.6e}`, `L_peak={item['L_peak']}`"
                for item in best_by_case
            ]
        )
        + "\n"
    )
    (args.output_dir / "run_config.json").write_text(
        json.dumps(
            {
                "a_unit_csv": str(args.a_unit_csv),
                "args": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
                "best_by_case": best_by_case,
            },
            indent=2,
        )
    )
    print(f"Saved outputs to: {args.output_dir}")
    for item in best_by_case:
        print(
            f"{item['case_name']}: m_best={item['m_best_eV']:.6e} eV, "
            f"Rtau_max_unit={item['Rtau_max_unit']:.6e}, L_peak={item['L_peak']}"
        )

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
