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

from patchy_reionization import Cl_phiphi_toy, Cl_tautau_toy
from patchy_reionization import SolveConfig, compute_A_unit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m-min", type=float, default=1.0e-35)
    parser.add_argument("--m-max", type=float, default=1.0e-26)
    parser.add_argument("--num-mass", type=int, default=70)
    parser.add_argument("--z-ini", type=float, default=1.0e7)
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
    parser.add_argument("--output-dir", type=Path, default=Path("results/04-rtau_phi_needed"))
    parser.add_argument("--no-show", action="store_true")
    return parser


def make_case_specs(args: argparse.Namespace) -> list[dict]:
    cases = []
    for nphi in args.nphi_values:
        cases.append(
            {
                "name": f"nphi_{nphi:g}_no_cut",
                "nphi": nphi,
                "Lcut": None,
            }
        )
        if args.include_cutoff:
            cases.append(
                {
                    "name": f"nphi_{nphi:g}_Lcut_{args.Lcut:g}",
                    "nphi": nphi,
                    "Lcut": args.Lcut,
                }
            )
    return cases


def save_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def append_csv_row(path: Path, fieldnames: list[str], row: dict) -> None:
    needs_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if needs_header:
            writer.writeheader()
        writer.writerow(row)


def write_summary_md(
    path: Path,
    *,
    args: argparse.Namespace,
    cases: list[dict],
    best_by_case: list[dict],
) -> None:
    lines = [
        "# 04-scan_rtau_needed summary",
        "",
        "## Run configuration",
        "",
        f"- mass range: `{args.m_min:.3e}` to `{args.m_max:.3e}` eV",
        f"- number of masses: `{args.num_mass}`",
        f"- z_ini: `{args.z_ini:.3e}`",
        f"- L range: `2` to `{args.L_max}`",
        f"- targets: `{', '.join(f'{t:g}' for t in args.targets)}`",
        "",
        "## Toy spectrum cases",
        "",
    ]
    for case in cases:
        lcut_label = "None" if case["Lcut"] is None else f"{case['Lcut']:g}"
        lines.append(f"- `{case['name']}`: `nphi={case['nphi']:g}`, `Lcut={lcut_label}`")

    lines.extend(["", "## Best mass by case", ""])
    for item in best_by_case:
        lines.append(
            f"- `{item['case_name']}`: "
            f"`m_best={item['m_best_eV']:.6e} eV`, "
            f"`A_unit={item['A_unit_at_best']:.6e}`, "
            f"`Rtau_max_unit={item['Rtau_max_unit']:.6e}`, "
            f"`L_peak={item['L_peak']}`"
        )
        for target_key, value in item["phi_needed"].items():
            lines.append(f"  target `{target_key}` -> `phi_needed={value:.6e}`")

    path.write_text("\n".join(lines) + "\n")


def safe_phi_needed(target: float, rtau_max_unit: float) -> float:
    if not np.isfinite(rtau_max_unit) or rtau_max_unit <= 0.0:
        return float("nan")
    return float(np.sqrt(target / rtau_max_unit))


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    config = SolveConfig(z_ini=args.z_ini, method="DOP853", rtol=1.0e-9, atol_rel=1.0e-12, dense_output=True)
    cases = make_case_specs(args)

    m_grid = np.logspace(np.log10(args.m_min), np.log10(args.m_max), args.num_mass)
    L = np.arange(2, args.L_max + 1)
    Cl_tau = Cl_tautau_toy(L, Ctau0=args.ctau0, Lc=args.Lc, sigmaL=args.sigmaL)
    a_unit_csv = args.output_dir / "A_unit_scan.csv"
    if a_unit_csv.exists():
        a_unit_csv.unlink()
    a_rows = []
    A_unit_vals_list = []
    total_mass = len(m_grid)
    for index, mass in enumerate(m_grid, start=1):
        A_unit = compute_A_unit(m_eV=mass, config=config)
        A_unit_vals_list.append(A_unit)
        row = {"mass_eV": float(mass), "A_unit": float(A_unit)}
        a_rows.append(row)
        append_csv_row(a_unit_csv, ["mass_eV", "A_unit"], row)
        print(f"[A_unit] {index}/{total_mass}: m={mass:.6e} eV, A_unit={A_unit:.6e}")

    A_unit_vals = np.array(A_unit_vals_list)

    all_rows: list[dict] = []
    best_by_case: list[dict] = []

    plt.figure(figsize=(7, 5))
    plt.loglog(m_grid, np.abs(A_unit_vals))
    plt.xlabel(r"$m_a \ [{\rm eV}]$")
    plt.ylabel(r"$|A_{\rm unit}(m_a)|$")
    plt.title("Unit-amplitude response scan")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output_dir / "A_unit_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(7, 5))
    scan_csv = args.output_dir / "Rtau_phi_needed_scan.csv"
    if scan_csv.exists():
        scan_csv.unlink()
    fieldnames = [
        "case_name",
        "mass_eV",
        "A_unit",
        "Rtau_max_unit",
        "L_peak",
        "nphi",
        "Lcut",
        *[f"phi_needed_target_{target:g}" for target in args.targets],
    ]
    for case in cases:
        rtau_max_vals = []
        peak_ls = []
        print(f"[case] {case['name']} start")

        for index, (mass, A_unit) in enumerate(zip(m_grid, A_unit_vals, strict=True), start=1):
            Cl_phi = Cl_phiphi_toy(
                L,
                Cphi0=args.cphi0,
                Lstar=args.Lstar,
                nphi=case["nphi"],
                Lcut=case["Lcut"],
            )
            rtau_curve = (A_unit**2) * (Cl_tau / Cl_phi)
            idx = int(np.nanargmax(np.abs(rtau_curve)))
            rtau_max_unit = float(np.nanmax(np.abs(rtau_curve)))
            L_peak = int(L[idx])

            row = {
                "case_name": case["name"],
                "mass_eV": float(mass),
                "A_unit": float(A_unit),
                "Rtau_max_unit": rtau_max_unit,
                "L_peak": L_peak,
                "nphi": float(case["nphi"]),
                "Lcut": "" if case["Lcut"] is None else float(case["Lcut"]),
            }
            for target in args.targets:
                key = f"phi_needed_target_{target:g}"
                row[key] = safe_phi_needed(target, rtau_max_unit)

            all_rows.append(row)
            append_csv_row(scan_csv, fieldnames, row)
            rtau_max_vals.append(rtau_max_unit)
            peak_ls.append(L_peak)
            print(
                f"[case] {case['name']} {index}/{total_mass}: "
                f"m={mass:.6e} eV, Rtau_max_unit={rtau_max_unit:.6e}, L_peak={L_peak}"
            )

        rtau_max_vals_arr = np.array(rtau_max_vals)
        peak_ls_arr = np.array(peak_ls)
        best_idx = int(np.nanargmax(rtau_max_vals_arr))
        best_item = {
            "case_name": case["name"],
            "m_best_eV": float(m_grid[best_idx]),
            "A_unit_at_best": float(A_unit_vals[best_idx]),
            "Rtau_max_unit": float(rtau_max_vals_arr[best_idx]),
            "L_peak": int(peak_ls_arr[best_idx]),
            "phi_needed": {
                f"{target:g}": safe_phi_needed(target, float(rtau_max_vals_arr[best_idx]))
                for target in args.targets
            },
        }
        best_by_case.append(best_item)

        plt.loglog(m_grid, rtau_max_vals_arr, label=case["name"])
        write_summary_md(args.output_dir / "summary.md", args=args, cases=cases, best_by_case=best_by_case)
        print(
            f"[case] {case['name']} done: "
            f"m_best={best_item['m_best_eV']:.6e} eV, "
            f"Rtau_max_unit={best_item['Rtau_max_unit']:.6e}"
        )

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
            case_rows = [row for row in all_rows if row["case_name"] == case["name"]]
            masses = np.array([row["mass_eV"] for row in case_rows])
            phi_needed = np.array([row[f"phi_needed_target_{target:g}"] for row in case_rows])
            plt.loglog(masses, phi_needed, label=case["name"])
        plt.xlabel(r"$m_a \ [{\rm eV}]$")
        plt.ylabel(r"$\phi_{\rm needed}(m_a)$")
        plt.title(rf"Required amplitude for $R_{{\tau,\max}}={target:g}$")
        plt.grid(True, which="both", alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(args.output_dir / f"phi_needed_target_{target:g}.png", dpi=160, bbox_inches="tight")
        plt.close()

    config_payload = {
        "solve_config": asdict(config),
        "run_args": {
            key: (str(value) if isinstance(value, Path) else value)
            for key, value in vars(args).items()
        },
        "cases": cases,
        "best_by_case": best_by_case,
    }
    (args.output_dir / "run_config.json").write_text(json.dumps(config_payload, indent=2))
    write_summary_md(args.output_dir / "summary.md", args=args, cases=cases, best_by_case=best_by_case)

    print(f"Saved outputs to: {args.output_dir}")
    for item in best_by_case:
        print(
            f"{item['case_name']}: "
            f"m_best={item['m_best_eV']:.6e} eV, "
            f"Rtau_max_unit={item['Rtau_max_unit']:.6e}, "
            f"L_peak={item['L_peak']}"
        )

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
