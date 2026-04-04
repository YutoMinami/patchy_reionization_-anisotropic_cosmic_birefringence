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
    parser.add_argument(
        "--matched-csv",
        type=Path,
        default=Path("results/11-matched-scan-global-split/matched_scan.csv"),
    )
    parser.add_argument("--L-max", type=int, default=1200)
    parser.add_argument("--ctau0", type=float, default=1.0e-10)
    parser.add_argument("--cphi0", type=float, default=1.0e-10)
    parser.add_argument("--Lc", type=float, default=300.0)
    parser.add_argument("--sigmaL", type=float, default=100.0)
    parser.add_argument("--Lstar", type=float, default=300.0)
    parser.add_argument("--nphi-values", type=float, nargs="+", default=[1.0, 2.0, 3.0])
    parser.add_argument("--include-cutoff", action="store_true", default=True)
    parser.add_argument("--Lcut", type=float, default=800.0)
    parser.add_argument("--targets", type=float, nargs="+", default=[0.1, 1.0, 10.0])
    parser.add_argument("--output-dir", type=Path, default=Path("results/12-ratio-with-matched-scan"))
    parser.add_argument("--no-show", action="store_true")
    return parser


def make_case_specs(args: argparse.Namespace) -> list[dict[str, float | str | None]]:
    cases: list[dict[str, float | str | None]] = []
    for nphi in args.nphi_values:
        cases.append({"name": f"nphi_{nphi:g}_no_cut", "nphi": nphi, "Lcut": None})
        if args.include_cutoff:
            cases.append({"name": f"nphi_{nphi:g}_Lcut_{args.Lcut:g}", "nphi": nphi, "Lcut": args.Lcut})
    return cases


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, float | str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def safe_phi_needed(target: float, rtau_max_unit: float) -> float:
    if not np.isfinite(rtau_max_unit) or rtau_max_unit <= 0.0:
        return float("nan")
    return float(np.sqrt(target / rtau_max_unit))


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matched_rows = read_csv(args.matched_csv)
    masses = np.array([float(row["mass_eV"]) for row in matched_rows])
    a_units = np.array([float(row["A_unit"]) for row in matched_rows])
    phi_max = np.array([float(row["phi_amp_max"]) for row in matched_rows])

    cases = make_case_specs(args)
    L = np.arange(2, args.L_max + 1)
    Cl_tau = Cl_tautau_toy(L, Ctau0=args.ctau0, Lc=args.Lc, sigmaL=args.sigmaL)

    rows: list[dict[str, float | str]] = []
    best_summary: list[dict[str, float | str]] = []

    for case in cases:
        ratio_curve = []
        for mass, A_unit, phi_amp in zip(masses, a_units, phi_max, strict=True):
            Cl_phi = Cl_phiphi_toy(
                L,
                Cphi0=args.cphi0,
                Lstar=args.Lstar,
                nphi=float(case["nphi"]),
                Lcut=case["Lcut"],
            )
            spectrum_ratio_curve = Cl_tau / Cl_phi
            peak_idx = int(np.nanargmax(np.abs(spectrum_ratio_curve)))
            spectrum_ratio_at_peak = float(np.abs(spectrum_ratio_curve[peak_idx]))
            rtau_max_unit = float((A_unit**2) * spectrum_ratio_at_peak)
            phi_needed_R1 = safe_phi_needed(1.0, rtau_max_unit)
            ratio = phi_needed_R1 / phi_amp if phi_amp != 0.0 else float("nan")
            ratio_curve.append(ratio)

            row = {
                "case_name": case["name"],
                "mass_eV": float(mass),
                "A_unit": float(A_unit),
                "phi_amp_max": float(phi_amp),
                "Rtau_max_unit": rtau_max_unit,
                "L_peak": int(L[peak_idx]),
                "nphi": float(case["nphi"]),
                "Lcut": "" if case["Lcut"] is None else float(case["Lcut"]),
                "spectrum_ratio_at_peak": spectrum_ratio_at_peak,
                "phi_needed_target_1": phi_needed_R1,
                "ratio_phi_needed_over_phi_amp_max": ratio,
            }
            for target in args.targets:
                row[f"phi_needed_target_{target:g}"] = safe_phi_needed(target, rtau_max_unit)
            rows.append(row)

        ratio_curve = np.array(ratio_curve)
        best_idx = int(np.nanargmin(ratio_curve))
        worst_idx = int(np.nanargmax(ratio_curve))
        best_summary.append(
            {
                "case_name": case["name"],
                "min_ratio_mass_eV": float(masses[best_idx]),
                "min_ratio": float(ratio_curve[best_idx]),
                "max_ratio_mass_eV": float(masses[worst_idx]),
                "max_ratio": float(ratio_curve[worst_idx]),
            }
        )

    rows.sort(key=lambda row: (str(row["case_name"]), float(row["mass_eV"])))
    fieldnames = list(rows[0].keys())
    write_csv(args.output_dir / "matched_ratio_scan.csv", fieldnames, rows)

    case_names = sorted({str(row["case_name"]) for row in rows})

    plt.figure(figsize=(8, 5.5))
    plt.axhline(1.0, color="black", ls="--", lw=1.0)
    for case_name in case_names:
        case_rows = [row for row in rows if row["case_name"] == case_name]
        case_m = np.array([float(row["mass_eV"]) for row in case_rows])
        case_ratio = np.array([float(row["ratio_phi_needed_over_phi_amp_max"]) for row in case_rows])
        plt.loglog(case_m, case_ratio, label=case_name)
    plt.xlabel(r"$m_a\ [{\rm eV}]$")
    plt.ylabel(r"$\phi_{\rm needed}/\phi_{\rm amp,max}$")
    plt.title(r"Matched-scan viability ratio")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "matched_ratio_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharex=True)
    for case_name in case_names:
        case_rows = [row for row in rows if row["case_name"] == case_name]
        case_m = np.array([float(row["mass_eV"]) for row in case_rows])
        case_A = np.array([abs(float(row["A_unit"])) for row in case_rows])
        case_spec = np.array([float(row["spectrum_ratio_at_peak"]) for row in case_rows])
        axes[0].loglog(case_m, case_A, label=case_name)
        axes[1].loglog(case_m, np.sqrt(case_spec), label=case_name)

    axes[0].set_xlabel(r"$m_a\ [{\rm eV}]$")
    axes[0].set_ylabel(r"$|A_{\rm unit}(m_a)|$")
    axes[0].set_title(r"Matched $A_{\rm unit}$")
    axes[0].grid(True, which="both", alpha=0.3)

    axes[1].set_xlabel(r"$m_a\ [{\rm eV}]$")
    axes[1].set_ylabel(r"$\sqrt{\max_L(C_L^{\tau\tau}/C_L^{\phi\phi})}$")
    axes[1].set_title(r"Toy spectrum contribution")
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "matched_breakdown.png", dpi=160, bbox_inches="tight")
    plt.close()

    summary_lines = [
        "# 12-ratio_with_matched_scan summary",
        "",
        f"- matched csv: `{args.matched_csv}`",
        f"- targets: `{args.targets}`",
        "",
        "## Min and max ratio by case",
        "",
    ]
    for item in best_summary:
        summary_lines.append(
            f"- `{item['case_name']}`: "
            f"`min_ratio={item['min_ratio']:.6e}` at `m={item['min_ratio_mass_eV']:.6e} eV`, "
            f"`max_ratio={item['max_ratio']:.6e}` at `m={item['max_ratio_mass_eV']:.6e} eV`"
        )
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps(
            {
                "matched_csv": str(args.matched_csv),
                "args": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
                "best_summary": best_summary,
            },
            indent=2,
        )
    )

    print(f"Saved outputs to: {args.output_dir}")
    for item in best_summary:
        print(
            f"{item['case_name']}: min_ratio={item['min_ratio']:.6e} at m={item['min_ratio_mass_eV']:.6e} eV, "
            f"max_ratio={item['max_ratio']:.6e} at m={item['max_ratio_mass_eV']:.6e} eV"
        )

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
