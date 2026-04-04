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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rtau-csv",
        type=Path,
        default=Path("results/04b-rtau-needed/Rtau_phi_needed_scan.csv"),
    )
    parser.add_argument(
        "--phi-max-csv",
        type=Path,
        default=Path("results/07a-phi-amp-max-split/phi_amp_max_scan.csv"),
    )
    parser.add_argument("--target", type=float, default=1.0)
    parser.add_argument("--output-dir", type=Path, default=Path("results/08-ratio-sanity-check"))
    parser.add_argument("--no-show", action="store_true")
    return parser


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def safe_float(value: str) -> float:
    if value == "":
        return float("nan")
    return float(value)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, float | str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rtau_rows = read_csv(args.rtau_csv)
    phi_max_rows = read_csv(args.phi_max_csv)
    phi_needed_key = f"phi_needed_target_{args.target:g}"

    phi_max_by_mass = {float(row["mass_eV"]): float(row["phi_amp_max"]) for row in phi_max_rows}

    rows: list[dict[str, float | str]] = []
    for row in rtau_rows:
        mass = float(row["mass_eV"])
        phi_amp_max = phi_max_by_mass.get(mass)
        if phi_amp_max is None:
            continue

        phi_needed = float(row[phi_needed_key])
        A_unit = float(row["A_unit"])
        rtau_max_unit = float(row["Rtau_max_unit"])
        spectrum_ratio_at_peak = rtau_max_unit / (A_unit**2) if A_unit != 0.0 else float("nan")
        ratio = phi_needed / phi_amp_max if phi_amp_max != 0.0 else float("nan")

        rows.append(
            {
                "case_name": row["case_name"],
                "mass_eV": mass,
                "A_unit": A_unit,
                "A_unit_abs": abs(A_unit),
                "Rtau_max_unit": rtau_max_unit,
                "sqrt_Rtau_max_unit": float(np.sqrt(rtau_max_unit)),
                "spectrum_ratio_at_peak": spectrum_ratio_at_peak,
                "sqrt_spectrum_ratio_at_peak": float(np.sqrt(spectrum_ratio_at_peak))
                if np.isfinite(spectrum_ratio_at_peak) and spectrum_ratio_at_peak >= 0.0
                else float("nan"),
                "phi_needed": phi_needed,
                "phi_amp_max": phi_amp_max,
                "ratio_phi_needed_over_phi_amp_max": ratio,
                "log10_ratio": float(np.log10(ratio)) if np.isfinite(ratio) and ratio > 0.0 else float("nan"),
                "L_peak": int(float(row["L_peak"])),
                "nphi": float(row["nphi"]),
                "Lcut": safe_float(row["Lcut"]),
            }
        )

    if not rows:
        raise RuntimeError("No overlapping masses found between rtau and phi_amp_max CSV inputs.")

    rows.sort(key=lambda item: (str(item["case_name"]), float(item["mass_eV"])))
    fieldnames = list(rows[0].keys())
    write_csv(args.output_dir / "ratio_sanity_check.csv", fieldnames, rows)

    case_names = sorted({str(row["case_name"]) for row in rows})
    best_summary: list[dict[str, float | str]] = []

    plt.figure(figsize=(8, 5.5))
    plt.axhline(1.0, color="black", ls="--", lw=1.0)
    for case_name in case_names:
        case_rows = [row for row in rows if row["case_name"] == case_name]
        masses = np.array([float(row["mass_eV"]) for row in case_rows])
        ratios = np.array([float(row["ratio_phi_needed_over_phi_amp_max"]) for row in case_rows])
        plt.loglog(masses, ratios, label=case_name)

        best_row = min(case_rows, key=lambda item: float(item["ratio_phi_needed_over_phi_amp_max"]))
        worst_row = max(case_rows, key=lambda item: float(item["ratio_phi_needed_over_phi_amp_max"]))
        best_summary.append(
            {
                "case_name": case_name,
                "min_ratio_mass_eV": float(best_row["mass_eV"]),
                "min_ratio": float(best_row["ratio_phi_needed_over_phi_amp_max"]),
                "max_ratio_mass_eV": float(worst_row["mass_eV"]),
                "max_ratio": float(worst_row["ratio_phi_needed_over_phi_amp_max"]),
            }
        )

    plt.xlabel(r"$m_a\ [{\rm eV}]$")
    plt.ylabel(r"$\phi_{\rm needed}/\phi_{\rm amp,max}$")
    plt.title(rf"Sanity check ratio for $R_{{\tau,\max}}={args.target:g}$")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "ratio_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharex=True)
    for case_name in case_names:
        case_rows = [row for row in rows if row["case_name"] == case_name]
        masses = np.array([float(row["mass_eV"]) for row in case_rows])
        A_abs = np.array([float(row["A_unit_abs"]) for row in case_rows])
        sqrt_spec = np.array([float(row["sqrt_spectrum_ratio_at_peak"]) for row in case_rows])
        axes[0].loglog(masses, A_abs, label=case_name)
        axes[1].loglog(masses, sqrt_spec, label=case_name)

    axes[0].set_xlabel(r"$m_a\ [{\rm eV}]$")
    axes[0].set_ylabel(r"$|A_{\rm unit}(m_a)|$")
    axes[0].set_title(r"Response amplitude")
    axes[0].grid(True, which="both", alpha=0.3)

    axes[1].set_xlabel(r"$m_a\ [{\rm eV}]$")
    axes[1].set_ylabel(r"$\sqrt{\max_L(C_L^{\tau\tau}/C_L^{\phi\phi})}$")
    axes[1].set_title(r"Toy spectrum contribution at peak")
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(args.output_dir / "aunit_and_spectrum_breakdown.png", dpi=160, bbox_inches="tight")
    plt.close()

    summary_lines = [
        "# 08-sanity_check_ratio summary",
        "",
        f"- target: `{args.target:g}`",
        f"- rtau csv: `{args.rtau_csv}`",
        f"- phi_max csv: `{args.phi_max_csv}`",
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
                "target": args.target,
                "rtau_csv": str(args.rtau_csv),
                "phi_max_csv": str(args.phi_max_csv),
                "num_rows": len(rows),
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
