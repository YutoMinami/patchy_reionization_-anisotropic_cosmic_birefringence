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
        "--normalization-csv",
        type=Path,
        default=Path("results/09-normalization-check/normalization_check.csv"),
    )
    parser.add_argument("--high-mass-threshold", type=float, default=1.0e-31)
    parser.add_argument("--output-dir", type=Path, default=Path("results/10-constant-scaling-check"))
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


def fit_constant(values: np.ndarray) -> tuple[float, float]:
    mean = float(np.mean(values))
    rms_frac = float(np.sqrt(np.mean(((values / mean) - 1.0) ** 2)))
    return mean, rms_frac


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    raw_rows = read_csv(args.normalization_csv)
    rows: list[dict[str, float]] = []
    for row in raw_rows:
        mass = float(row["mass_eV"])
        saved = float(row["A_unit_saved"])
        recon = float(row["A_unit_reconstructed"])
        ratio = recon / saved if saved != 0.0 else float("nan")
        rows.append(
            {
                "mass_eV": mass,
                "A_unit_saved": saved,
                "A_unit_reconstructed": recon,
                "scaling_ratio": ratio,
                "high_mass_flag": 1.0 if mass >= args.high_mass_threshold else 0.0,
            }
        )

    masses = np.array([row["mass_eV"] for row in rows])
    ratios = np.array([row["scaling_ratio"] for row in rows])
    high_mask = masses >= args.high_mass_threshold
    low_mask = ~high_mask

    all_mean, all_rms = fit_constant(ratios)
    high_mean, high_rms = fit_constant(ratios[high_mask])
    low_mean, low_rms = fit_constant(ratios[low_mask])

    for row in rows:
        ratio = row["scaling_ratio"]
        row["residual_all"] = ratio / all_mean - 1.0
        row["residual_high"] = ratio / high_mean - 1.0 if row["high_mass_flag"] else float("nan")
        row["residual_low"] = ratio / low_mean - 1.0 if not row["high_mass_flag"] else float("nan")

    write_csv(args.output_dir / "constant_scaling_check.csv", list(rows[0].keys()), rows)

    plt.figure(figsize=(7.5, 5.2))
    plt.semilogx(masses, ratios, marker="o", ms=3, lw=1.2, label="A_reconstructed / A_saved")
    plt.axhline(all_mean, color="black", ls="--", lw=1.0, label=f"all mean = {all_mean:.3f}")
    plt.axhline(high_mean, color="tab:red", ls=":", lw=1.2, label=f"high-mass mean = {high_mean:.3f}")
    plt.axvline(args.high_mass_threshold, color="gray", ls="--", lw=1.0)
    plt.xlabel(r"$m_a\ [{\rm eV}]$")
    plt.ylabel(r"$A_{\rm reconstructed}/A_{\rm saved}$")
    plt.title("Constant scaling check")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "scaling_ratio_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(7.5, 5.2))
    plt.semilogx(masses, np.abs([row["residual_all"] for row in rows]), label="all-mass constant")
    plt.semilogx(
        masses[high_mask],
        np.abs([row["residual_high"] for row in rows if row["high_mass_flag"]]),
        label="high-mass constant",
    )
    plt.axvline(args.high_mass_threshold, color="gray", ls="--", lw=1.0)
    plt.xlabel(r"$m_a\ [{\rm eV}]$")
    plt.ylabel(r"fractional residual")
    plt.title("Residual after constant rescaling")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "constant_residual_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    summary_lines = [
        "# 10-check_constant_scaling summary",
        "",
        f"- input normalization csv: `{args.normalization_csv}`",
        f"- high-mass threshold: `{args.high_mass_threshold:.6e} eV`",
        "",
        f"- all-mass mean scaling ratio: `{all_mean:.6e}` with rms fractional scatter `{all_rms:.6e}`",
        f"- low-mass mean scaling ratio: `{low_mean:.6e}` with rms fractional scatter `{low_rms:.6e}`",
        f"- high-mass mean scaling ratio: `{high_mean:.6e}` with rms fractional scatter `{high_rms:.6e}`",
        "",
        "Interpretation:",
        "- If the rms scatter is small, a constant rescaling is a good approximation in that mass range.",
        "- If the rms scatter is large, the mismatch is mass-dependent and should not be fixed by one global factor.",
    ]
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps(
            {
                "normalization_csv": str(args.normalization_csv),
                "high_mass_threshold": args.high_mass_threshold,
                "all_mean": all_mean,
                "all_rms_fractional_scatter": all_rms,
                "low_mean": low_mean,
                "low_rms_fractional_scatter": low_rms,
                "high_mean": high_mean,
                "high_rms_fractional_scatter": high_rms,
            },
            indent=2,
        )
    )

    print(f"Saved outputs to: {args.output_dir}")
    print(f"all-mass mean = {all_mean:.6e}, rms scatter = {all_rms:.6e}")
    print(f"low-mass mean = {low_mean:.6e}, rms scatter = {low_rms:.6e}")
    print(f"high-mass mean = {high_mean:.6e}, rms scatter = {high_rms:.6e}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
