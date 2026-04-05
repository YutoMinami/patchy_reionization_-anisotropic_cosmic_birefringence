from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine the fine mres budget (27b) with the amplitude-scaling sensitivity "
            "(29) to estimate the best-case observable fraction achievable within the "
            "currently scanned surrogate parameter range."
        )
    )
    parser.add_argument(
        "--budget-summary",
        type=Path,
        default=Path("results/27b-natural-unit-budget-mres201/natural_unit_budget_mres_summary.csv"),
    )
    parser.add_argument(
        "--scaling-summary",
        type=Path,
        default=Path("results/29-amplitude-scaling-sensitivity/amplitude_scaling_sensitivity.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/30-bestcase-scaled-budget"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    budget = pd.read_csv(args.budget_summary)
    scaling = pd.read_csv(args.scaling_summary)

    rows = []
    for _, brow in budget.iterrows():
        g = float(brow["g_GeV_inv"])
        current_frac = float(brow["max_frac_eff"])
        for _, srow in scaling.iterrows():
            boost = float(srow["max_template_boost"])
            rows.append(
                {
                    "g_GeV_inv": g,
                    "current_max_frac_eff": current_frac,
                    "best_mass_ratio_m_over_mres": float(brow["mass_over_mres"]),
                    "best_m_eV": float(brow["mass_eV"]),
                    "best_L": float(brow["L_peak_frac_eff"]),
                    "p_amp": float(srow["amp_exponent"]),
                    "p_b": float(srow["bias_exponent"]),
                    "max_boost": boost,
                    "tau": float(srow["tau_at_max"]),
                    "delta_y": float(srow["delta_y_at_max"]),
                    "b": float(srow["b_at_max"]),
                    "A_model": float(srow["A_model_at_max"]),
                    "bestcase_frac_eff": current_frac * boost,
                }
            )

    df = pd.DataFrame(rows).sort_values(
        ["g_GeV_inv", "bestcase_frac_eff"], ascending=[True, False]
    )
    df.to_csv(outdir / "bestcase_scaled_budget.csv", index=False)

    summary_rows = []
    lines = [
        "# 30-bestcase_scaled_budget summary",
        "",
        "This run combines the fine `27b` budget with the `29` scaling sensitivity to estimate",
        "the largest anisotropic-CB budget fraction reachable within the currently scanned",
        "surrogate parameter range.",
        "",
        "It assumes that the visibility-weighted best point from `27b` can be multiplied by the",
        "maximum template boost available in each `(p_amp, p_b)` scaling scenario.",
        "",
    ]

    for g, sub in df.groupby("g_GeV_inv"):
        best = sub.iloc[0]
        summary_rows.append(best.to_dict())
        lines.extend(
            [
                f"- `g = {g:.2e} GeV^-1`",
                f"  - current best fraction: `{best['current_max_frac_eff']:.6e}`",
                f"  - best-case scaled fraction within current parameter range: `{best['bestcase_frac_eff']:.6e}`",
                f"  - achieved with `p_amp = {best['p_amp']:.0f}`, `p_b = {best['p_b']:.0f}`",
                f"  - corresponding surrogate point: `tau = {best['tau']:.3f}`, `Delta y = {best['delta_y']:.1f}`, `b = {best['b']:.1f}`, `A = {best['A_model']:.6f}`",
                "",
            ]
        )

    pd.DataFrame(summary_rows).to_csv(
        outdir / "bestcase_scaled_budget_summary.csv", index=False
    )
    (outdir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
