from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--budget-summary-csv",
        type=Path,
        default=Path("results/27b-natural-unit-budget-mres201/natural_unit_budget_mres_summary.csv"),
    )
    p.add_argument("--targets", type=float, nargs="+", default=[1.0e-2, 1.0e-1, 1.0])
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/28-required-template-boost"),
    )
    return p


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, float | str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = read_csv(args.budget_summary_csv)
    out_rows: list[dict[str, float | str]] = []

    for row in rows:
        current_frac = float(row["max_frac_eff"])
        for target in args.targets:
            boost = target / current_frac if current_frac > 0.0 else float("inf")
            out_rows.append(
                {
                    "g_GeV_inv": float(row["g_GeV_inv"]),
                    "mass_eV": float(row["mass_eV"]),
                    "mass_over_mres": float(row["mass_over_mres"]),
                    "A_tau_eff_physical": float(row["A_tau_eff_physical"]),
                    "current_max_frac_eff": current_frac,
                    "target_frac": target,
                    "required_Dpeak_boost": boost,
                }
            )

    write_csv(args.output_dir / "required_template_boost.csv", list(out_rows[0].keys()), out_rows)

    best_by_g = {}
    for row in rows:
        g = float(row["g_GeV_inv"])
        cur = float(row["max_frac_eff"])
        if g not in best_by_g or cur > float(best_by_g[g]["max_frac_eff"]):
            best_by_g[g] = row

    lines = [
        "# 28-required_template_boost summary",
        "",
        "This run estimates how much the current surrogate `D_L^{tau tau}` normalization would need to be boosted",
        "to make the visibility-weighted patchy birefringence budget reach selected fractions of the anisotropic-CB limit.",
        "",
        "Because the budget scales linearly with the template normalization, the required boost is simply",
        "`target_frac / current_max_frac_eff`.",
        "",
    ]

    for g, row in sorted(best_by_g.items()):
        current = float(row["max_frac_eff"])
        lines.append(f"- `g = {g:.2e} GeV^-1`, best current point at `m/m_res = {float(row['mass_over_mres']):.6f}` has `max_frac_eff = {current:.6e}`")
        for target in args.targets:
            boost = target / current if current > 0.0 else float('inf')
            lines.append(f"  - to reach `{target:.2e}` of the limit, need template boost `{boost:.6e}`")
        lines.append("")

    (args.output_dir / "summary.md").write_text("\n".join(lines))
    (args.output_dir / "run_config.json").write_text(
        json.dumps({k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}, indent=2)
    )
    print(f"Saved outputs to: {args.output_dir}")


if __name__ == "__main__":
    main()
