from __future__ import annotations

import argparse
import csv
import json
from importlib import util
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_module_from_path(name: str, path: Path):
    spec = util.spec_from_file_location(name, path)
    module = util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


MOD15 = load_module_from_path("surrogate15_for_29", ROOT / "scripts" / "15-a_reff_surrogate_bound.py")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tau-values", type=float, nargs="+", default=[0.068, 0.11])
    p.add_argument("--delta-y-values", type=float, nargs="+", default=[14.7, 19.0])
    p.add_argument("--b-values", type=float, nargs="+", default=[1.0, 2.0])
    p.add_argument("--tau-fid", type=float, default=0.068)
    p.add_argument("--delta-y-fid", type=float, default=19.0)
    p.add_argument("--b-fid", type=float, default=1.0)
    p.add_argument("--a-fid-target", type=float, default=0.07)
    p.add_argument("--z-max", type=float, default=40.0)
    p.add_argument("--n-z", type=int, default=5000)
    p.add_argument("--z-min-a", type=float, default=0.1)
    p.add_argument("--amp-exponents", type=float, nargs="+", default=[1.0, 2.0, 3.0])
    p.add_argument("--bias-exponents", type=float, nargs="+", default=[1.0, 2.0, 3.0])
    p.add_argument(
        "--required-boost-csv",
        type=Path,
        default=Path("results/28-required-template-boost/required_template_boost.csv"),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/29-amplitude-scaling-sensitivity"),
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

    z_grid = np.linspace(0.0, args.z_max, args.n_z)
    a_cache: dict[tuple[float, float], tuple[float, float]] = {}
    for tau, delta_y in product(args.tau_values, args.delta_y_values):
        a_cache[(tau, delta_y)] = MOD15.a_raw_of_tau_delta(
            tau_target=tau,
            delta_y=delta_y,
            z_grid=z_grid,
            z_min_a=args.z_min_a,
        )
    a_raw_fid, _ = a_cache[(args.tau_fid, args.delta_y_fid)]

    required_rows = read_csv(args.required_boost_csv)
    best_required_by_g_target: dict[tuple[float, float], float] = {}
    for row in required_rows:
        key = (float(row["g_GeV_inv"]), float(row["target_frac"]))
        boost = float(row["required_Dpeak_boost"])
        if key not in best_required_by_g_target or boost < best_required_by_g_target[key]:
            best_required_by_g_target[key] = boost

    rows: list[dict[str, float | str]] = []
    for p_amp, p_b in product(args.amp_exponents, args.bias_exponents):
        max_boost = 0.0
        max_case = None
        for tau, delta_y, b in product(args.tau_values, args.delta_y_values, args.b_values):
            a_raw, z_re = a_cache[(tau, delta_y)]
            a_model = args.a_fid_target * (a_raw / a_raw_fid)
            boost = (a_model / args.a_fid_target) ** p_amp * (b / args.b_fid) ** p_b
            if boost > max_boost:
                max_boost = boost
                max_case = (tau, delta_y, b, a_model, z_re)

        tau, delta_y, b, a_model, z_re = max_case
        row = {
            "amp_exponent": p_amp,
            "bias_exponent": p_b,
            "max_template_boost": max_boost,
            "tau_at_max": tau,
            "delta_y_at_max": delta_y,
            "b_at_max": b,
            "A_model_at_max": a_model,
            "z_re_at_max": z_re,
        }
        for (g, target), required in sorted(best_required_by_g_target.items()):
            row[f"required_boost_g{g:.1e}_target{target:.1e}"] = required
            row[f"meets_g{g:.1e}_target{target:.1e}"] = max_boost >= required
        rows.append(row)

    write_csv(args.output_dir / "amplitude_scaling_sensitivity.csv", list(rows[0].keys()), rows)

    lines = [
        "# 29-amplitude_scaling_sensitivity summary",
        "",
        "This run tests how much the surrogate template normalization could increase if the Dvorkin-Smith-inspired amplitude scaling were steeper than the current linear choice.",
        "",
        "The script assumes a family",
        "`D_peak ∝ (A/A_fid)^p_amp (b/b_fid)^p_b`",
        "and reports the largest boost available within the currently scanned `(tau, Delta y, b)` ranges.",
        "",
        "Best required boosts from `28`:",
    ]
    for (g, target), required in sorted(best_required_by_g_target.items()):
        lines.append(f"- `g={g:.2e} GeV^-1`, target `{target:.2e}` of the limit -> required boost `{required:.6e}`")

    lines += ["", "Scaling scenarios:"]
    for row in rows:
        lines.append(
            f"- `p_amp={row['amp_exponent']:.1f}`, `p_b={row['bias_exponent']:.1f}` -> "
            f"`max_boost={row['max_template_boost']:.6e}` "
            f"at `tau={row['tau_at_max']:.3f}`, `Delta y={row['delta_y_at_max']:.1f}`, `b={row['b_at_max']:.1f}`, `A={row['A_model_at_max']:.6f}`"
        )

    (args.output_dir / "summary.md").write_text("\n".join(lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps({k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}, indent=2)
    )
    print(f"Saved outputs to: {args.output_dir}")


if __name__ == "__main__":
    main()
