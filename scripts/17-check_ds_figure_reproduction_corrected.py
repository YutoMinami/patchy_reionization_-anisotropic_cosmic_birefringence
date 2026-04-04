from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from importlib import util

import matplotlib

if "--no-show" in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_module_from_path(name: str, path: Path):
    spec = util.spec_from_file_location(name, path)
    module = util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


MOD15 = load_module_from_path(
    "surrogate15",
    ROOT / "scripts" / "15-a_reff_surrogate_bound.py",
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L-min", type=int, default=2)
    parser.add_argument("--L-max", type=int, default=3000)
    parser.add_argument("--z-max", type=float, default=40.0)
    parser.add_argument("--n-z", type=int, default=5000)
    parser.add_argument("--z-min-a", type=float, default=0.1)
    parser.add_argument("--a-fid-target", type=float, default=0.07)
    parser.add_argument("--tau-fid", type=float, default=0.084)
    parser.add_argument("--delta-y-fid", type=float, default=19.0)
    parser.add_argument("--dpeak-fid-fig16", type=float, default=2.5e-5)
    parser.add_argument("--dpeak-fid-fig18", type=float, default=2.9e-5)
    parser.add_argument("--siglnl-scale", type=float, default=1.0)
    parser.add_argument("--output-dir", type=Path, default=Path("results/17-ds-figure-check-corrected"))
    parser.add_argument("--no-show", action="store_true")
    return parser


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, float | str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    L = np.arange(args.L_min, args.L_max + 1)
    z_grid = np.linspace(0.0, args.z_max, args.n_z)

    fig16_cases = [
        {"tau": 0.084, "delta_y": 19.0, "color": "black", "paper_A": 0.07, "label": r"$\tau=0.084,\ \Delta y=19.0$"},
        {"tau": 0.11, "delta_y": 14.7, "color": "red", "paper_A": 0.07, "label": r"$\tau=0.11,\ \Delta y=14.7$"},
        {"tau": 0.11, "delta_y": 19.0, "color": "royalblue", "paper_A": 0.09, "label": r"$\tau=0.11,\ \Delta y=19.0$"},
    ]
    fig18_cases = [
        {
            "Rbar": 5.0,
            "siglnR": 0.693,
            "color": "black",
            "label": r"$\bar R=5,\ \sigma_{\ln R}=0.693$",
        },
        {
            "Rbar": 6.0,
            "siglnR": 0.659,
            "color": "red",
            "label": r"$\bar R=6,\ \sigma_{\ln R}=0.659$",
        },
        {
            "Rbar": 6.0,
            "siglnR": 0.693,
            "color": "royalblue",
            "label": r"$\bar R=6,\ \sigma_{\ln R}=0.693$",
        },
    ]

    cache = {}
    for case in fig16_cases:
        cache[(case["tau"], case["delta_y"])] = MOD15.a_raw_of_tau_delta(
            tau_target=case["tau"],
            delta_y=case["delta_y"],
            z_grid=z_grid,
            z_min_a=args.z_min_a,
        )

    a_raw_fid, _ = cache[(args.tau_fid, args.delta_y_fid)]
    rows: list[dict[str, float | str]] = []

    plt.figure(figsize=(8, 5.5))
    for case in fig16_cases:
        tau = case["tau"]
        delta_y = case["delta_y"]
        a_raw, z_re = cache[(tau, delta_y)]
        a_model = args.a_fid_target * (a_raw / a_raw_fid)
        reff = 5.0 * np.exp(4.0 * 0.693**2)
        chi_re = float(np.interp(z_re, z_grid, MOD15.cumulative_trapezoid_zero(MOD15.C_KM_S / MOD15.H_of_z(z_grid), z_grid)))
        l_peak = chi_re / reff
        d_peak = args.dpeak_fid_fig16 * (a_model / args.a_fid_target)
        d_shape = d_peak * MOD15.d_shape_lognormal(L, L_peak=l_peak, sigma_lnL=args.siglnl_scale * 0.693)
        label = case["label"] + rf", $A_{{\rm model}}={a_model:.3f}$, $A_{{\rm paper}}={case['paper_A']:.2f}$"
        plt.plot(L, d_shape, lw=2.0, color=case["color"], label=label)
        rows.append(
            {
                "figure": "fig16_like_corrected",
                "tau": tau,
                "delta_y": delta_y,
                "A_model": a_model,
                "A_paper": case["paper_A"],
                "z_re": z_re,
                "Reff_Mpc": reff,
                "L_peak_surrogate": l_peak,
                "D_peak_surrogate": d_peak,
            }
        )
    plt.xlabel(r"$L$")
    plt.ylabel(r"$D_L^{\tau\tau}$ surrogate")
    plt.title("Corrected Fig.16-like check")
    plt.xlim(args.L_min, args.L_max)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.output_dir / "fig16_like_corrected.png", dpi=160, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(8, 5.5))
    tau = args.tau_fid
    delta_y = args.delta_y_fid
    a_raw, z_re = cache[(tau, delta_y)]
    a_model = args.a_fid_target * (a_raw / a_raw_fid)
    chi_re = float(np.interp(z_re, z_grid, MOD15.cumulative_trapezoid_zero(MOD15.C_KM_S / MOD15.H_of_z(z_grid), z_grid)))
    for case in fig18_cases:
        reff = case["Rbar"] * np.exp(4.0 * case["siglnR"] ** 2)
        l_peak = chi_re / reff
        d_peak = args.dpeak_fid_fig18
        d_shape = d_peak * MOD15.d_shape_lognormal(L, L_peak=l_peak, sigma_lnL=args.siglnl_scale * case["siglnR"])
        plt.plot(
            L,
            d_shape,
            lw=2.0,
            color=case["color"],
            label=case["label"] + rf", $R_{{\rm eff}}={reff:.1f}$",
        )
        rows.append(
            {
                "figure": "fig18_like_corrected",
                "tau": tau,
                "delta_y": delta_y,
                "A_model": a_model,
                "A_paper": args.a_fid_target,
                "z_re": z_re,
                "Reff_Mpc": reff,
                "L_peak_surrogate": l_peak,
                "D_peak_surrogate": d_peak,
            }
        )
    plt.xlabel(r"$L$")
    plt.ylabel(r"$D_L^{\tau\tau}$ surrogate")
    plt.title("Corrected Fig.18-like check")
    plt.xlim(args.L_min, args.L_max)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.output_dir / "fig18_like_corrected.png", dpi=160, bbox_inches="tight")
    plt.close()

    write_csv(args.output_dir / "figure_reproduction_summary.csv", list(rows[0].keys()), rows)
    summary_lines = [
        "# 17-check_ds_figure_reproduction_corrected summary",
        "",
        "Corrected Fig.16-like reproduction using the paper labels",
        "- fiducial anchor: (tau, Delta y) = (0.084, 19.0) -> A = 0.07",
        "- comparison anchors: (0.11, 14.7) -> A ~ 0.07, (0.11, 19.0) -> A ~ 0.09",
        "- Fig.18 anchors: (R, sigma_lnR) = (5, 0.693), (6, 0.659), (6, 0.693)",
        f"- adopted peak fiducials: Fig.16 -> D_peak={args.dpeak_fid_fig16:.2e}, Fig.18 -> D_peak={args.dpeak_fid_fig18:.2e}",
        "- In Fig.18-like reproduction, the peak height is held fixed and only the characteristic scale is changed.",
        "",
    ]
    for row in rows:
        if row["figure"] == "fig16_like_corrected":
            summary_lines.append(
                f"- tau={row['tau']:.3f}, Delta y={row['delta_y']:.1f}: "
                f"A_model={row['A_model']:.4f}, A_paper={row['A_paper']:.4f}, "
                f"z_re={row['z_re']:.3f}, L_peak={row['L_peak_surrogate']:.2f}"
            )
        if row["figure"] == "fig18_like_corrected":
            summary_lines.append(
                f"- Fig18 case: R_eff={row['Reff_Mpc']:.2f} Mpc, "
                f"L_peak={row['L_peak_surrogate']:.2f}, D_peak={row['D_peak_surrogate']:.3e}"
            )
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps({k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}, indent=2)
    )
    print(f"Saved outputs to: {args.output_dir}")


if __name__ == "__main__":
    main()
