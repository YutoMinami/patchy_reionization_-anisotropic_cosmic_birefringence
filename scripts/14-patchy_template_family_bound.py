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
        "--matched-csv",
        type=Path,
        default=Path("results/11-matched-scan-global-split/matched_scan.csv"),
    )
    parser.add_argument("--A-cb-limit", type=float, default=1.0e-4)
    parser.add_argument("--L-min", type=int, default=2)
    parser.add_argument("--L-max", type=int, default=100)
    parser.add_argument(
        "--selected-masses",
        type=float,
        nargs="+",
        default=[2.0309176209047306e-27, 5.878016072274924e-27, 1.0e-26],
    )
    parser.add_argument(
        "--Lpeak-values",
        type=float,
        nargs="+",
        default=[20.0, 50.0, 100.0, 300.0, 600.0],
    )
    parser.add_argument(
        "--sigma-log-values",
        type=float,
        nargs="+",
        default=[0.3, 0.5, 0.8],
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/14-patchy-template-family"),
    )
    parser.add_argument("--no-show", action="store_true")
    return parser


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, float | str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def nearest_rows(rows: list[dict[str, str]], targets: list[float]) -> list[dict[str, str]]:
    masses = np.array([float(row["mass_eV"]) for row in rows])
    selected = []
    for target in targets:
        idx = int(np.argmin(np.abs(masses - target)))
        selected.append(rows[idx])
    return selected


def d_patchy_lognormal(L: np.ndarray, L_peak: float, sigma_log: float) -> np.ndarray:
    shape = np.exp(-0.5 * (np.log(L / L_peak) / sigma_log) ** 2)
    return shape


def c_from_d(L: np.ndarray, D: np.ndarray) -> np.ndarray:
    return 2.0 * np.pi * D / (L * (L + 1.0))


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matched_rows = read_csv(args.matched_csv)
    selected_rows = nearest_rows(matched_rows, args.selected_masses)

    L = np.arange(args.L_min, args.L_max + 1)
    template_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | str]] = []

    plt.figure(figsize=(8, 5.5))
    for L_peak in args.Lpeak_values:
        for sigma_log in args.sigma_log_values:
            D_shape = d_patchy_lognormal(L, L_peak=L_peak, sigma_log=sigma_log)
            label = rf"$L_{{\rm peak}}={L_peak:g},\,\sigma_{{\ln L}}={sigma_log:g}$"
            plt.semilogx(L, D_shape, lw=1.8, label=label)
            for ell, d_val in zip(L, D_shape, strict=True):
                template_rows.append(
                    {
                        "template_name": f"Lpeak_{L_peak:g}_siglog_{sigma_log:g}",
                        "L_peak": L_peak,
                        "sigma_log": sigma_log,
                        "L": int(ell),
                        "D_shape_unit_peak": float(d_val),
                        "C_shape_unit_peak": float(c_from_d(np.array([ell]), np.array([d_val]))[0]),
                    }
                )
    plt.xlabel(r"$L$")
    plt.ylabel(r"unit-peak $D_L^{\tau\tau}$ shape")
    plt.title("Literature-inspired patchy template family")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(args.output_dir / "template_family_shapes.png", dpi=160, bbox_inches="tight")
    plt.close()

    for selected in selected_rows:
        mass = float(selected["mass_eV"])
        A_unit = float(selected["A_unit"])
        phi_amp_max = float(selected["phi_amp_max"])
        A_tau_eff = A_unit * phi_amp_max

        xs = []
        ys = []
        labels = []

        for L_peak in args.Lpeak_values:
            for sigma_log in args.sigma_log_values:
                D_shape = d_patchy_lognormal(L, L_peak=L_peak, sigma_log=sigma_log)
                idx_peak_in_window = int(np.nanargmax(D_shape))
                L_at_window_peak = int(L[idx_peak_in_window])
                max_shape_in_window = float(D_shape[idx_peak_in_window])
                D_peak_limit = args.A_cb_limit / (A_tau_eff**2 * max_shape_in_window)
                C_peak_limit_at_window = float(
                    c_from_d(np.array([L_at_window_peak]), np.array([D_peak_limit]))[0]
                )

                summary_rows.append(
                    {
                        "template_name": f"Lpeak_{L_peak:g}_siglog_{sigma_log:g}",
                        "mass_eV": mass,
                        "A_unit": A_unit,
                        "phi_amp_max": phi_amp_max,
                        "A_tau_eff": A_tau_eff,
                        "L_peak_template": L_peak,
                        "sigma_log": sigma_log,
                        "L_at_window_peak": L_at_window_peak,
                        "max_shape_in_window": max_shape_in_window,
                        "D_peak_limit": D_peak_limit,
                        "C_limit_at_window_peak": C_peak_limit_at_window,
                    }
                )
                xs.append(L_peak)
                ys.append(D_peak_limit)
                labels.append(f"{sigma_log:g}")

        plt.figure(figsize=(8, 5.5))
        for sigma_log in args.sigma_log_values:
            mask = [row for row in summary_rows if row["mass_eV"] == mass and row["sigma_log"] == sigma_log]
            xvals = [float(row["L_peak_template"]) for row in mask]
            yvals = [float(row["D_peak_limit"]) for row in mask]
            plt.loglog(xvals, yvals, marker="o", lw=1.8, label=rf"$\sigma_{{\ln L}}={sigma_log:g}$")
        plt.xlabel(r"template $L_{\rm peak}$")
        plt.ylabel(r"allowed peak $D_L^{\tau\tau}$")
        plt.title(rf"Allowed unit-peak template normalization at $m_a={mass:.2e}\,{{\rm eV}}$")
        plt.grid(True, which="both", alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(args.output_dir / f"dpeak_limit_vs_lpeak_{mass:.3e}.png", dpi=160, bbox_inches="tight")
        plt.close()

    write_csv(args.output_dir / "template_family_shapes.csv", list(template_rows[0].keys()), template_rows)
    write_csv(args.output_dir / "template_family_bound_summary.csv", list(summary_rows[0].keys()), summary_rows)

    summary_lines = [
        "# 14-patchy_template_family_bound summary",
        "",
        "This run replaces the old toy Gaussian C_L^{tau tau} with a lightweight literature-inspired family.",
        "The family is a unit-peak lognormal bump in D_L^{tau tau} = L(L+1)C_L^{tau tau}/2pi.",
        "",
        "Interpretation:",
        "- `D_peak_limit` is the largest allowed peak amplitude of the template in D_L^{tau tau}.",
        "- The dependence on template family enters through how much of the bump falls inside the constrained multipole window.",
        "- High-L-peaked templates are less constrained in this simplified reinterpretation because the published anisotropic-CB amplitude mainly constrains low to intermediate L.",
        "",
    ]
    for mass in sorted({float(row["mass_eV"]) for row in summary_rows}):
        subset = [row for row in summary_rows if float(row["mass_eV"]) == mass]
        best = min(subset, key=lambda row: float(row["D_peak_limit"]))
        loosest = max(subset, key=lambda row: float(row["D_peak_limit"]))
        summary_lines.append(
            f"- `m={mass:.6e} eV`: strongest limit `D_peak<{float(best['D_peak_limit']):.6e}` "
            f"for template `(L_peak={float(best['L_peak_template']):g}, sigma_log={float(best['sigma_log']):g})`; "
            f"weakest limit `D_peak<{float(loosest['D_peak_limit']):.6e}` "
            f"for template `(L_peak={float(loosest['L_peak_template']):g}, sigma_log={float(loosest['sigma_log']):g})`."
        )
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps(
            {
                "matched_csv": str(args.matched_csv),
                "A_cb_limit": args.A_cb_limit,
                "args": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
                "selected_rows": selected_rows,
            },
            indent=2,
        )
    )

    print(f"Saved outputs to: {args.output_dir}")
    for mass in sorted({float(row['mass_eV']) for row in summary_rows}):
        subset = [row for row in summary_rows if float(row["mass_eV"]) == mass]
        best = min(subset, key=lambda row: float(row["D_peak_limit"]))
        loosest = max(subset, key=lambda row: float(row["D_peak_limit"]))
        print(
            f"m={mass:.6e} eV: strongest D_peak_limit={float(best['D_peak_limit']):.6e} "
            f"at L_peak={float(best['L_peak_template']):g}, sigma_log={float(best['sigma_log']):g}; "
            f"weakest D_peak_limit={float(loosest['D_peak_limit']):.6e} "
            f"at L_peak={float(loosest['L_peak_template']):g}, sigma_log={float(loosest['sigma_log']):g}"
        )

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
