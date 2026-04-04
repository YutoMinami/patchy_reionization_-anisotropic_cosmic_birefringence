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
    parser.add_argument("--A-cb-limit", type=float, default=1.0e-4)
    parser.add_argument("--L-min", type=int, default=2)
    parser.add_argument("--L-max", type=int, default=100)
    parser.add_argument("--ctau0", type=float, default=1.0e-10)
    parser.add_argument("--cphi0", type=float, default=1.0e-10)
    parser.add_argument("--Lc", type=float, default=300.0)
    parser.add_argument("--sigmaL", type=float, default=100.0)
    parser.add_argument("--Lstar", type=float, default=300.0)
    parser.add_argument("--nphi-values", type=float, nargs="+", default=[1.0, 2.0, 3.0])
    parser.add_argument("--include-cutoff", action="store_true", default=True)
    parser.add_argument("--Lcut", type=float, default=800.0)
    parser.add_argument(
        "--selected-masses",
        type=float,
        nargs="+",
        default=[2.0309176209047306e-27, 5.878016072274924e-27, 1.0e-26],
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/13-acb-reinterpretation"))
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


def make_case_specs(args: argparse.Namespace) -> list[dict[str, float | str | None]]:
    cases: list[dict[str, float | str | None]] = []
    for nphi in args.nphi_values:
        cases.append({"name": f"nphi_{nphi:g}_no_cut", "nphi": nphi, "Lcut": None})
        if args.include_cutoff:
            cases.append({"name": f"nphi_{nphi:g}_Lcut_{args.Lcut:g}", "nphi": nphi, "Lcut": args.Lcut})
    return cases


def cl_alphaalpha_limit(L: np.ndarray, A_cb_limit: float) -> np.ndarray:
    return 2.0 * np.pi * A_cb_limit / (L * (L + 1.0))


def nearest_rows(rows: list[dict[str, str]], targets: list[float]) -> list[dict[str, str]]:
    masses = np.array([float(row["mass_eV"]) for row in rows])
    selected = []
    for target in targets:
        idx = int(np.argmin(np.abs(masses - target)))
        selected.append(rows[idx])
    return selected


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matched_rows = read_csv(args.matched_csv)
    selected_rows = nearest_rows(matched_rows, args.selected_masses)
    cases = make_case_specs(args)

    L = np.arange(args.L_min, args.L_max + 1)
    claa_lim = cl_alphaalpha_limit(L, args.A_cb_limit)
    cl_tau_template = Cl_tautau_toy(L, Ctau0=args.ctau0, Lc=args.Lc, sigmaL=args.sigmaL)

    budget_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | str]] = []

    plt.figure(figsize=(8, 5.5))
    plt.loglog(L, claa_lim, color="black", lw=2, label=r"$C_L^{\alpha\alpha,\rm lim}$")
    plt.xlabel(r"$L$")
    plt.ylabel(r"$C_L^{\alpha\alpha,\rm lim}$")
    plt.title(rf"Scale-invariant limit from $A_{{\rm CB}}<{args.A_cb_limit:.1e}$")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "claa_limit.png", dpi=160, bbox_inches="tight")
    plt.close()

    for selected in selected_rows:
        mass = float(selected["mass_eV"])
        A_unit = float(selected["A_unit"])
        phi_amp_max = float(selected["phi_amp_max"])
        A_tau_eff = A_unit * phi_amp_max

        plt.figure(figsize=(8, 5.5))
        for case in cases:
            cl_phi = Cl_phiphi_toy(
                L,
                Cphi0=args.cphi0,
                Lstar=args.Lstar,
                nphi=float(case["nphi"]),
                Lcut=case["Lcut"],
            )
            patchy_power = (A_tau_eff**2) * cl_tau_template
            genuine_shape = cl_phi
            template_scale_limit = claa_lim / patchy_power
            budget_fraction = patchy_power / claa_lim

            for ell, cl_lim, patchy, frac, scale_lim in zip(
                L, claa_lim, patchy_power, budget_fraction, template_scale_limit, strict=True
            ):
                budget_rows.append(
                    {
                        "case_name": case["name"],
                        "mass_eV": mass,
                        "A_unit": A_unit,
                        "phi_amp_max": phi_amp_max,
                        "A_tau_eff": A_tau_eff,
                        "L": int(ell),
                        "Claa_limit": float(cl_lim),
                        "patchy_power_if_phi_amp_max": float(patchy),
                        "patchy_budget_fraction": float(frac),
                        "template_scale_limit": float(scale_lim),
                        "genuine_shape_toy": float(genuine_shape[ell - args.L_min]),
                    }
                )

            plt.loglog(L, template_scale_limit, label=case["name"])

            summary_rows.append(
                {
                    "case_name": case["name"],
                    "mass_eV": mass,
                    "A_tau_eff": A_tau_eff,
                    "max_patchy_budget_fraction": float(np.nanmax(budget_fraction)),
                    "L_at_max_fraction": int(L[int(np.nanargmax(budget_fraction))]),
                    "min_template_scale_limit": float(np.nanmin(template_scale_limit)),
                    "L_at_min_scale_limit": int(L[int(np.nanargmin(template_scale_limit))]),
                }
            )

        plt.axhline(1.0, color="black", ls="--", lw=1.0)
        plt.xlabel(r"$L$")
        plt.ylabel(r"allowed template rescaling")
        plt.title(rf"Allowed patchy-template rescaling at $m_a={mass:.2e}\,{{\rm eV}}$")
        plt.grid(True, which="both", alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(args.output_dir / f"template_scale_limit_{mass:.3e}.png", dpi=160, bbox_inches="tight")
        plt.close()

    write_csv(args.output_dir / "acb_reinterpretation_scan.csv", list(budget_rows[0].keys()), budget_rows)
    write_csv(args.output_dir / "acb_reinterpretation_summary.csv", list(summary_rows[0].keys()), summary_rows)

    summary_lines = [
        "# 13-acb_reinterpretation_with_matched summary",
        "",
        f"- matched csv: `{args.matched_csv}`",
        f"- A_CB limit: `{args.A_cb_limit:.3e}`",
        "",
        "Interpretation:",
        "- `patchy_budget_fraction > 1` means the toy patchy term alone would exceed the anisotropic CB budget at that L.",
        "- `template_scale_limit < 1` means the current toy tau-template normalization would need to be reduced to fit within the A_CB budget.",
        "- This is benchmark-dependent because it uses the toy C_L^{tau tau} template.",
        "",
    ]
    for item in summary_rows:
        summary_lines.append(
            f"- `{item['case_name']}`, `m={item['mass_eV']:.6e} eV`: "
            f"`max_fraction={item['max_patchy_budget_fraction']:.6e}` at `L={item['L_at_max_fraction']}`, "
            f"`min_scale_limit={item['min_template_scale_limit']:.6e}` at `L={item['L_at_min_scale_limit']}`"
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
    for item in summary_rows:
        print(
            f"{item['case_name']}, m={item['mass_eV']:.6e} eV: "
            f"max_fraction={item['max_patchy_budget_fraction']:.6e}, "
            f"min_scale_limit={item['min_template_scale_limit']:.6e}"
        )

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
