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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-unit-csv", type=Path, default=Path("results/04a-a-unit/A_unit_scan.csv"))
    parser.add_argument("--L-min", type=int, default=2)
    parser.add_argument("--L-max", type=int, default=100)
    parser.add_argument("--A-cb-limit", type=float, default=1.0e-4)
    parser.add_argument(
        "--best-masses",
        type=float,
        nargs="+",
        default=[2.0309176209047306e-27, 5.878016072274924e-27, 1.0e-26],
    )
    parser.add_argument(
        "--A-tau-values",
        type=float,
        nargs="*",
        default=None,
        help="Optional absolute A_tau values. If provided, also compute C_L^{tau tau} limits.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/06-cltautau-bound"))
    parser.add_argument("--no-show", action="store_true")
    return parser


def read_a_unit_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    masses = []
    avals = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            masses.append(float(row["mass_eV"]))
            avals.append(float(row["A_unit"]))
    return np.array(masses), np.array(avals)


def nearest_a_unit(masses: np.ndarray, avals: np.ndarray, target_mass: float) -> tuple[float, float]:
    idx = int(np.argmin(np.abs(masses - target_mass)))
    return float(masses[idx]), float(avals[idx])


def cl_alphaalpha_limit(L: np.ndarray, A_cb_limit: float) -> np.ndarray:
    return 2.0 * np.pi * A_cb_limit / (L * (L + 1.0))


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    masses, A_unit_vals = read_a_unit_csv(args.a_unit_csv)
    L = np.arange(args.L_min, args.L_max + 1)
    claa_lim = cl_alphaalpha_limit(L, args.A_cb_limit)

    selected = [nearest_a_unit(masses, A_unit_vals, mass) for mass in args.best_masses]

    rows = []
    for mass, a_unit in selected:
        row = {"mass_eV": mass, "A_unit": a_unit}
        rows.append(row)

    with (args.output_dir / "selected_masses.json").open("w") as f:
        json.dump(
            {
                "A_cb_limit": args.A_cb_limit,
                "selected": [{"mass_eV": mass, "A_unit": a_unit} for mass, a_unit in selected],
                "note": "A_cb limit interpreted as scale-invariant birefringence amplitude with C_L^{aa} = 2pi A_CB / [L(L+1)].",
            },
            f,
            indent=2,
        )

    plt.figure(figsize=(8, 5.5))
    plt.loglog(L, claa_lim, color="black", lw=2, label=r"$C_L^{\alpha\alpha,{\rm lim}}$")
    plt.xlabel(r"$L$")
    plt.ylabel(r"$C_L^{\alpha\alpha,{\rm lim}}$")
    plt.title(rf"Anisotropic CB Limit from $A_{{CB}} < {args.A_cb_limit:.1e}$")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "Claa_limit_from_Acb.png", dpi=160, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(8, 5.5))
    for mass, a_unit in selected:
        # Patchy-only limit: A_tau^2 * C_L^{tau tau} <= C_L^{aa,lim}
        plt.loglog(L, claa_lim, label=rf"$m_a={mass:.2e}$ eV")
    plt.xlabel(r"$L$")
    plt.ylabel(r"$[A_\tau^2 C_L^{\tau\tau}]_{\rm lim}$")
    plt.title(r"Patchy-Only Bound on $A_\tau^2 C_L^{\tau\tau}$")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "Atau2Cltautau_limit.png", dpi=160, bbox_inches="tight")
    plt.close()

    limit_rows = []
    for mass, a_unit in selected:
        for ell, cl_lim in zip(L, claa_lim, strict=True):
            row = {
                "mass_eV": mass,
                "A_unit": a_unit,
                "L": int(ell),
                "Claa_limit": float(cl_lim),
                "Atau2_Cltautau_limit": float(cl_lim),
            }
            if args.A_tau_values:
                for atau in args.A_tau_values:
                    key = f"Cltautau_limit_Atau_{atau:.3e}"
                    row[key] = float(cl_lim / (atau**2))
            limit_rows.append(row)

    fieldnames = list(limit_rows[0].keys()) if limit_rows else []
    with (args.output_dir / "cltautau_bound.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(limit_rows)

    if args.A_tau_values:
        for atau in args.A_tau_values:
            plt.figure(figsize=(8, 5.5))
            for mass, a_unit in selected:
                values = []
                for ell in L:
                    cl_lim = cl_alphaalpha_limit(np.array([ell]), args.A_cb_limit)[0]
                    values.append(cl_lim / (atau**2))
                plt.loglog(L, values, label=rf"$m_a={mass:.2e}$ eV")
            plt.xlabel(r"$L$")
            plt.ylabel(r"$C_L^{\tau\tau,{\rm lim}}$")
            plt.title(rf"$C_L^{{\tau\tau}}$ limit for $|A_\tau|={atau:.2e}$")
            plt.grid(True, which="both", alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig(args.output_dir / f"Cltautau_limit_Atau_{atau:.3e}.png", dpi=160, bbox_inches="tight")
            plt.close()

    summary_lines = [
        "# 06-bound_cltautau_from_acb summary",
        "",
        "## Inputs",
        "",
        f"- A_CB upper limit: `{args.A_cb_limit:.3e}`",
        f"- multipole range: `{args.L_min}` to `{args.L_max}`",
        f"- A_unit CSV: `{args.a_unit_csv}`",
        "",
        "## Interpretation",
        "",
        "- This step does not yet assume a physical model for tau fluctuations.",
        "- It converts an observational anisotropic birefringence amplitude bound into a bound on `A_tau^2 C_L^{tau tau}`.",
        "- If absolute `A_tau` values are supplied, it also converts that into a `C_L^{tau tau}` bound.",
        "",
        "## Selected masses",
        "",
    ]
    for mass, a_unit in selected:
        summary_lines.append(f"- `m_a = {mass:.6e} eV`, `A_unit = {a_unit:.6e}`")
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")

    print(f"Saved outputs to: {args.output_dir}")
    print(f"Using A_CB upper limit = {args.A_cb_limit:.3e}")
    for mass, a_unit in selected:
        print(f"selected mass = {mass:.6e} eV, A_unit = {a_unit:.6e}")

    if not args.no-show:
        plt.show()


if __name__ == "__main__":
    main()
