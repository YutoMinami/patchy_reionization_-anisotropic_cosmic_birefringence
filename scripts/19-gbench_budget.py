from __future__ import annotations

import argparse
import csv
import json
from importlib import util
from pathlib import Path
import sys

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


MOD15 = load_module_from_path("surrogate15", ROOT / "scripts" / "15-a_reff_surrogate_bound.py")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matched-csv",
        type=Path,
        default=Path("results/11-matched-scan-global-split/matched_scan.csv"),
    )
    parser.add_argument("--mass", type=float, default=5.878016072274924e-27)
    parser.add_argument("--L-min", type=int, default=2)
    parser.add_argument("--L-max", type=int, default=1000)
    parser.add_argument("--A-cb-limit", type=float, default=1.0e-4)
    parser.add_argument("--genuine-limit-fraction", type=float, default=0.3)
    parser.add_argument("--tau", type=float, default=0.084)
    parser.add_argument("--delta-y", type=float, default=19.0)
    parser.add_argument("--Rbar", type=float, default=5.0)
    parser.add_argument("--siglnR", type=float, default=0.693)
    parser.add_argument("--dpeak-fid", type=float, default=2.9e-5)
    parser.add_argument("--siglnl-scale", type=float, default=1.0)
    parser.add_argument("--z-max", type=float, default=40.0)
    parser.add_argument("--n-z", type=int, default=5000)
    parser.add_argument("--z-min-a", type=float, default=0.1)
    parser.add_argument("--delta-tau-values", type=float, nargs="+", default=[0.05, 0.084, 0.11])
    parser.add_argument("--g-values", type=float, nargs="+", default=[1.4e-12, 4.0e-12])
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/19-gbench-budget"),
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


def nearest_row(rows: list[dict[str, str]], target: float) -> dict[str, str]:
    masses = np.array([float(row["mass_eV"]) for row in rows])
    idx = int(np.argmin(np.abs(masses - target)))
    return rows[idx]


def claa_limit(L: np.ndarray, A_cb_limit: float) -> np.ndarray:
    return 2.0 * np.pi * A_cb_limit / (L * (L + 1.0))


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matched_row = nearest_row(read_csv(args.matched_csv), args.mass)
    mass = float(matched_row["mass_eV"])
    A_unit = float(matched_row["A_unit"])
    phi_amp_max = float(matched_row["phi_amp_max"])
    A_tau_raw = A_unit * phi_amp_max

    L = np.arange(args.L_min, args.L_max + 1)
    z_grid = np.linspace(0.0, args.z_max, args.n_z)
    _, z_re = MOD15.a_raw_of_tau_delta(args.tau, args.delta_y, z_grid, args.z_min_a)
    chi_re = float(np.interp(z_re, z_grid, MOD15.cumulative_trapezoid_zero(MOD15.C_KM_S / MOD15.H_of_z(z_grid), z_grid)))
    reff = args.Rbar * np.exp(4.0 * args.siglnR**2)
    l_peak = chi_re / reff
    sigma_lnL = args.siglnl_scale * args.siglnR
    d_fid = args.dpeak_fid * MOD15.d_shape_lognormal(L, l_peak, sigma_lnL)
    c_fid = MOD15.c_from_d(L, d_fid)
    cl_lim = claa_limit(L, args.A_cb_limit)
    cl_phi_bench = args.genuine_limit_fraction * cl_lim

    rows: list[dict[str, float | str]] = []
    delta_tau_fid = 0.084
    for delta_tau in args.delta_tau_values:
        amp_scale = (delta_tau / delta_tau_fid) ** 2
        c_tau = amp_scale * c_fid
        d_tau = amp_scale * d_fid
        cl_alpha_tau_raw = (A_tau_raw**2) * c_tau
        for g in args.g_values:
            a_tau_physical = 0.5 * g * A_tau_raw
            cl_alpha_tau_phys = (a_tau_physical**2) * c_tau
            cl_total_phys = cl_phi_bench + cl_alpha_tau_phys
            for ell, d_val, c_val, raw_val, phys_val, total_val, lim_val in zip(
                L, d_tau, c_tau, cl_alpha_tau_raw, cl_alpha_tau_phys, cl_total_phys, cl_lim, strict=True
            ):
                rows.append(
                    {
                        "mass_eV": mass,
                        "A_unit": A_unit,
                        "phi_amp_max": phi_amp_max,
                        "A_tau_raw": A_tau_raw,
                        "g_GeV_inv": g,
                        "A_tau_physical": a_tau_physical,
                        "tau": args.tau,
                        "delta_y": args.delta_y,
                        "Rbar_Mpc": args.Rbar,
                        "siglnR": args.siglnR,
                        "Reff_Mpc": reff,
                        "z_re": z_re,
                        "L_peak_surrogate": l_peak,
                        "delta_tau": delta_tau,
                        "delta_tau_amp_scale": amp_scale,
                        "L": int(ell),
                        "D_tau_tau": float(d_val),
                        "C_tau_tau": float(c_val),
                        "C_alphaalpha_tau_raw": float(raw_val),
                        "C_alphaalpha_tau_physical": float(phys_val),
                        "C_alphaalpha_phi_bench": float(cl_phi_bench[ell - args.L_min]),
                        "C_alphaalpha_total_physical": float(total_val),
                        "C_alphaalpha_limit": float(lim_val),
                    }
                )

    write_csv(args.output_dir / "gbench_budget_scan.csv", list(rows[0].keys()), rows)

    # raw plot
    plt.figure(figsize=(8, 5.5))
    for delta_tau in args.delta_tau_values:
        mask = [row for row in rows if float(row["delta_tau"]) == delta_tau and float(row["g_GeV_inv"]) == args.g_values[0]]
        x = np.array([int(row["L"]) for row in mask])
        y = np.array([float(row["C_alphaalpha_tau_raw"]) for row in mask])
        plt.plot(x, y, lw=2.0, label=rf"raw, $\Delta\tau={delta_tau:g}$")
    plt.plot(L, cl_lim, color="black", lw=2.0, ls="--", label=r"$C_L^{\alpha\alpha,\rm lim}$")
    plt.xlabel(r"$L$")
    plt.ylabel(r"$C_L^{\alpha\alpha}$")
    plt.title("Raw patchy response before inserting $g_{a\\gamma}$")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(args.output_dir / "raw_patchy_alpha_vs_limit.png", dpi=160, bbox_inches="tight")
    plt.close()

    # physical plots per g
    for g in args.g_values:
        plt.figure(figsize=(8, 5.5))
        for delta_tau in args.delta_tau_values:
            mask = [row for row in rows if float(row["delta_tau"]) == delta_tau and float(row["g_GeV_inv"]) == g]
            x = np.array([int(row["L"]) for row in mask])
            y = np.array([float(row["C_alphaalpha_tau_physical"]) for row in mask])
            plt.plot(x, y, lw=2.0, label=rf"patchy, $\Delta\tau={delta_tau:g}$")
        plt.plot(L, cl_lim, color="black", lw=2.0, ls="--", label=r"$C_L^{\alpha\alpha,\rm lim}$")
        plt.xlabel(r"$L$")
        plt.ylabel(r"$C_L^{\alpha\alpha}$")
        plt.title(rf"Physical patchy birefringence with $g={g:.1e}\,{{\rm GeV}}^{{-1}}$")
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=9)
        plt.tight_layout()
        plt.savefig(args.output_dir / f"physical_patchy_alpha_vs_limit_g_{g:.1e}.png", dpi=160, bbox_inches="tight")
        plt.close()

        plt.figure(figsize=(8, 5.5))
        plt.plot(L, cl_phi_bench, color="gray", lw=2.0, label=rf"genuine benchmark = {args.genuine_limit_fraction:.2f} limit")
        for delta_tau in args.delta_tau_values:
            mask = [row for row in rows if float(row["delta_tau"]) == delta_tau and float(row["g_GeV_inv"]) == g]
            x = np.array([int(row["L"]) for row in mask])
            y = np.array([float(row["C_alphaalpha_total_physical"]) for row in mask])
            plt.plot(x, y, lw=2.0, label=rf"total, $\Delta\tau={delta_tau:g}$")
        plt.plot(L, cl_lim, color="black", lw=2.0, ls="--", label=r"$C_L^{\alpha\alpha,\rm lim}$")
        plt.xlabel(r"$L$")
        plt.ylabel(r"$C_L^{\alpha\alpha}$")
        plt.title(rf"Total physical birefringence budget with $g={g:.1e}\,{{\rm GeV}}^{{-1}}$")
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=9)
        plt.tight_layout()
        plt.savefig(args.output_dir / f"physical_total_alpha_vs_limit_g_{g:.1e}.png", dpi=160, bbox_inches="tight")
        plt.close()

    summary_lines = [
        "# 19-gbench_budget summary",
        "",
        "This run separates the raw patchy response from the physical birefringence after inserting a benchmark ALP-photon coupling.",
        "",
        "Coupling benchmark:",
        "- Chandra / NGC1275 bound: `g_{a\\gamma\\gamma} \\lesssim 1.4 - 4.0 x 10^-12 GeV^-1` for `m_a \\lesssim 10^-12 eV`.",
        "- Source: Berg et al., arXiv:1605.01043.",
        "",
        f"- matched mass: `{mass:.6e} eV`",
        f"- raw `A_tau_raw = phi_amp_max * A_unit = {A_tau_raw:.6e}`",
        f"- template anchor: `tau={args.tau:.3f}`, `Delta y={args.delta_y:.1f}`, `Rbar={args.Rbar:.1f}`, `siglnR={args.siglnR:.3f}`",
        f"- surrogate: `z_re={z_re:.3f}`, `R_eff={reff:.3f} Mpc`, `L_peak={l_peak:.3f}`",
        "",
        "Interpretation:",
        "- `raw_patchy_alpha_vs_limit.png` is intentionally unphysical as an absolute observable because it omits the factor `g_{a\\gamma}/2`.",
        "- The `physical_*` plots insert `A_tau = (g/2) * A_tau_raw`.",
        "- This assumes the current field normalization can be interpreted canonically in GeV units; if not, an additional normalization check is needed.",
    ]
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps({k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}, indent=2)
    )
    print(f"Saved outputs to: {args.output_dir}")


if __name__ == "__main__":
    main()
