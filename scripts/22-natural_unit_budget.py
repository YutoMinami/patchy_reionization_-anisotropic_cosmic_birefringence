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
        "--natural-unit-csv",
        type=Path,
        default=Path("results/21-natural-unit-normalization/natural_unit_summary.csv"),
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
        default=Path("results/22-natural-unit-budget"),
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


def nearest_row(rows: list[dict[str, str]], target: float, g: float) -> dict[str, str]:
    candidates = [row for row in rows if abs(float(row["g_GeV_inv"]) - g) < 1e-30]
    if not candidates:
        raise ValueError(f"No row found for g={g}")
    masses = np.array([float(row["mass_eV"]) for row in candidates])
    idx = int(np.argmin(np.abs(masses - target)))
    return candidates[idx]


def claa_limit(L: np.ndarray, A_cb_limit: float) -> np.ndarray:
    return 2.0 * np.pi * A_cb_limit / (L * (L + 1.0))


def d_from_c(L: np.ndarray, C: np.ndarray) -> np.ndarray:
    return L * (L + 1.0) * C / (2.0 * np.pi)


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    nat_rows = read_csv(args.natural_unit_csv)
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
    dl_lim = d_from_c(L, cl_lim)
    cl_phi_bench = args.genuine_limit_fraction * cl_lim
    dl_phi_bench = d_from_c(L, cl_phi_bench)

    rows: list[dict[str, float | str]] = []
    delta_tau_fid = 0.084

    for g in args.g_values:
        nat_row = nearest_row(nat_rows, args.mass, g)
        mass = float(nat_row["mass_eV"])
        phi_amp_max_nat = float(nat_row["phi_amp_max_nat_GeV"])
        a_unit_nat = float(nat_row["A_unit_nat_GeV"])
        a_tau_raw_nat = float(nat_row["A_tau_raw_nat_GeV"])
        a_tau_phys = float(nat_row["A_tau_physical_dimless"])

        for delta_tau in args.delta_tau_values:
            amp_scale = (delta_tau / delta_tau_fid) ** 2
            c_tau = amp_scale * c_fid
            d_tau = amp_scale * d_fid
            cl_alpha_tau = (a_tau_phys**2) * c_tau
            dl_alpha_tau = d_from_c(L, cl_alpha_tau)
            cl_total = cl_phi_bench + cl_alpha_tau
            dl_total = d_from_c(L, cl_total)

            for ell, d_tau_val, c_tau_val, cl_tau_val, dl_tau_val, cl_total_val, dl_total_val, cl_lim_val, dl_lim_val in zip(
                L, d_tau, c_tau, cl_alpha_tau, dl_alpha_tau, cl_total, dl_total, cl_lim, dl_lim, strict=True
            ):
                rows.append(
                    {
                        "mass_eV": mass,
                        "g_GeV_inv": g,
                        "phi_amp_max_nat_GeV": phi_amp_max_nat,
                        "A_unit_nat_GeV": a_unit_nat,
                        "A_tau_raw_nat_GeV": a_tau_raw_nat,
                        "A_tau_physical": a_tau_phys,
                        "tau": args.tau,
                        "delta_y": args.delta_y,
                        "Rbar_Mpc": args.Rbar,
                        "siglnR": args.siglnR,
                        "Reff_Mpc": reff,
                        "z_re": z_re,
                        "L_peak_surrogate": l_peak,
                        "sigma_lnL": sigma_lnL,
                        "delta_tau": delta_tau,
                        "delta_tau_fid": delta_tau_fid,
                        "delta_tau_amp_scale": amp_scale,
                        "L": int(ell),
                        "D_tau_tau": float(d_tau_val),
                        "C_tau_tau": float(c_tau_val),
                        "C_alphaalpha_tau": float(cl_tau_val),
                        "D_alphaalpha_tau": float(dl_tau_val),
                        "C_alphaalpha_phi_bench": float(cl_phi_bench[ell - args.L_min]),
                        "D_alphaalpha_phi_bench": float(dl_phi_bench[ell - args.L_min]),
                        "C_alphaalpha_total": float(cl_total_val),
                        "D_alphaalpha_total": float(dl_total_val),
                        "C_alphaalpha_limit": float(cl_lim_val),
                        "D_alphaalpha_limit": float(dl_lim_val),
                    }
                )

    write_csv(args.output_dir / "natural_unit_budget_scan.csv", list(rows[0].keys()), rows)

    plt.figure(figsize=(8, 5.5))
    for delta_tau in args.delta_tau_values:
        mask = [row for row in rows if float(row["delta_tau"]) == delta_tau and abs(float(row["g_GeV_inv"]) - args.g_values[0]) < 1e-30]
        x = np.array([int(row["L"]) for row in mask])
        y = np.array([float(row["D_tau_tau"]) for row in mask])
        plt.plot(x, y, lw=2.0, label=rf"$\Delta\tau={delta_tau:g}$, $\Delta y={args.delta_y:g}$")
    plt.xlabel(r"$L$")
    plt.ylabel(r"$D_L^{\tau\tau}$")
    plt.title(rf"Patchy template family at $m_a={mass:.2e}\,{{\rm eV}}$")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "deltatau_variation_dtautau.png", dpi=160, bbox_inches="tight")
    plt.close()

    for g in args.g_values:
        plt.figure(figsize=(8, 5.5))
        for delta_tau in args.delta_tau_values:
            mask = [row for row in rows if float(row["delta_tau"]) == delta_tau and abs(float(row["g_GeV_inv"]) - g) < 1e-30]
            x = np.array([int(row["L"]) for row in mask])
            y = np.array([float(row["D_alphaalpha_tau"]) for row in mask])
            plt.plot(x, y, lw=2.0, label=rf"patchy, $\Delta\tau={delta_tau:g}$")
        plt.plot(L, dl_lim, color="black", lw=2.0, ls="--", label=r"$D_L^{\alpha\alpha,{\rm lim}}$")
        plt.xlabel(r"$L$")
        plt.ylabel(r"$D_L^{\alpha\alpha}$")
        plt.title(rf"Physical patchy birefringence with $g={g:.1e}\,{{\rm GeV}}^{{-1}}$")
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=9)
        plt.tight_layout()
        plt.savefig(args.output_dir / f"physical_patchy_dalpha_vs_limit_g_{g:.1e}.png", dpi=160, bbox_inches="tight")
        plt.close()

        plt.figure(figsize=(8, 5.5))
        plt.plot(L, dl_phi_bench, color="gray", lw=2.0, label=rf"genuine benchmark = {args.genuine_limit_fraction:.2f} limit")
        for delta_tau in args.delta_tau_values:
            mask = [row for row in rows if float(row["delta_tau"]) == delta_tau and abs(float(row["g_GeV_inv"]) - g) < 1e-30]
            x = np.array([int(row["L"]) for row in mask])
            y = np.array([float(row["D_alphaalpha_total"]) for row in mask])
            plt.plot(x, y, lw=2.0, label=rf"total, $\Delta\tau={delta_tau:g}$")
        plt.plot(L, dl_lim, color="black", lw=2.0, ls="--", label=r"$D_L^{\alpha\alpha,{\rm lim}}$")
        plt.xlabel(r"$L$")
        plt.ylabel(r"$D_L^{\alpha\alpha}$")
        plt.title(rf"Total physical birefringence budget with $g={g:.1e}\,{{\rm GeV}}^{{-1}}$")
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=9)
        plt.tight_layout()
        plt.savefig(args.output_dir / f"physical_total_dalpha_vs_limit_g_{g:.1e}.png", dpi=160, bbox_inches="tight")
        plt.close()

    mass = float(rows[0]["mass_eV"])
    summary_lines = [
        "# 22-natural_unit_budget summary",
        "",
        "This run rebuilds the birefringence budget using the natural-unit quantities from `21`.",
        "",
        f"- matched mass: `{mass:.6e} eV`",
        f"- template anchor: `tau={args.tau:.3f}`, `Delta y={args.delta_y:.1f}`, `Rbar={args.Rbar:.1f}`, `siglnR={args.siglnR:.3f}`",
        f"- surrogate: `z_re={z_re:.3f}`, `R_eff={reff:.3f} Mpc`, `L_peak={l_peak:.3f}`",
        f"- genuine benchmark: `{args.genuine_limit_fraction:.2f}` times the anisotropic-CB limit",
        "",
        "Natural-unit inputs:",
    ]
    for g in args.g_values:
        nat_row = nearest_row(nat_rows, args.mass, g)
        summary_lines.append(
            f"- `g={g:.2e} GeV^-1`: `phi_amp_max_nat={float(nat_row['phi_amp_max_nat_GeV']):.6e} GeV`, "
            f"`A_unit_nat={float(nat_row['A_unit_nat_GeV']):.6e} GeV`, "
            f"`A_tau_physical={float(nat_row['A_tau_physical_dimless']):.6e}`"
        )
    summary_lines += [
        "",
        "Interpretation:",
        "- Unlike `18` and the raw part of `19`, these curves use the natural-unit `A_tau` from `21`.",
        "- The main comparison is done in `D_L^{alpha alpha}` so the plotted amplitudes can be compared to the paper-style `D_L` language directly.",
        "- If these curves now sit near or below the anisotropic-CB limit, the previous huge amplitudes were mainly a unit mismatch.",
    ]
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps({k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}, indent=2)
    )
    print(f"Saved outputs to: {args.output_dir}")


if __name__ == "__main__":
    main()
