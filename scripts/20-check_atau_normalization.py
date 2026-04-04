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
    parser.add_argument("--tau", type=float, default=0.084)
    parser.add_argument("--delta-y", type=float, default=19.0)
    parser.add_argument("--Rbar", type=float, default=5.0)
    parser.add_argument("--siglnR", type=float, default=0.693)
    parser.add_argument("--dpeak-fid", type=float, default=2.9e-5)
    parser.add_argument("--siglnl-scale", type=float, default=1.0)
    parser.add_argument("--delta-tau", type=float, default=0.084)
    parser.add_argument("--delta-tau-fid", type=float, default=0.084)
    parser.add_argument("--z-max", type=float, default=40.0)
    parser.add_argument("--n-z", type=int, default=5000)
    parser.add_argument("--z-min-a", type=float, default=0.1)
    parser.add_argument("--g-values", type=float, nargs="+", default=[1.4e-12, 4.0e-12])
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/20-atau-normalization-check"),
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

    matched = nearest_row(read_csv(args.matched_csv), args.mass)
    mass = float(matched["mass_eV"])
    A_unit = float(matched["A_unit"])
    phi_amp_max = float(matched["phi_amp_max"])
    A_tau_raw = A_unit * phi_amp_max

    L = np.arange(args.L_min, args.L_max + 1)
    z_grid = np.linspace(0.0, args.z_max, args.n_z)
    _, z_re = MOD15.a_raw_of_tau_delta(args.tau, args.delta_y, z_grid, args.z_min_a)
    chi_re = float(np.interp(z_re, z_grid, MOD15.cumulative_trapezoid_zero(MOD15.C_KM_S / MOD15.H_of_z(z_grid), z_grid)))
    reff = args.Rbar * np.exp(4.0 * args.siglnR**2)
    l_peak = chi_re / reff
    sigma_lnL = args.siglnl_scale * args.siglnR
    amp_scale = (args.delta_tau / args.delta_tau_fid) ** 2

    d_tau = amp_scale * args.dpeak_fid * MOD15.d_shape_lognormal(L, l_peak, sigma_lnL)
    c_tau = MOD15.c_from_d(L, d_tau)
    c_lim = claa_limit(L, args.A_cb_limit)

    rows: list[dict[str, float | str]] = []
    peak_idx = int(np.argmax(d_tau))
    L_peak_grid = int(L[peak_idx])
    c_tau_peak = float(c_tau[peak_idx])
    c_lim_peak = float(c_lim[peak_idx])

    for g in args.g_values:
        # physical if we naively assume the current field unit is the canonical one that g should multiply
        a_tau_naive = 0.5 * g * A_tau_raw
        c_alpha_peak_naive = (a_tau_naive**2) * c_tau_peak

        # If the ODE field is actually phi_code and phi_phys = U_phi * phi_code,
        # then A_tau_phys = (g/2) * U_phi * A_tau_raw.
        # The U_phi needed to saturate the limit at the peak is:
        uphi_limit_peak = 2.0 * np.sqrt(c_lim_peak / (g**2 * A_tau_raw**2 * c_tau_peak))

        rows.append(
            {
                "mass_eV": mass,
                "g_GeV_inv": g,
                "A_unit": A_unit,
                "phi_amp_max_code_units": phi_amp_max,
                "A_tau_raw_code_units": A_tau_raw,
                "tau": args.tau,
                "delta_y": args.delta_y,
                "delta_tau": args.delta_tau,
                "Reff_Mpc": reff,
                "z_re": z_re,
                "L_peak_surrogate": l_peak,
                "L_peak_grid": L_peak_grid,
                "C_tau_tau_peak": c_tau_peak,
                "C_alphaalpha_limit_peak": c_lim_peak,
                "A_tau_naive_if_uphi_eq_1": a_tau_naive,
                "C_alphaalpha_peak_naive_if_uphi_eq_1": c_alpha_peak_naive,
                "Uphi_limit_peak": uphi_limit_peak,
            }
        )

    write_csv(args.output_dir / "atau_normalization_summary.csv", list(rows[0].keys()), rows)

    # show how the required field-unit factor varies with L
    plt.figure(figsize=(8, 5.5))
    for g in args.g_values:
        uphi_vs_L = 2.0 * np.sqrt(c_lim / (g**2 * A_tau_raw**2 * c_tau))
        plt.semilogy(L, uphi_vs_L, lw=2.0, label=rf"$g={g:.1e}\,{{\rm GeV}}^{{-1}}$")
    plt.xlabel(r"$L$")
    plt.ylabel(r"required field-unit factor $U_\phi$")
    plt.title(r"Field-unit factor needed so $(g/2)^2 U_\phi^2 A_{\tau,\rm raw}^2 C_L^{\tau\tau}$ matches the limit")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir / "uphi_required_vs_L.png", dpi=160, bbox_inches="tight")
    plt.close()

    summary_lines = [
        "# 20-check_atau_normalization summary",
        "",
        "This run asks whether the current ODE field normalization can be directly multiplied by `g_{a\\gamma}`.",
        "",
        "Definitions:",
        "- code/raw quantity: `A_tau_raw = phi_amp_max * A_unit`",
        "- naive physical quantity if the code field were already canonical: `A_tau = (g/2) * A_tau_raw`",
        "- more generally, if `phi_phys = U_phi * phi_code`, then `A_tau = (g/2) * U_phi * A_tau_raw`",
        "",
        f"- matched mass: `{mass:.6e} eV`",
        f"- `A_unit = {A_unit:.6e}`",
        f"- `phi_amp_max` (code units) = `{phi_amp_max:.6e}`",
        f"- `A_tau_raw` (code units) = `{A_tau_raw:.6e}`",
        f"- template anchor: `tau={args.tau:.3f}`, `Delta y={args.delta_y:.1f}`, `delta_tau={args.delta_tau:.3f}`, `R_eff={reff:.3f} Mpc`",
        f"- peak at `L={L_peak_grid}` with `C_tau_tau={c_tau_peak:.6e}` and `C_alphaalpha_limit={c_lim_peak:.6e}`",
        "",
        "Interpretation:",
        "- If `U_phi_limit_peak` is extremely small or extremely large, then the present code-field normalization is unlikely to be the canonical field normalization that should be used with `g` directly.",
        "- This does not by itself prove the physics is wrong; it isolates a missing unit-conversion / canonical-normalization step.",
        "",
    ]
    for row in rows:
        summary_lines.append(
            f"- `g={row['g_GeV_inv']:.2e} GeV^-1`: "
            f"`A_tau_naive={row['A_tau_naive_if_uphi_eq_1']:.6e}`, "
            f"`C_peak_naive={row['C_alphaalpha_peak_naive_if_uphi_eq_1']:.6e}`, "
            f"`U_phi_limit_peak={row['Uphi_limit_peak']:.6e}`"
        )
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps({k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}, indent=2)
    )
    print(f"Saved outputs to: {args.output_dir}")


if __name__ == "__main__":
    main()
