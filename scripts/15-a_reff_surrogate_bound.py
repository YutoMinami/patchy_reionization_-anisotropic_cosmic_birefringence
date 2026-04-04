from __future__ import annotations

import argparse
import csv
import json
from itertools import product
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

from patchy_reionization import F_H, H0_KM_S_MPC, H0_SI, M_P, OMEGA_L, OMEGA_M, OMEGA_R, RHO_B0, SIGMA_T


C_KM_S = 299792.458
C_SI = 299792458.0
N_H0_SI = F_H * RHO_B0 / M_P


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
    parser.add_argument("--tau-values", type=float, nargs="+", default=[0.068, 0.11])
    parser.add_argument("--delta-y-values", type=float, nargs="+", default=[14.7, 19.0])
    parser.add_argument("--b-values", type=float, nargs="+", default=[1.0, 2.0])
    parser.add_argument("--rbar-values", type=float, nargs="+", default=[5.0, 6.0, 8.0])
    parser.add_argument("--siglnr-values", type=float, nargs="+", default=[0.409, 0.541, 0.693])
    parser.add_argument("--z-max", type=float, default=40.0)
    parser.add_argument("--n-z", type=int, default=5000)
    parser.add_argument("--z-min-a", type=float, default=0.1)
    parser.add_argument("--a-fid-target", type=float, default=0.07)
    parser.add_argument("--tau-fid", type=float, default=0.068)
    parser.add_argument("--delta-y-fid", type=float, default=19.0)
    parser.add_argument("--b-fid", type=float, default=1.0)
    parser.add_argument("--dpeak-fid", type=float, default=2.5e-5)
    parser.add_argument("--siglnl-scale", type=float, default=1.0)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/15-a-reff-surrogate"),
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


def E_of_z(z: np.ndarray) -> np.ndarray:
    return np.sqrt(OMEGA_R * (1.0 + z) ** 4 + OMEGA_M * (1.0 + z) ** 3 + OMEGA_L)


def H_of_z(z: np.ndarray) -> np.ndarray:
    return H0_KM_S_MPC * E_of_z(z)


def y_of_z(z: np.ndarray | float) -> np.ndarray | float:
    return (1.0 + np.asarray(z)) ** 1.5


def xe_tanh(z: np.ndarray, z_re: float, delta_y: float) -> np.ndarray:
    y = y_of_z(z)
    y_re = y_of_z(z_re)
    return 0.5 * (1.0 + np.tanh((y_re - y) / delta_y))


def cumulative_trapezoid_zero(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    dx = np.diff(x)
    area = 0.5 * (y[1:] + y[:-1]) * dx
    return np.concatenate([[0.0], np.cumsum(area)])


def tau_of_zre(z_re: float, delta_y: float, z_grid: np.ndarray) -> float:
    xe = xe_tanh(z_grid, z_re=z_re, delta_y=delta_y)
    integrand = C_SI * SIGMA_T * N_H0_SI * (1.0 + z_grid) ** 2 * xe / (H0_SI * E_of_z(z_grid))
    return float(np.trapezoid(integrand, z_grid))


def solve_zre_for_tau(tau_target: float, delta_y: float, z_grid: np.ndarray) -> float:
    low, high = 3.0, 30.0
    tau_low = tau_of_zre(low, delta_y, z_grid)
    tau_high = tau_of_zre(high, delta_y, z_grid)
    if not (tau_low <= tau_target <= tau_high):
        # Fall back to clipping if the proxy range misses the requested value.
        return low if abs(tau_target - tau_low) < abs(tau_target - tau_high) else high
    for _ in range(80):
        mid = 0.5 * (low + high)
        tau_mid = tau_of_zre(mid, delta_y, z_grid)
        if tau_mid < tau_target:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def a_raw_of_tau_delta(
    tau_target: float,
    delta_y: float,
    z_grid: np.ndarray,
    z_min_a: float,
) -> tuple[float, float]:
    z_re = solve_zre_for_tau(tau_target=tau_target, delta_y=delta_y, z_grid=z_grid)
    xe = xe_tanh(z_grid, z_re=z_re, delta_y=delta_y)
    chi = cumulative_trapezoid_zero(C_KM_S / H_of_z(z_grid), z_grid)
    mask = z_grid >= z_min_a
    z = z_grid[mask]
    chi_m = chi[mask]
    xe_m = xe[mask]
    integrand = ((1.0 + z) ** 4 / H_of_z(z)) * xe_m * (1.0 - xe_m) / (chi_m**2)
    a_raw = float(np.trapezoid(integrand, z))
    return a_raw, z_re


def d_shape_lognormal(L: np.ndarray, L_peak: float, sigma_lnL: float) -> np.ndarray:
    sigma_lnL = max(sigma_lnL, 1.0e-3)
    return np.exp(-0.5 * (np.log(L / L_peak) / sigma_lnL) ** 2)


def c_from_d(L: np.ndarray, D: np.ndarray) -> np.ndarray:
    return 2.0 * np.pi * D / (L * (L + 1.0))


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matched_rows = read_csv(args.matched_csv)
    selected_rows = nearest_rows(matched_rows, args.selected_masses)
    L = np.arange(args.L_min, args.L_max + 1)
    claa_lim = 2.0 * np.pi * args.A_cb_limit / (L * (L + 1.0))

    z_grid = np.linspace(0.0, args.z_max, args.n_z)
    a_cache: dict[tuple[float, float], tuple[float, float]] = {}
    for tau, delta_y in product(args.tau_values, args.delta_y_values):
        a_cache[(tau, delta_y)] = a_raw_of_tau_delta(
            tau_target=tau,
            delta_y=delta_y,
            z_grid=z_grid,
            z_min_a=args.z_min_a,
        )

    a_raw_fid, z_re_fid = a_cache[(args.tau_fid, args.delta_y_fid)]

    shape_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | str]] = []

    for tau, delta_y, b, rbar, siglnr in product(
        args.tau_values,
        args.delta_y_values,
        args.b_values,
        args.rbar_values,
        args.siglnr_values,
    ):
        a_raw, z_re = a_cache[(tau, delta_y)]
        a_model = args.a_fid_target * (a_raw / a_raw_fid)
        reff_mpc = rbar * np.exp(4.0 * siglnr**2)
        chi_re_mpc = float(np.interp(z_re, z_grid, cumulative_trapezoid_zero(C_KM_S / H_of_z(z_grid), z_grid)))
        l_peak = chi_re_mpc / reff_mpc
        sigma_lnL = args.siglnl_scale * siglnr
        d_peak = args.dpeak_fid * (b / args.b_fid) * (a_model / args.a_fid_target)
        d_shape = d_peak * d_shape_lognormal(L, L_peak=l_peak, sigma_lnL=sigma_lnL)
        c_shape = c_from_d(L, d_shape)

        for ell, d_val, c_val in zip(L, d_shape, c_shape, strict=True):
            shape_rows.append(
                {
                    "tau": tau,
                    "delta_y": delta_y,
                    "b": b,
                    "Rbar_Mpc": rbar,
                    "siglnR": siglnr,
                    "A_model": a_model,
                    "z_re": z_re,
                    "Reff_Mpc": reff_mpc,
                    "L_peak_surrogate": l_peak,
                    "siglnL_surrogate": sigma_lnL,
                    "D_peak_surrogate": d_peak,
                    "L": int(ell),
                    "D_shape": float(d_val),
                    "C_shape": float(c_val),
                }
            )

        for selected in selected_rows:
            mass = float(selected["mass_eV"])
            A_unit = float(selected["A_unit"])
            phi_amp_max = float(selected["phi_amp_max"])
            A_tau_eff = A_unit * phi_amp_max
            patchy_power = (A_tau_eff**2) * c_shape
            budget_fraction = patchy_power / claa_lim
            summary_rows.append(
                {
                    "mass_eV": mass,
                    "tau": tau,
                    "delta_y": delta_y,
                    "b": b,
                    "Rbar_Mpc": rbar,
                    "siglnR": siglnr,
                    "A_model": a_model,
                    "z_re": z_re,
                    "Reff_Mpc": reff_mpc,
                    "L_peak_surrogate": l_peak,
                    "siglnL_surrogate": sigma_lnL,
                    "D_peak_surrogate": d_peak,
                    "A_tau_eff": A_tau_eff,
                    "max_patchy_budget_fraction": float(np.nanmax(budget_fraction)),
                    "L_at_max_fraction": int(L[int(np.nanargmax(budget_fraction))]),
                    "min_scale_factor_to_fit": float(1.0 / np.nanmax(budget_fraction)),
                }
            )

    write_csv(args.output_dir / "surrogate_shapes.csv", list(shape_rows[0].keys()), shape_rows)
    write_csv(args.output_dir / "surrogate_bound_summary.csv", list(summary_rows[0].keys()), summary_rows)

    plt.figure(figsize=(8, 5.5))
    plotted = 0
    for row in summary_rows:
        if float(row["mass_eV"]) != float(selected_rows[0]["mass_eV"]):
            continue
        if plotted >= 8:
            break
        D = d_shape_lognormal(L, float(row["L_peak_surrogate"]), float(row["siglnL_surrogate"])) * float(row["D_peak_surrogate"])
        label = rf"$\tau={float(row['tau']):.3f},\,\Delta y={float(row['delta_y']):g},\,R_{{\rm eff}}={float(row['Reff_Mpc']):.1f}$"
        plt.loglog(L, D, lw=1.6, label=label)
        plotted += 1
    plt.xlabel(r"$L$")
    plt.ylabel(r"$D_L^{\tau\tau}$ surrogate")
    plt.title("Representative A + R_eff surrogate shapes")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.output_dir / "surrogate_shapes_examples.png", dpi=160, bbox_inches="tight")
    plt.close()

    for mass in sorted({float(row["mass_eV"]) for row in summary_rows}):
        subset = [row for row in summary_rows if float(row["mass_eV"]) == mass]
        plt.figure(figsize=(8, 5.5))
        xs = [float(row["Reff_Mpc"]) for row in subset]
        ys = [float(row["min_scale_factor_to_fit"]) for row in subset]
        cs = [float(row["A_model"]) for row in subset]
        sc = plt.scatter(xs, ys, c=cs, cmap="viridis", s=45)
        plt.yscale("log")
        plt.xlabel(r"$R_{\rm eff}\ [{\rm Mpc}]$")
        plt.ylabel("required down-scaling to fit A_CB")
        plt.title(rf"A + $R_{{\rm eff}}$ surrogate at $m_a={mass:.2e}\,{{\rm eV}}$")
        plt.grid(True, which="both", alpha=0.3)
        cbar = plt.colorbar(sc)
        cbar.set_label(r"$A$ surrogate")
        plt.tight_layout()
        plt.savefig(args.output_dir / f"scale_factor_vs_reff_{mass:.3e}.png", dpi=160, bbox_inches="tight")
        plt.close()

    summary_lines = [
        "# 15-a_reff_surrogate_bound summary",
        "",
        "This run replaces the free unit-peak normalization of `14` with a Dvorkin-Smith-inspired surrogate.",
        "",
        "Approximation choices:",
        "- The observable reionization-amplitude combination is modeled by Eq. (77)-inspired `A(tau, Delta y)` using a tanh-like ionization history.",
        "- The characteristic angular scale is modeled through `R_eff = Rbar * exp(4 sigma_lnR^2)` from Eq. (78).",
        "- The shape remains a lognormal bump in `D_L^{tau tau}`, centered at `L_peak ~= chi(z_re)/R_eff`.",
        f"- The surrogate amplitude is calibrated to `A_fid={args.a_fid_target:.3f}` and `D_peak_fid={args.dpeak_fid:.3e}` at `(tau, Delta y)=({args.tau_fid:.3f}, {args.delta_y_fid:.1f})`.",
        "",
    ]
    for mass in sorted({float(row["mass_eV"]) for row in summary_rows}):
        subset = [row for row in summary_rows if float(row["mass_eV"]) == mass]
        best = min(subset, key=lambda row: float(row["min_scale_factor_to_fit"]))
        loosest = max(subset, key=lambda row: float(row["min_scale_factor_to_fit"]))
        summary_lines.append(
            f"- `m={mass:.6e} eV`: strongest tension `scale<{float(best['min_scale_factor_to_fit']):.6e}` "
            f"for `(tau={float(best['tau']):.3f}, Delta y={float(best['delta_y']):g}, b={float(best['b']):g}, "
            f"Rbar={float(best['Rbar_Mpc']):g}, siglnR={float(best['siglnR']):.3f})`; "
            f"weakest tension `scale<{float(loosest['min_scale_factor_to_fit']):.6e}` "
            f"for `(tau={float(loosest['tau']):.3f}, Delta y={float(loosest['delta_y']):g}, b={float(loosest['b']):g}, "
            f"Rbar={float(loosest['Rbar_Mpc']):g}, siglnR={float(loosest['siglnR']):.3f})`."
        )
    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps(
            {
                "matched_csv": str(args.matched_csv),
                "A_cb_limit": args.A_cb_limit,
                "args": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
                "a_cache": {
                    f"tau_{tau:.6f}_dy_{dy:.6f}": {"A_raw": a_raw, "z_re": z_re}
                    for (tau, dy), (a_raw, z_re) in a_cache.items()
                },
            },
            indent=2,
        )
    )

    print(f"Saved outputs to: {args.output_dir}")
    for mass in sorted({float(row['mass_eV']) for row in summary_rows}):
        subset = [row for row in summary_rows if float(row["mass_eV"]) == mass]
        best = min(subset, key=lambda row: float(row["min_scale_factor_to_fit"]))
        loosest = max(subset, key=lambda row: float(row["min_scale_factor_to_fit"]))
        print(
            f"m={mass:.6e} eV: strongest scale<{float(best['min_scale_factor_to_fit']):.6e}, "
            f"weakest scale<{float(loosest['min_scale_factor_to_fit']):.6e}"
        )

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
