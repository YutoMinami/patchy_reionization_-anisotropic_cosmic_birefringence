from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import cumulative_trapezoid

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patchy_reionization import (
    EV_TO_J,
    F_H,
    HBAR,
    H_of_a,
    M_P,
    MPC_TO_M,
    OMEGA_B,
    RHO_CRIT0,
    SIGMA_T,
    SolveConfig,
    Z_REI,
    d_eta_d_tau,
    solve_phi_background,
)

C_SI = 299792458.0
N_H0 = F_H * OMEGA_B * RHO_CRIT0 / M_P


def x_e_tanh(z_arr: np.ndarray, z_rei: float, delta_y: float) -> np.ndarray:
    y = (1.0 + z_arr) ** 1.5
    y_rei = (1.0 + z_rei) ** 1.5
    return 0.5 * (1.0 + np.tanh((y_rei - y) / delta_y))


def build_visibility_grid(z_rei: float, delta_y: float, z_max: float, n_z: int):
    z_arr = np.linspace(0.0, z_max, n_z)
    a_arr = 1.0 / (1.0 + z_arr)
    x_arr = np.log(a_arr)
    H_arr = H_of_a(a_arr)
    dchi_dz = C_SI / (H_arr * MPC_TO_M)
    chi_arr = np.zeros_like(z_arr)
    chi_arr[1:] = cumulative_trapezoid(dchi_dz, z_arr)

    xe_arr = x_e_tanh(z_arr, z_rei=z_rei, delta_y=delta_y)
    n_e_arr = xe_arr * N_H0 * (1.0 + z_arr) ** 3
    dtau_dz_arr = n_e_arr * SIGMA_T * C_SI / H_arr / (1.0 + z_arr)
    tau_arr = np.zeros_like(z_arr)
    tau_arr[1:] = cumulative_trapezoid(dtau_dz_arr, z_arr)
    g_z_arr = dtau_dz_arr * np.exp(-tau_arr)
    norm_g = np.trapezoid(g_z_arr, z_arr)

    return {
        "z_arr": z_arr,
        "a_arr": a_arr,
        "x_arr": x_arr,
        "H_arr": H_arr,
        "chi_arr": chi_arr,
        "xe_arr": xe_arr,
        "g_z_arr": g_z_arr,
        "norm_g": norm_g,
    }


def lambda_osc_Mpc(m_eV: float, z_rei: float) -> float:
    a_rei = 1.0 / (1.0 + z_rei)
    m_SI = m_eV * EV_TO_J / HBAR
    return 2.0 * np.pi * C_SI / (m_SI * a_rei) / MPC_TO_M


def build_mass_list(args: argparse.Namespace) -> tuple[np.ndarray, float, float, float]:
    a_rei = 1.0 / (1.0 + args.z_rei)
    omega_factor = 2.0 * np.pi * HBAR * C_SI / (EV_TO_J * a_rei)
    m_res_fid = omega_factor / (args.r_eff * MPC_TO_M)
    m_res_min = omega_factor / (args.r_eff_max * MPC_TO_M)
    m_res_max = omega_factor / (args.r_eff_min * MPC_TO_M)

    if args.masses:
        masses = np.array(sorted(set(float(m) for m in args.masses)))
    else:
        m_scan = np.logspace(np.log10(m_res_min / args.m_res_pad), np.log10(m_res_max * args.m_res_pad), args.num_m_res)
        masses = np.array(sorted(set(np.concatenate(([args.m_best, m_res_fid, m_res_min, m_res_max], m_scan)).tolist())))
    if args.mass_min is not None:
        masses = masses[masses >= args.mass_min]
    if args.mass_max is not None:
        masses = masses[masses <= args.mass_max]
    return masses, m_res_fid, m_res_min, m_res_max


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--z-rei", type=float, default=Z_REI)
    p.add_argument("--delta-y", type=float, default=19.0)
    p.add_argument("--z-max", type=float, default=25.0)
    p.add_argument("--n-z", type=int, default=4000)
    p.add_argument("--z-ini", type=float, default=1.0e7)
    p.add_argument("--m-best", type=float, default=5.878016072274924e-27)
    p.add_argument("--r-eff", type=float, default=34.139)
    p.add_argument("--r-eff-min", type=float, default=5.0)
    p.add_argument("--r-eff-max", type=float, default=50.0)
    p.add_argument("--m-res-pad", type=float, default=3.0)
    p.add_argument("--num-m-res", type=int, default=25)
    p.add_argument("--masses", type=float, nargs="+", default=None)
    p.add_argument("--mass-min", type=float, default=None)
    p.add_argument("--mass-max", type=float, default=None)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--output-dir", type=Path, default=Path("results/24-visibility-averaged-aeff"))
    return p


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "visibility_averaged_aeff_scan.csv"
    if csv_path.exists() and not args.resume:
        csv_path.unlink()

    grid = build_visibility_grid(args.z_rei, args.delta_y, args.z_max, args.n_z)
    masses, m_res_fid, m_res_min, m_res_max = build_mass_list(args)

    config = SolveConfig(
        z_ini=args.z_ini,
        method="DOP853",
        rtol=1.0e-9,
        atol_rel=1.0e-12,
        dense_output=True,
    )

    x_eval = float(np.log(1.0 / (1.0 + args.z_rei)))
    rows: list[dict[str, float]] = []
    fieldnames: list[str] | None = None
    completed_masses: set[float] = set()

    if args.resume and csv_path.exists():
        with csv_path.open(newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            for row in reader:
                rows.append(row)
                completed_masses.add(float(row["mass_eV"]))

    if completed_masses:
        masses = np.array([m for m in masses if not any(np.isclose(m, done, rtol=0.0, atol=1.0e-40) for done in completed_masses)])
        print(f"[24] resume: skipping {len(completed_masses)} completed masses")

    if masses.size == 0:
        print(f"No remaining masses. Existing CSV is already complete: {csv_path}")
        print(f"m_res fiducial = {m_res_fid:.6e} eV")
        print(f"m_res range    = {m_res_min:.6e} -- {m_res_max:.6e} eV")
        return

    start_idx = len(completed_masses)
    total_count = start_idx + len(masses)
    for offset, mass in enumerate(masses, start=1):
        sol = solve_phi_background(m_eV=mass, phi_ini=1.0, config=config)
        state_grid = np.asarray(sol.sol(grid["x_arr"]))
        dphi_dx_arr = state_grid[1]
        state_eval = np.asarray(sol.sol(x_eval))
        phi_eval = float(state_eval[0])
        dphi_dx_eval = float(state_eval[1])

        dotphi_conf_arr = grid["a_arr"] * grid["H_arr"] * dphi_dx_arr
        deta_dtau_arr = d_eta_d_tau(grid["a_arr"], x_e=grid["xe_arr"])
        response_arr = dotphi_conf_arr * deta_dtau_arr

        A_unit = float(grid["a_arr"][np.argmin(np.abs(grid["x_arr"] - x_eval))] * 0.0)
        A_unit = float(np.exp(x_eval) * H_of_a(np.exp(x_eval)) * dphi_dx_eval * d_eta_d_tau(np.exp(x_eval)))
        A_eff = float(np.trapezoid(grid["g_z_arr"] * response_arr, grid["z_arr"]) / grid["norm_g"])
        ratio = A_eff / A_unit if A_unit != 0.0 else np.nan

        lam = lambda_osc_Mpc(mass, args.z_rei)
        delta_chi_patchy = (C_SI / (H_of_a(np.exp(x_eval)) * MPC_TO_M)) * 3.0
        N_osc_patchy = delta_chi_patchy / lam

        row = {
            "mass_eV": float(mass),
            "z_rei": args.z_rei,
            "delta_y": args.delta_y,
            "R_eff_Mpc": args.r_eff,
            "m_best_eV": args.m_best,
            "m_res_fid_eV": float(m_res_fid),
            "m_res_min_eV": float(m_res_min),
            "m_res_max_eV": float(m_res_max),
            "phi_unit_at_zrei": phi_eval,
            "dphi_dx_unit_at_zrei": dphi_dx_eval,
            "A_unit_thinshell": A_unit,
            "A_eff_visibility": A_eff,
            "Aeff_over_Aunit": float(ratio),
            "lambda_osc_Mpc": float(lam),
            "N_osc_patchy": float(N_osc_patchy),
            "g_norm_z": float(grid["norm_g"]),
        }
        rows.append(row)
        if fieldnames is None:
            fieldnames = list(row.keys())
            with csv_path.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
        with csv_path.open("a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writerow(row)
        print(f"[24] {start_idx + offset}/{total_count}: m={mass:.6e} eV, Aeff/Aunit={ratio:.6e}")

    print(f"Saved: {csv_path}")
    print(f"m_res fiducial = {m_res_fid:.6e} eV")
    print(f"m_res range    = {m_res_min:.6e} -- {m_res_max:.6e} eV")


if __name__ == "__main__":
    main()
