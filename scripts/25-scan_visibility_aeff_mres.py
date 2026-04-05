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


def build_visibility_grid(z_rei: float, delta_y: float, z_max: float, n_z: int) -> dict[str, np.ndarray | float]:
    z_arr = np.linspace(0.0, z_max, n_z)
    a_arr = 1.0 / (1.0 + z_arr)
    x_arr = np.log(a_arr)
    H_arr = H_of_a(a_arr)

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
        "xe_arr": xe_arr,
        "g_z_arr": g_z_arr,
        "norm_g": norm_g,
    }


def lambda_osc_Mpc(m_eV: float, z_rei: float) -> float:
    a_rei = 1.0 / (1.0 + z_rei)
    m_SI = m_eV * EV_TO_J / HBAR
    return 2.0 * np.pi * C_SI / (m_SI * a_rei) / MPC_TO_M


def m_res_fiducial(z_rei: float, r_eff: float) -> float:
    a_rei = 1.0 / (1.0 + z_rei)
    omega_factor = 2.0 * np.pi * HBAR * C_SI / (EV_TO_J * a_rei)
    return omega_factor / (r_eff * MPC_TO_M)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--z-rei", type=float, default=Z_REI)
    p.add_argument("--delta-y", type=float, default=19.0)
    p.add_argument("--z-max", type=float, default=25.0)
    p.add_argument("--n-z", type=int, default=4000)
    p.add_argument("--z-ini", type=float, default=1.0e7)
    p.add_argument("--r-eff", type=float, default=34.139)
    p.add_argument("--factor-min", type=float, default=0.5)
    p.add_argument("--factor-max", type=float, default=2.0)
    p.add_argument("--num-mass", type=int, default=31)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--output-dir", type=Path, default=Path("results/25-visibility-mres-fine"))
    return p


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "visibility_aeff_mres_scan.csv"
    if csv_path.exists() and not args.resume:
        csv_path.unlink()

    m_res = m_res_fiducial(args.z_rei, args.r_eff)
    masses = np.logspace(
        np.log10(args.factor_min * m_res),
        np.log10(args.factor_max * m_res),
        args.num_mass,
    )

    grid = build_visibility_grid(args.z_rei, args.delta_y, args.z_max, args.n_z)
    config = SolveConfig(
        z_ini=args.z_ini,
        method="DOP853",
        rtol=1.0e-9,
        atol_rel=1.0e-12,
        dense_output=True,
    )
    x_eval = float(np.log(1.0 / (1.0 + args.z_rei)))

    fieldnames: list[str] | None = None
    completed_masses: set[float] = set()

    if args.resume and csv_path.exists():
        with csv_path.open(newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            for row in reader:
                completed_masses.add(float(row["mass_eV"]))

    if completed_masses:
        masses = np.array([m for m in masses if not any(np.isclose(m, done, rtol=0.0, atol=1.0e-40) for done in completed_masses)])
        print(f"[25] resume: skipping {len(completed_masses)} completed masses")

    if masses.size == 0:
        print(f"No remaining masses. Existing CSV is already complete: {csv_path}")
        print(f"m_res fiducial = {m_res:.6e} eV")
        return

    start_idx = len(completed_masses)
    total_count = start_idx + len(masses)
    delta_chi_patchy = (C_SI / (H_of_a(np.exp(x_eval)) * MPC_TO_M)) * 3.0

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

        a_eval = np.exp(x_eval)
        A_unit = float(a_eval * H_of_a(a_eval) * dphi_dx_eval * d_eta_d_tau(a_eval))
        A_eff = float(np.trapezoid(grid["g_z_arr"] * response_arr, grid["z_arr"]) / grid["norm_g"])
        ratio = A_eff / A_unit if A_unit != 0.0 else np.nan
        lam = lambda_osc_Mpc(mass, args.z_rei)

        row = {
            "mass_eV": float(mass),
            "m_res_fid_eV": float(m_res),
            "mass_over_mres": float(mass / m_res),
            "z_rei": args.z_rei,
            "delta_y": args.delta_y,
            "R_eff_Mpc": args.r_eff,
            "phi_unit_at_zrei": phi_eval,
            "dphi_dx_unit_at_zrei": dphi_dx_eval,
            "A_unit_thinshell": A_unit,
            "A_eff_visibility": A_eff,
            "Aeff_over_Aunit": float(ratio),
            "lambda_osc_Mpc": float(lam),
            "N_osc_patchy": float(delta_chi_patchy / lam),
            "g_norm_z": float(grid["norm_g"]),
        }

        if fieldnames is None:
            fieldnames = list(row.keys())
            with csv_path.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
        with csv_path.open("a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writerow(row)
        print(f"[25] {start_idx + offset}/{total_count}: m/mres={mass/m_res:.6f}, Aeff/Aunit={ratio:.6e}")

    print(f"Saved: {csv_path}")
    print(f"m_res fiducial = {m_res:.6e} eV")


if __name__ == "__main__":
    main()
