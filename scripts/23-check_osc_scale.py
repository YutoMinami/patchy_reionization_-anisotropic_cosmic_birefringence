"""23-check_osc_scale.py

Step 1 of the Caveat 2 verification (HANDOFF.md):
Compute the ALP oscillation comoving wavelength lambda_osc(m_a) at z_rei
and compare with the reionization bubble scale R_eff and the
visibility-function comoving width Delta_chi_vis.

Key quantities:
  lambda_osc(m_a) = 2pi * hbar * c / (m_a [eV] * eV_to_J / hbar * a_rei)
  N_osc = Delta_chi_vis / lambda_osc   (number of oscillations in visibility window)

If N_osc >> 1, many ALP oscillations fit in the visibility window and the
thin-shell approximation may significantly overestimate A_eff. If N_osc << 1,
the thin-shell is accurate.
"""
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
    Z_REI,
)

C_SI = 299792458.0

# comoving number density of H nucleons [m^{-3}]
N_H0 = F_H * OMEGA_B * RHO_CRIT0 / M_P


# ── Oscillation wavelength ─────────────────────────────────────────────────

def lambda_osc_Mpc(m_eV: float | np.ndarray, z_rei: float = Z_REI) -> float | np.ndarray:
    """Comoving ALP oscillation wavelength at z_rei [Mpc]. Vectorised."""
    a_rei = 1.0 / (1.0 + z_rei)
    m_SI = np.asarray(m_eV) * EV_TO_J / HBAR  # rad / s
    return 2.0 * np.pi * C_SI / (m_SI * a_rei) / MPC_TO_M


# ── Visibility function ────────────────────────────────────────────────────

def x_e_tanh(z_arr: np.ndarray, z_rei: float, delta_y: float) -> np.ndarray:
    y = (1.0 + z_arr) ** 1.5
    y_rei = (1.0 + z_rei) ** 1.5
    return 0.5 * (1.0 + np.tanh((y_rei - y) / delta_y))


def build_chi_grid(
    z_max: float = 25.0, n_z: int = 4000
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (z_arr, chi_arr [Mpc], H_arr [SI]) over a uniform z grid.

    H_arr is returned so callers can reuse it without a second loop.
    """
    z_arr = np.linspace(0.0, z_max, n_z)
    a_arr = 1.0 / (1.0 + z_arr)
    H_arr = H_of_a(a_arr)                          # vectorised — one call
    dchi_dz = C_SI / (H_arr * MPC_TO_M)
    chi_arr = np.zeros(n_z)
    chi_arr[1:] = cumulative_trapezoid(dchi_dz, z_arr)
    return z_arr, chi_arr, H_arr


def compute_visibility_width(
    z_rei: float = Z_REI,
    delta_y: float = 19.0,
    z_max: float = 25.0,
    n_z: int = 4000,
) -> dict[str, float]:
    """Compute comoving width of the visibility function from tanh model.

    Returns dict with chi_rei, sigma_chi, fwhm_chi, tau_total [Mpc / dimensionless].
    """
    z_arr, chi_arr, H_arr = build_chi_grid(z_max=z_max, n_z=n_z)

    xe = x_e_tanh(z_arr, z_rei=z_rei, delta_y=delta_y)
    n_e = xe * N_H0 * (1.0 + z_arr) ** 3          # physical n_e [m^{-3}]

    # dτ/dz = n_e σ_T c / (H(z) (1+z))
    dtau_dz_arr = n_e * SIGMA_T * C_SI / H_arr / (1.0 + z_arr)

    tau_arr = np.zeros(n_z)
    tau_arr[1:] = cumulative_trapezoid(dtau_dz_arr, z_arr)

    g_z = dtau_dz_arr * np.exp(-tau_arr)

    idx_rei = np.argmin(np.abs(z_arr - z_rei))
    chi_rei = chi_arr[idx_rei]

    norm = np.trapezoid(g_z, z_arr)
    sigma_chi = np.sqrt(np.trapezoid((chi_arr - chi_rei) ** 2 * g_z, z_arr) / norm)
    fwhm_chi = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma_chi

    return {
        "chi_rei": chi_rei,
        "sigma_chi": sigma_chi,
        "fwhm_chi": fwhm_chi,
        "tau_total": float(tau_arr[-1]),
    }


# ── Argument parsing ───────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--z-rei", type=float, default=Z_REI)
    p.add_argument("--delta-y", type=float, default=19.0,
                   help="Tanh ionization history width parameter Delta_y")
    p.add_argument("--r-eff", type=float, default=34.139,
                   help="Representative bubble effective radius R_eff [Mpc] (from run 22)")
    p.add_argument("--r-eff-min", type=float, default=5.0,
                   help="Lower end of R_eff range for shading [Mpc]")
    p.add_argument("--r-eff-max", type=float, default=50.0,
                   help="Upper end of R_eff range for shading [Mpc]")
    p.add_argument("--m-min", type=float, default=1e-33,
                   help="Minimum ALP mass [eV]")
    p.add_argument("--m-max", type=float, default=1e-24,
                   help="Maximum ALP mass [eV]")
    p.add_argument("--num-mass", type=int, default=200,
                   help="Number of mass grid points (log-spaced)")
    p.add_argument("--m-best", type=float, default=5.878016e-27,
                   help="Reference mass from run 04a/11 [eV]")
    p.add_argument("--delta-z-patchy", type=float, default=3.0,
                   help="Width of the patchy reionization epoch [Delta z]")
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/23-osc-scale"))
    p.add_argument("--no-show", action="store_true")
    return p


# ── Main ───────────────────────────────────────────────────────────────────

def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Computing visibility function width...")
    vis = compute_visibility_width(z_rei=args.z_rei, delta_y=args.delta_y)
    chi_rei = vis["chi_rei"]
    sigma_chi = vis["sigma_chi"]
    fwhm_chi = vis["fwhm_chi"]
    tau_total = vis["tau_total"]
    print(f"  chi_rei   = {chi_rei:.1f} Mpc")
    print(f"  sigma_chi = {sigma_chi:.1f} Mpc  (1-sigma width of g(chi))")
    print(f"  fwhm_chi  = {fwhm_chi:.1f} Mpc")
    print(f"  tau_total = {tau_total:.4f}")

    delta_chi_vis = 2.0 * sigma_chi

    # Patchy epoch comoving width
    a_rei = 1.0 / (1.0 + args.z_rei)
    dchi_dz_rei = C_SI / (H_of_a(a_rei) * MPC_TO_M)
    delta_chi_patchy = dchi_dz_rei * args.delta_z_patchy
    delta_chi_patchy_m = delta_chi_patchy * MPC_TO_M

    # Temporal crossover: N_osc = 1 within patchy epoch
    m_crossover = 2.0 * np.pi * HBAR * C_SI / (EV_TO_J * a_rei * delta_chi_patchy_m)

    # Spatial resonance masses: lambda_osc = R_eff
    omega_factor = 2.0 * np.pi * HBAR * C_SI / (EV_TO_J * a_rei)
    m_res_fiducial = omega_factor / (args.r_eff * MPC_TO_M)
    m_res_min = omega_factor / (args.r_eff_max * MPC_TO_M)
    m_res_max = omega_factor / (args.r_eff_min * MPC_TO_M)

    # Mass grid — vectorised lambda_osc
    m_arr = np.logspace(np.log10(args.m_min), np.log10(args.m_max), args.num_mass)
    lam_arr = lambda_osc_Mpc(m_arr, z_rei=args.z_rei)
    N_osc_arr = delta_chi_vis / lam_arr
    N_osc_patchy_arr = delta_chi_patchy / lam_arr

    # Values at reference mass
    lam_best = lambda_osc_Mpc(args.m_best, z_rei=args.z_rei)
    N_osc_best = delta_chi_vis / lam_best
    N_osc_patchy_best = delta_chi_patchy / lam_best

    print(f"  delta_chi_vis (2sigma)        = {delta_chi_vis:.1f} Mpc")
    print(f"  dchi/dz at z_rei              = {dchi_dz_rei:.1f} Mpc/dz")
    print(f"  delta_chi_patchy (Dz={args.delta_z_patchy})    = {delta_chi_patchy:.1f} Mpc")
    print(f"  m_crossover (N_osc=1)         = {m_crossover:.3e} eV")
    print(f"  m_res (lambda=R_eff fiducial) = {m_res_fiducial:.3e} eV")
    print(f"  m_res range                   = {m_res_min:.2e}–{m_res_max:.2e} eV")
    print(f"  lambda_osc at m_best          = {lam_best:.4f} Mpc = {lam_best*1e3:.1f} kpc")
    print(f"  N_osc at m_best (full vis)    = {N_osc_best:.1f}")
    print(f"  N_osc at m_best (patchy)      = {N_osc_patchy_best:.1f}")

    # ── Save CSV ───────────────────────────────────────────────────────────
    csv_path = args.output_dir / "osc_scale_scan.csv"
    fieldnames = ["mass_eV", "lambda_osc_Mpc", "N_osc_full", "N_osc_patchy"]
    rows = [
        {"mass_eV": m, "lambda_osc_Mpc": lam, "N_osc_full": Nf, "N_osc_patchy": Np}
        for m, lam, Nf, Np in zip(m_arr, lam_arr, N_osc_arr, N_osc_patchy_arr)
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {csv_path}")

    # ── Plot ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(8, 9), sharex=True)
    dz_p = args.delta_z_patchy

    ax = axes[0]
    ax.loglog(m_arr, lam_arr, color="C0", lw=2, label=r"$\lambda_{\rm osc}(m_a)$")
    ax.axhspan(args.r_eff_min, args.r_eff_max, alpha=0.15, color="C1",
               label=rf"$R_{{\rm eff}}$ range ({args.r_eff_min}–{args.r_eff_max} Mpc)")
    ax.axhline(args.r_eff, color="C1", ls="--", lw=1.5,
               label=rf"$R_{{\rm eff}} = {args.r_eff:.1f}$ Mpc (fiducial)")
    ax.axhline(delta_chi_patchy, color="C3", ls="-", lw=1.5,
               label=rf"$\Delta\chi_{{\rm patchy}} = {delta_chi_patchy:.0f}$ Mpc ($\Delta z={dz_p}$)")
    ax.axhline(delta_chi_vis, color="C2", ls=":", lw=1.5,
               label=rf"$\Delta\chi_{{\rm vis}} = {delta_chi_vis:.0f}$ Mpc (full tanh $2\sigma$)")
    ax.axvline(args.m_best, color="gray", ls="--", lw=1,
               label=rf"$m_{{\rm best}} = {args.m_best:.2e}$ eV")
    ax.axvline(m_crossover, color="red", ls=":", lw=1.5,
               label=rf"$m_* = {m_crossover:.1e}$ eV ($N_{{\rm osc}}=1$, temporal)")
    ax.axvspan(m_res_min, m_res_max, alpha=0.12, color="C4",
               label=rf"$m_{{\rm res}}$ band ({m_res_min:.0e}–{m_res_max:.0e} eV)")
    ax.axvline(m_res_fiducial, color="C4", ls="-.", lw=1.5,
               label=rf"$m_{{\rm res}} = {m_res_fiducial:.1e}$ eV (spatial resonance)")
    ax.set_ylabel(r"$\lambda_{\rm osc}$ [Mpc]")
    ax.set_title(r"ALP oscillation scale vs bubble / visibility scales")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[1]
    ax.loglog(m_arr, N_osc_arr, color="C2", lw=1.5, ls=":",
              label=r"$N_{\rm osc}$ (full tanh $2\sigma$)")
    ax.loglog(m_arr, N_osc_patchy_arr, color="C3", lw=2,
              label=rf"$N_{{\rm osc}}$ (patchy $\Delta z={dz_p}$)")
    ax.axhline(1.0, color="red", ls=":", lw=1.5, label=r"$N_{\rm osc} = 1$")
    ax.axvline(args.m_best, color="gray", ls="--", lw=1,
               label=rf"$m_{{\rm best}} = {args.m_best:.2e}$ eV")
    ax.axvspan(m_res_min, m_res_max, alpha=0.12, color="C4",
               label=r"$m_{\rm res}$ band")
    ax.scatter([args.m_best], [N_osc_patchy_best], color="C3", zorder=5, s=60)
    ax.annotate(
        rf"$N_{{\rm osc}} = {N_osc_patchy_best:.0f}$",
        xy=(args.m_best, N_osc_patchy_best),
        xytext=(args.m_best * 8, N_osc_patchy_best * 4),
        fontsize=9,
        arrowprops=dict(arrowstyle="->", color="C3"),
    )
    ax.set_xlabel(r"$m_a$ [eV]")
    ax.set_ylabel(r"$N_{\rm osc}$")
    ax.set_title(r"Number of ALP oscillations within reionization window")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    plot_path = args.output_dir / "osc_scale_vs_m.png"
    plt.savefig(plot_path, dpi=160, bbox_inches="tight")
    if not args.no_show:
        plt.show()
    plt.close()
    print(f"Saved: {plot_path}")

    # ── Summary ────────────────────────────────────────────────────────────
    lines = [
        "# 23-check_osc_scale summary",
        "",
        "Step 1 of Caveat 2 verification: compare lambda_osc(m_a) to R_eff and Delta_chi_vis.",
        "",
        "## Visibility function (tanh model)",
        f"- `z_rei = {args.z_rei}`, `delta_y = {args.delta_y}`",
        f"- `chi_rei = {chi_rei:.1f} Mpc`",
        f"- `sigma_chi = {sigma_chi:.1f} Mpc`  (1-sigma comoving width of g(chi))",
        f"- `fwhm_chi = {fwhm_chi:.1f} Mpc`",
        f"- `tau_total = {tau_total:.4f}`",
        f"- `delta_chi_vis (2 sigma, full tanh) = {delta_chi_vis:.1f} Mpc`",
        "",
        "## Patchy reionization epoch width",
        f"- `delta_z_patchy = {args.delta_z_patchy}`",
        f"- `dchi/dz at z_rei = {dchi_dz_rei:.1f} Mpc/dz`",
        f"- `delta_chi_patchy = {delta_chi_patchy:.1f} Mpc`",
        "",
        "## Bubble scale reference",
        f"- `R_eff (fiducial, from run 22) = {args.r_eff:.1f} Mpc`",
        f"- `R_eff range = {args.r_eff_min}--{args.r_eff_max} Mpc`",
        "",
        "## At the reference mass",
        f"- `m_best = {args.m_best:.6e} eV`",
        f"- `lambda_osc = {lam_best:.4f} Mpc = {lam_best * 1e3:.1f} kpc`",
        f"- `N_osc (full vis)   = {N_osc_best:.1f}`",
        f"- `N_osc (patchy)     = {N_osc_patchy_best:.1f}`",
        "",
        "## Crossover mass (N_osc = 1, patchy epoch)",
        f"- `m_crossover = {m_crossover:.3e} eV`",
        "",
        "## Spatial resonance mass (lambda_osc = R_eff)",
        f"- `m_res (fiducial, R_eff={args.r_eff} Mpc) = {m_res_fiducial:.3e} eV`",
        f"- `m_res range ({args.r_eff_min}--{args.r_eff_max} Mpc) = {m_res_min:.2e}--{m_res_max:.2e} eV`",
        "  m_res is ~3 orders of magnitude below m_best.",
        "",
        "## Interpretation",
        "- At `m_best`, `lambda_osc << R_eff`: no spatial resonance at the reference mass.",
        "- Spatial resonance occurs at `m_res ~ 1e-29 eV`; A_unit there is comparable to m_best.",
        "- `N_osc >> 1` at `m_best`: dominant concern is temporal phase averaging.",
        "  The full visibility-weighted A_eff is quantified in script 24.",
    ]
    (args.output_dir / "summary.md").write_text("\n".join(lines) + "\n")
    print(f"Saved: {args.output_dir / 'summary.md'}")

    (args.output_dir / "run_config.json").write_text(
        json.dumps(
            {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
            indent=2,
        )
    )
    print(f"Saved: {args.output_dir / 'run_config.json'}")


if __name__ == "__main__":
    main()
