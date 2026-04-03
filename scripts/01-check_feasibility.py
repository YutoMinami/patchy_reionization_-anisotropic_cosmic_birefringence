from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patchy_reionization import Cl_phiphi_toy, Cl_tautau_toy, R_tau_toy
from patchy_reionization import SolveConfig, compute_A_unit, rescale_A_from_unit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m-min", type=float, default=1.0e-35)
    parser.add_argument("--m-max", type=float, default=1.0e-27)
    parser.add_argument("--num-mass", type=int, default=60)
    parser.add_argument("--z-ini", type=float, default=1.0e7)
    parser.add_argument("--phi-amp", type=float, default=1.0)
    parser.add_argument("--L-max", type=int, default=1200)
    parser.add_argument("--nphi", type=float, default=2.0)
    parser.add_argument("--Lcut", type=float, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--no-show", action="store_true")
    return parser


def maybe_savefig(output_dir: Path | None, name: str) -> None:
    if output_dir is None:
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / name, dpi=160, bbox_inches="tight")


def main() -> None:
    args = build_parser().parse_args()
    config = SolveConfig(z_ini=args.z_ini)

    m_grid = np.logspace(np.log10(args.m_min), np.log10(args.m_max), args.num_mass)
    A_unit_vals = np.array([compute_A_unit(m_eV=m, config=config) for m in m_grid])
    A_vals = rescale_A_from_unit(A_unit_vals, args.phi_amp)

    plt.figure(figsize=(7, 5))
    plt.loglog(m_grid, np.abs(A_unit_vals), label=r"$|A_{\rm unit}(m_a)|$")
    if args.phi_amp != 1.0:
        plt.loglog(m_grid, np.abs(A_vals), ls="--", label=rf"$|A|$ for $\phi_{{\rm amp}}={args.phi_amp:g}$")
    plt.xlabel(r"$m_a \ [{\rm eV}]$")
    plt.ylabel(r"$|\mathcal{A}| = |\dot{\phi}_{\rm conf} d\eta/d\tau|$")
    plt.title("Patchy-reionization response scan")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    maybe_savefig(args.output_dir, "A_unit_vs_m.png")

    idx = np.nanargmax(np.abs(A_vals))
    m_pick = m_grid[idx]
    A_pick = A_vals[idx]

    L = np.arange(2, args.L_max + 1)
    Cl_tau = Cl_tautau_toy(L)
    Cl_phi = Cl_phiphi_toy(L, nphi=args.nphi, Lcut=args.Lcut)
    R = R_tau_toy(L, A_pick, nphi=args.nphi, Lcut=args.Lcut)

    plt.figure(figsize=(7, 5))
    plt.loglog(L, Cl_tau, label=r"$C_L^{\tau\tau}$")
    plt.loglog(L, Cl_phi, label=r"$C_L^{\phi\phi}$")
    plt.xlabel(r"$L$")
    plt.ylabel(r"$C_L$")
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    maybe_savefig(args.output_dir, "toy_cls.png")

    plt.figure(figsize=(7, 5))
    plt.semilogy(L, np.abs(R))
    plt.axhline(1.0, ls="--", color="black", alpha=0.6)
    plt.xlabel(r"$L$")
    plt.ylabel(r"$R_\tau(L)$")
    plt.title(r"$R_\tau(L)=\mathcal{A}^2 C_L^{\tau\tau}/C_L^{\phi\phi}$")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    maybe_savefig(args.output_dir, "Rtau_vs_L.png")

    print(f"m_pick = {m_pick:.16e} eV")
    print(f"A_unit(m_pick) = {A_unit_vals[idx]:.16e}")
    print(f"A(m_pick; phi_amp={args.phi_amp:g}) = {A_pick:.16e}")
    print(f"max R_tau = {np.nanmax(np.abs(R)):.16e}")
    print(f"L at max R_tau = {L[np.nanargmax(np.abs(R))]}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
