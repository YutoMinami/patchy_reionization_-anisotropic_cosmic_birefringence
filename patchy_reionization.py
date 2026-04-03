from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


# SI constants
HBAR = 1.054571817e-34
EV_TO_J = 1.602176634e-19
MPC_TO_M = 3.085677581491367e22
SIGMA_T = 6.6524587321e-29
G_NEWTON = 6.67430e-11
M_P = 1.67262192369e-27


# Planck-like cosmology
H0_KM_S_MPC = 67.4
H0_SI = H0_KM_S_MPC * 1000.0 / MPC_TO_M
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
OMEGA_L = 1.0 - OMEGA_M - OMEGA_R
OMEGA_B = 0.049


# Reionization
Z_REI = 7.7
X_REI = np.log(1.0 / (1.0 + Z_REI))
A_REI = np.exp(X_REI)
Y_HE = 0.24
F_H = 1.0 - Y_HE
X_E_REI = 1.0


RHO_CRIT0 = 3.0 * H0_SI**2 / (8.0 * np.pi * G_NEWTON)
RHO_B0 = OMEGA_B * RHO_CRIT0


@dataclass
class SolveConfig:
    z_ini: float = 1.0e7
    z_end: float = 0.0
    method: str = "DOP853"
    rtol: float = 1.0e-9
    atol_rel: float = 1.0e-12
    dense_output: bool = True
    n_eval: int | None = None


def e2_of_a(a: float | np.ndarray) -> float | np.ndarray:
    return OMEGA_R * a**-4 + OMEGA_M * a**-3 + OMEGA_L


def H_of_a(a: float | np.ndarray) -> float | np.ndarray:
    return H0_SI * np.sqrt(e2_of_a(a))


def dlnH_dx(a: float | np.ndarray) -> float | np.ndarray:
    dE2_dx = -4.0 * OMEGA_R * a**-4 - 3.0 * OMEGA_M * a**-3
    return 0.5 * dE2_dx / e2_of_a(a)


def n_e_of_a(a: float, x_e: float = X_E_REI) -> float:
    rho_b = RHO_B0 * a**-3
    n_b = rho_b / M_P
    n_H = F_H * n_b
    return x_e * n_H


def d_eta_d_tau(a: float, x_e: float = X_E_REI) -> float:
    return 1.0 / (a * n_e_of_a(a, x_e=x_e) * SIGMA_T)


def m_eV_to_omega_SI(m_eV: float) -> float:
    return (m_eV * EV_TO_J) / HBAR


def alp_ode_x(x: float, y: np.ndarray, m_eV: float) -> list[float]:
    a = np.exp(x)
    H = H_of_a(a)
    phi, dphi_dx = y
    mass_ratio = m_eV_to_omega_SI(m_eV) / H
    d2phi_dx2 = -(3.0 + dlnH_dx(a)) * dphi_dx - mass_ratio**2 * phi
    return [dphi_dx, d2phi_dx2]


def solve_phi_background(
    m_eV: float,
    phi_ini: float = 1.0,
    dphi_dx_ini: float = 0.0,
    config: SolveConfig | None = None,
):
    config = config or SolveConfig()
    x_ini = np.log(1.0 / (1.0 + config.z_ini))
    x_end = np.log(1.0 / (1.0 + config.z_end))

    kwargs = {
        "fun": lambda x, y: alp_ode_x(x, y, m_eV),
        "t_span": (x_ini, x_end),
        "y0": [phi_ini, dphi_dx_ini],
        "rtol": config.rtol,
        "atol": [
            config.atol_rel * max(abs(phi_ini), 1.0e-30),
            config.atol_rel * max(abs(phi_ini), abs(dphi_dx_ini), 1.0e-30),
        ],
        "method": config.method,
        "dense_output": config.dense_output,
    }

    if config.n_eval is not None:
        kwargs["t_eval"] = np.linspace(x_ini, x_end, config.n_eval)

    return solve_ivp(**kwargs)


def evaluate_state_at_x(sol, x_target: float = X_REI) -> np.ndarray:
    if sol.sol is not None:
        return np.asarray(sol.sol(x_target))

    y0 = np.interp(x_target, sol.t, sol.y[0])
    y1 = np.interp(x_target, sol.t, sol.y[1])
    return np.array([y0, y1])


def compute_A_from_solution(sol, x_target: float = X_REI) -> float:
    _, dphi_dx = evaluate_state_at_x(sol, x_target=x_target)
    a_target = np.exp(x_target)
    dotphi_phys = H_of_a(a_target) * dphi_dx
    dotphi_conf = a_target * dotphi_phys
    return dotphi_conf * d_eta_d_tau(a_target)


def compute_A_of_m(
    m_eV: float,
    phi_ini: float = 1.0,
    dphi_dx_ini: float = 0.0,
    config: SolveConfig | None = None,
):
    sol = solve_phi_background(
        m_eV=m_eV,
        phi_ini=phi_ini,
        dphi_dx_ini=dphi_dx_ini,
        config=config,
    )
    return compute_A_from_solution(sol)


def compute_A_unit(m_eV: float, config: SolveConfig | None = None) -> float:
    return compute_A_of_m(m_eV=m_eV, phi_ini=1.0, dphi_dx_ini=0.0, config=config)


def rescale_A_from_unit(A_unit: float, phi_amp: float) -> float:
    return A_unit * phi_amp


def Cl_tautau_toy(
    L: np.ndarray,
    Ctau0: float = 1.0e-10,
    Lc: float = 300.0,
    sigmaL: float = 100.0,
) -> np.ndarray:
    return Ctau0 * np.exp(-0.5 * ((L - Lc) / sigmaL) ** 2)


def Cl_phiphi_toy(
    L: np.ndarray,
    Cphi0: float = 1.0e-10,
    Lstar: float = 300.0,
    nphi: float = 2.0,
    Lcut: float | None = None,
) -> np.ndarray:
    Cl = Cphi0 * (L / Lstar) ** (-nphi)
    if Lcut is not None:
        Cl *= np.exp(-(L / Lcut) ** 2)
    return Cl


def R_tau_toy(
    L: np.ndarray,
    A_dimless: float,
    Ctau0: float = 1.0e-10,
    Lc: float = 300.0,
    sigmaL: float = 100.0,
    Cphi0: float = 1.0e-10,
    Lstar: float = 300.0,
    nphi: float = 2.0,
    Lcut: float | None = None,
) -> np.ndarray:
    Cl_tau = Cl_tautau_toy(L, Ctau0=Ctau0, Lc=Lc, sigmaL=sigmaL)
    Cl_phi = Cl_phiphi_toy(L, Cphi0=Cphi0, Lstar=Lstar, nphi=nphi, Lcut=Lcut)
    return (A_dimless**2) * (Cl_tau / Cl_phi)
