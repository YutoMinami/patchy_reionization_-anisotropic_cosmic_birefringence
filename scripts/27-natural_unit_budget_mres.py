from __future__ import annotations

import argparse
import csv
import json
from importlib import util
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_module_from_path(name: str, path: Path):
    spec = util.spec_from_file_location(name, path)
    module = util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


MOD15 = load_module_from_path("surrogate15_for_27", ROOT / "scripts" / "15-a_reff_surrogate_bound.py")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--natural-unit-aeff-csv",
        type=Path,
        default=Path("results/26-natural-unit-aeff-mres/natural_unit_aeff_mres_scan.csv"),
    )
    p.add_argument("--L-min", type=int, default=2)
    p.add_argument("--L-max", type=int, default=1000)
    p.add_argument("--A-cb-limit", type=float, default=1.0e-4)
    p.add_argument("--genuine-limit-fraction", type=float, default=0.3)
    p.add_argument("--tau", type=float, default=0.084)
    p.add_argument("--delta-y", type=float, default=19.0)
    p.add_argument("--Rbar", type=float, default=5.0)
    p.add_argument("--siglnR", type=float, default=0.693)
    p.add_argument("--dpeak-fid", type=float, default=2.9e-5)
    p.add_argument("--siglnl-scale", type=float, default=1.0)
    p.add_argument("--z-max", type=float, default=40.0)
    p.add_argument("--n-z", type=int, default=5000)
    p.add_argument("--z-min-a", type=float, default=0.1)
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/27-natural-unit-budget-mres"),
    )
    return p


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, float | str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def claa_limit(L: np.ndarray, A_cb_limit: float) -> np.ndarray:
    return 2.0 * np.pi * A_cb_limit / (L * (L + 1.0))


def d_from_c(L: np.ndarray, C: np.ndarray) -> np.ndarray:
    return L * (L + 1.0) * C / (2.0 * np.pi)


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    nat_rows = read_csv(args.natural_unit_aeff_csv)
    L = np.arange(args.L_min, args.L_max + 1)
    z_grid = np.linspace(0.0, args.z_max, args.n_z)
    _, z_re = MOD15.a_raw_of_tau_delta(args.tau, args.delta_y, z_grid, args.z_min_a)
    chi_re = float(np.interp(z_re, z_grid, MOD15.cumulative_trapezoid_zero(MOD15.C_KM_S / MOD15.H_of_z(z_grid), z_grid)))
    reff = args.Rbar * np.exp(4.0 * args.siglnR**2)
    l_peak = chi_re / reff
    sigma_lnL = args.siglnl_scale * args.siglnR
    d_tau = args.dpeak_fid * MOD15.d_shape_lognormal(L, l_peak, sigma_lnL)
    c_tau = MOD15.c_from_d(L, d_tau)

    cl_lim = claa_limit(L, args.A_cb_limit)
    dl_lim = d_from_c(L, cl_lim)
    cl_phi_bench = args.genuine_limit_fraction * cl_lim
    dl_phi_bench = d_from_c(L, cl_phi_bench)

    scan_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | str]] = []

    grouped: dict[tuple[float, float], list[dict[str, str]]] = {}
    for row in nat_rows:
        key = (float(row["g_GeV_inv"]), float(row["mass_eV"]))
        grouped.setdefault(key, []).append(row)

    for (g, mass), rows in sorted(grouped.items()):
        row0 = rows[0]
        a_tau_unit = float(row0["A_tau_unit_physical"])
        a_tau_eff = float(row0["A_tau_eff_physical"])
        mass_over_mres = float(row0["mass_over_mres"])
        aeff_over_aunit = float(row0["Aeff_over_Aunit"])

        cl_alpha_tau_unit = (a_tau_unit**2) * c_tau
        cl_alpha_tau_eff = (a_tau_eff**2) * c_tau
        dl_alpha_tau_unit = d_from_c(L, cl_alpha_tau_unit)
        dl_alpha_tau_eff = d_from_c(L, cl_alpha_tau_eff)
        cl_total_unit = cl_phi_bench + cl_alpha_tau_unit
        cl_total_eff = cl_phi_bench + cl_alpha_tau_eff
        dl_total_unit = d_from_c(L, cl_total_unit)
        dl_total_eff = d_from_c(L, cl_total_eff)

        frac_unit = dl_alpha_tau_unit / dl_lim
        frac_eff = dl_alpha_tau_eff / dl_lim
        max_idx_unit = int(np.argmax(frac_unit))
        max_idx_eff = int(np.argmax(frac_eff))

        summary_rows.append(
            {
                "mass_eV": mass,
                "mass_over_mres": mass_over_mres,
                "g_GeV_inv": g,
                "Aeff_over_Aunit": aeff_over_aunit,
                "A_tau_unit_physical": a_tau_unit,
                "A_tau_eff_physical": a_tau_eff,
                "max_frac_unit": float(frac_unit[max_idx_unit]),
                "L_peak_frac_unit": int(L[max_idx_unit]),
                "max_frac_eff": float(frac_eff[max_idx_eff]),
                "L_peak_frac_eff": int(L[max_idx_eff]),
            }
        )

        for ell, dtau_val, ctau_val, cl_u, dl_u, cl_e, dl_e, cl_tu, dl_tu, cl_te, dl_te, cl_lim_val, dl_lim_val in zip(
            L,
            d_tau,
            c_tau,
            cl_alpha_tau_unit,
            dl_alpha_tau_unit,
            cl_alpha_tau_eff,
            dl_alpha_tau_eff,
            cl_total_unit,
            dl_total_unit,
            cl_total_eff,
            dl_total_eff,
            cl_lim,
            dl_lim,
            strict=True,
        ):
            scan_rows.append(
                {
                    "mass_eV": mass,
                    "mass_over_mres": mass_over_mres,
                    "g_GeV_inv": g,
                    "Aeff_over_Aunit": aeff_over_aunit,
                    "A_tau_unit_physical": a_tau_unit,
                    "A_tau_eff_physical": a_tau_eff,
                    "L": int(ell),
                    "D_tau_tau": float(dtau_val),
                    "C_tau_tau": float(ctau_val),
                    "D_alphaalpha_tau_unit": float(dl_u),
                    "D_alphaalpha_tau_eff": float(dl_e),
                    "D_alphaalpha_phi_bench": float(dl_phi_bench[ell - args.L_min]),
                    "D_alphaalpha_total_unit": float(dl_tu),
                    "D_alphaalpha_total_eff": float(dl_te),
                    "D_alphaalpha_limit": float(dl_lim_val),
                    "frac_tau_unit_to_limit": float(dl_u / dl_lim_val),
                    "frac_tau_eff_to_limit": float(dl_e / dl_lim_val),
                }
            )

    write_csv(args.output_dir / "natural_unit_budget_mres_scan.csv", list(scan_rows[0].keys()), scan_rows)
    write_csv(args.output_dir / "natural_unit_budget_mres_summary.csv", list(summary_rows[0].keys()), summary_rows)

    best_eff = max(summary_rows, key=lambda r: float(r["max_frac_eff"]))
    lines = [
        "# 27-natural_unit_budget_mres summary",
        "",
        "This run rebuilds the natural-unit birefringence budget near `m_res` using the visibility-weighted `A_tau_eff` from `26`.",
        "",
        f"- source CSV: `{args.natural_unit_aeff_csv}`",
        f"- surrogate template: `tau={args.tau:.3f}`, `Delta y={args.delta_y:.1f}`, `Rbar={args.Rbar:.1f}`, `siglnR={args.siglnR:.3f}`",
        f"- `L_peak = {l_peak:.3f}`",
        f"- genuine benchmark = `{args.genuine_limit_fraction:.2f}` times anisotropic-CB limit",
        "",
        "Most constraining visibility-weighted point in this scan:",
        f"- `g = {float(best_eff['g_GeV_inv']):.2e} GeV^-1`",
        f"- `m/m_res = {float(best_eff['mass_over_mres']):.6f}`",
        f"- `m = {float(best_eff['mass_eV']):.6e} eV`",
        f"- `Aeff/Aunit = {float(best_eff['Aeff_over_Aunit']):.6e}`",
        f"- `A_tau_eff_physical = {float(best_eff['A_tau_eff_physical']):.6e}`",
        f"- `max D_tau_eff / D_limit = {float(best_eff['max_frac_eff']):.6e}` at `L = {int(best_eff['L_peak_frac_eff'])}`",
    ]
    (args.output_dir / "summary.md").write_text("\n".join(lines) + "\n")
    (args.output_dir / "run_config.json").write_text(
        json.dumps({k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}, indent=2)
    )
    print(f"Saved outputs to: {args.output_dir}")


if __name__ == "__main__":
    main()
