from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
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

from patchy_reionization import SolveConfig, compute_A_unit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m-min", type=float, default=1.0e-35)
    parser.add_argument("--m-max", type=float, default=1.0e-26)
    parser.add_argument("--num-mass", type=int, default=40)
    parser.add_argument("--z-ini", type=float, default=1.0e7)
    parser.add_argument("--output-dir", type=Path, default=Path("results/04a-a-unit"))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    return parser


def append_csv_row(path: Path, fieldnames: list[str], row: dict) -> None:
    needs_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if needs_header:
            writer.writeheader()
        writer.writerow(row)


def read_existing_csv(path: Path) -> dict[float, float]:
    if not path.exists():
        return {}
    out: dict[float, float] = {}
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            out[float(row["mass_eV"])] = float(row["A_unit"])
    return out


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    config = SolveConfig(z_ini=args.z_ini, method="DOP853", rtol=1.0e-9, atol_rel=1.0e-12, dense_output=True)
    m_grid = np.logspace(np.log10(args.m_min), np.log10(args.m_max), args.num_mass)

    csv_path = args.output_dir / "A_unit_scan.csv"
    if csv_path.exists() and not args.resume:
        csv_path.unlink()

    existing = read_existing_csv(csv_path) if args.resume else {}
    A_vals = []
    for idx, mass in enumerate(m_grid, start=1):
        mass_key = float(mass)
        if mass_key in existing:
            A_unit = existing[mass_key]
            A_vals.append(A_unit)
            print(f"[A_unit] {idx}/{len(m_grid)}: m={mass:.6e} eV, A_unit={A_unit:.6e} [resume]")
            continue

        A_unit = compute_A_unit(m_eV=mass, config=config)
        A_vals.append(A_unit)
        append_csv_row(csv_path, ["mass_eV", "A_unit"], {"mass_eV": mass_key, "A_unit": float(A_unit)})
        print(f"[A_unit] {idx}/{len(m_grid)}: m={mass:.6e} eV, A_unit={A_unit:.6e}")

    plt.figure(figsize=(7, 5))
    plt.loglog(m_grid, np.abs(np.array(A_vals)))
    plt.xlabel(r"$m_a \ [{\rm eV}]$")
    plt.ylabel(r"$|A_{\rm unit}(m_a)|$")
    plt.title("Unit-amplitude response scan")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output_dir / "A_unit_vs_m.png", dpi=160, bbox_inches="tight")
    plt.close()

    payload = {
        "solve_config": asdict(config),
        "run_args": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
    }
    (args.output_dir / "run_config.json").write_text(json.dumps(payload, indent=2))
    (args.output_dir / "summary.md").write_text(
        "\n".join(
            [
                "# 04a-scan_a_unit summary",
                "",
                f"- mass range: `{args.m_min:.3e}` to `{args.m_max:.3e}` eV",
                f"- number of masses: `{args.num_mass}`",
                f"- z_ini: `{args.z_ini:.3e}`",
                f"- output csv: `{csv_path.name}`",
            ]
        )
        + "\n"
    )
    print(f"Saved outputs to: {args.output_dir}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
