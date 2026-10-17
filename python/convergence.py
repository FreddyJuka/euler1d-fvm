#!/usr/bin/env python3
"""Summarize a grid-convergence study of the density-wave test.

Reads every ``summary.csv`` under OUTPUT_ROOT (layout ``<variant>/N<N>/summary.csv``),
computes observed orders of accuracy between consecutive grids for the L1, L2 and Linf
density errors, and writes ``convergence.csv``, ``convergence.md`` and
``convergence_L1.png`` into OUTPUT_ROOT.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

NORMS = ("L1", "L2", "Linf")


def read_runs(root: Path) -> dict[str, list[dict[str, float]]]:
    """Return runs grouped by variant name, sorted by n_cells."""
    groups: dict[str, list[dict[str, float]]] = {}
    for path in sorted(root.rglob("summary.csv")):
        variant = path.relative_to(root).parts[0]
        with path.open(newline="") as f:
            row = next(csv.DictReader(f))
        groups.setdefault(variant, []).append(
            {
                "n_cells": int(row["n_cells"]),
                "dx": float(row["dx"]),
                "L1": float(row["L1_rho"]),
                "L2": float(row["L2_rho"]),
                "Linf": float(row["Linf_rho"]),
            }
        )
    for runs in groups.values():
        runs.sort(key=lambda r: r["n_cells"])
    return groups


def observed_order(e_coarse: float, e_fine: float, n_coarse: int, n_fine: int) -> float:
    """p = log(e_coarse/e_fine) / log(N_fine/N_coarse)."""
    if e_coarse <= 0.0 or e_fine <= 0.0:
        return math.nan
    return math.log(e_coarse / e_fine) / math.log(n_fine / n_coarse)


def add_orders(runs: list[dict[str, float]]) -> None:
    """Attach p_<norm> to each run (NaN for the coarsest grid)."""
    for i, run in enumerate(runs):
        for norm in NORMS:
            if i == 0:
                run[f"p_{norm}"] = math.nan
            else:
                prev = runs[i - 1]
                run[f"p_{norm}"] = observed_order(
                    prev[norm], run[norm], int(prev["n_cells"]), int(run["n_cells"])
                )


def fmt_order(p: float) -> str:
    return "" if math.isnan(p) else f"{p:.2f}"


def write_csv(root: Path, groups: dict[str, list[dict[str, float]]]) -> None:
    header = ["variant", "n_cells", "dx", *NORMS, *(f"p_{n}" for n in NORMS)]
    with (root / "convergence.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for variant, runs in groups.items():
            for r in runs:
                writer.writerow(
                    [variant, int(r["n_cells"]), f"{r['dx']:.16e}"]
                    + [f"{r[n]:.16e}" for n in NORMS]
                    + [fmt_order(r[f"p_{n}"]) for n in NORMS]
                )


def write_markdown(root: Path, groups: dict[str, list[dict[str, float]]]) -> None:
    lines = ["# Density-wave convergence", ""]
    for variant, runs in groups.items():
        lines += [
            f"## {variant}",
            "",
            "| N | dx | L1 | L2 | Linf | p_L1 | p_L2 | p_Linf |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for r in runs:
            lines.append(
                f"| {int(r['n_cells'])} | {r['dx']:.3e} | "
                + " | ".join(f"{r[n]:.3e}" for n in NORMS)
                + " | "
                + " | ".join(fmt_order(r[f"p_{n}"]) or "-" for n in NORMS)
                + " |"
            )
        lines.append("")
    (root / "convergence.md").write_text("\n".join(lines))


def plot_l1(root: Path, groups: dict[str, list[dict[str, float]]]) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for variant, runs in groups.items():
        dx = np.array([r["dx"] for r in runs])
        err = np.array([r["L1"] for r in runs])
        ax.loglog(dx, err, marker="o", label=variant)

    dx_all = [r["dx"] for runs in groups.values() for r in runs]
    # Anchor the reference slopes at the finest grid of the best-resolved variant.
    best = min(groups.values(), key=lambda runs: runs[-1]["L1"])
    dx0, e0 = best[-1]["dx"], best[-1]["L1"]
    span = np.array([min(dx_all), max(dx_all)])
    for order, style in ((1, "--"), (2, ":")):
        ax.loglog(span, e0 * (span / dx0) ** order, "k" + style, linewidth=1,
                  label=f"slope {order}")

    ax.set_xlabel("dx")
    ax.set_ylabel(r"$L_1$ error in density")
    ax.set_title("Density-wave grid convergence")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(root / "convergence_L1.png", dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_root", type=Path, help="directory containing <variant>/N<N>/")
    args = parser.parse_args()

    groups = read_runs(args.output_root)
    if not groups:
        parser.error(f"no summary.csv files found under {args.output_root}")
    for runs in groups.values():
        add_orders(runs)
    write_csv(args.output_root, groups)
    write_markdown(args.output_root, groups)
    plot_l1(args.output_root, groups)
    print(f"Wrote convergence.csv, convergence.md, convergence_L1.png to {args.output_root}")


if __name__ == "__main__":
    main()
