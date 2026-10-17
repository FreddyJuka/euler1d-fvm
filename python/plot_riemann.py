#!/usr/bin/env python3
"""Compare euler1d shock-tube runs against the exact Riemann solution.

Usage: plot_riemann.py --config CASE.cfg --run DIR:LABEL [--run DIR:LABEL ...] --out FIG.png
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from exact_riemann import solve  # noqa: E402
from riemann_config import load_config  # noqa: E402

MARKERS = ["o", "s", "^", "v", "D", "x", "+"]


def exact_on(cfg: dict, x: np.ndarray) -> tuple[np.ndarray, ...]:
    """Exact (rho, u, p, e) at points x for the case described by cfg."""
    g = float(cfg["gamma"])
    WL = (float(cfg["rho_L"]), float(cfg["u_L"]), float(cfg["p_L"]))
    WR = (float(cfg["rho_R"]), float(cfg["u_R"]), float(cfg["p_R"]))
    rho, u, p = solve(WL, WR, g, x, float(cfg["t_final"]), float(cfg["x0"]))
    return rho, u, p, p / ((g - 1.0) * rho)


def check_complete(run_dir: Path, row: dict[str, str]) -> None:
    """Exit with an error if the run in ``run_dir`` stopped before t_final."""
    t_final, t_reached = float(row["t_final"]), float(row["t_reached"])
    if abs(t_reached - t_final) > 1e-12 * abs(t_final):
        sys.exit(
            f"error: run in {run_dir} is truncated: t_reached = {t_reached!r} "
            f"differs from t_final = {t_final!r}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--run", action="append", required=True, metavar="DIR:LABEL")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    x_min, x_max = float(cfg["x_min"]), float(cfg["x_max"])
    x_ex = np.linspace(x_min, x_max, 2000)
    exact = exact_on(cfg, x_ex)

    labels = ["density", "velocity", "pressure", "specific internal energy"]
    symbols = [r"$\rho$", r"$u$", r"$p$", r"$e$"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    axes = axes.ravel()
    for ax, y, name, sym in zip(axes, exact, labels, symbols):
        ax.plot(x_ex, y, "k-", linewidth=1.2, label="exact")
        ax.set_ylabel(f"{name} {sym}")
        ax.grid(True, alpha=0.3)
    for ax in axes[2:]:
        ax.set_xlabel("x")

    metrics = []
    for k, spec in enumerate(args.run):
        run_dir, _, label = spec.partition(":")
        label = label or run_dir
        with (Path(run_dir) / "summary.csv").open(newline="") as f:
            check_complete(Path(run_dir), next(csv.DictReader(f)))
        data = np.genfromtxt(Path(run_dir) / "solution_final.csv", delimiter=",", names=True)
        x = data["x"]
        every = max(1, len(x) // 60)
        for ax, key in zip(axes, ["rho", "u", "p", "e"]):
            ax.plot(x, data[key], linestyle="none", marker=MARKERS[k % len(MARKERS)],
                    markersize=4, markevery=every, markerfacecolor="none", label=label)
        rho_ex = exact_on(cfg, x)[0]
        dx = (x_max - x_min) / len(x)
        l1 = np.sum(np.abs(data["rho"] - rho_ex)) * dx
        print(f"{args.config.stem} | {label} | L1_rho = {l1:.6e}")
        i_max = int(np.argmax(data["e"]))
        metrics.append([args.config.stem, label, len(x), f"{l1:.6e}", f"{data['e'][i_max]:.6e}",
                        f"{x[i_max]:.6e}", f"{exact[3].max():.6e}"])

    axes[0].legend()
    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    plt.close(fig)

    metrics_path = args.out.parent / f"{args.config.stem}_metrics.csv"
    with metrics_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["case", "label", "n_cells", "L1_rho", "e_max_numerical",
                         "x_at_e_max", "e_max_exact"])
        writer.writerows(metrics)


if __name__ == "__main__":
    main()
