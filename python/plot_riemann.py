#!/usr/bin/env python3
"""Compare euler1d shock-tube runs against the exact Riemann solution.

Usage: plot_riemann.py --config CASE.cfg --run DIR:LABEL [--run DIR:LABEL ...] --out FIG.png
"""
from __future__ import annotations

import argparse
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

    for k, spec in enumerate(args.run):
        run_dir, _, label = spec.partition(":")
        label = label or run_dir
        data = np.genfromtxt(Path(run_dir) / "solution_final.csv", delimiter=",", names=True)
        x = data["x"]
        every = max(1, len(x) // 60)
        for ax, key in zip(axes, ["rho", "u", "p", "e"]):
            ax.plot(x, data[key], linestyle="none", marker=MARKERS[k % len(MARKERS)],
                    markersize=4, markevery=every, markerfacecolor="none", label=label)
        rho_ex = exact_on(cfg, x)[0]
        dx = (x_max - x_min) / len(x)
        print(f"{label}: L1(rho) = {np.sum(np.abs(data['rho'] - rho_ex)) * dx:.6e}")

    axes[0].legend()
    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
