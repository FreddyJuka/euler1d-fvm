#!/usr/bin/env python3
"""Check that a shock-tube run's L1 density error does not exceed a threshold.

Usage: check_shock_tube.py METRICS_CSV --label TEXT --max-l1 VALUE
Exits with status 1 if the labelled run is missing or L1_rho > VALUE.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("csv", type=Path, help="<case>_metrics.csv from plot_riemann.py")
    parser.add_argument("--label", required=True)
    parser.add_argument("--max-l1", required=True, type=float)
    args = parser.parse_args()

    with args.csv.open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["label"] == args.label]
    if not rows:
        print(f"error: no run labelled '{args.label}' in {args.csv}", file=sys.stderr)
        return 1
    l1 = float(rows[0]["L1_rho"])
    if not l1 <= args.max_l1:
        print(f"error: {args.csv.stem} '{args.label}': L1_rho = {l1:.6e} > {args.max_l1:.6e}",
              file=sys.stderr)
        return 1
    print(f"ok: {args.csv.stem} '{args.label}': L1_rho = {l1:.6e} <= {args.max_l1:.6e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
