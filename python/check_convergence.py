#!/usr/bin/env python3
"""Check that the observed order of accuracy on the finest grid lies in an interval.

Usage: check_convergence.py CONVERGENCE_CSV --variant NAME --norm {L1,L2,Linf} --min P --max P
Exits with status 1 if the order is outside [min, max] or the variant is missing.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def finest_order(path: Path, variant: str, norm: str) -> float | None:
    """Return p_<norm> on the finest grid of ``variant``, or None if the variant is absent."""
    with path.open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["variant"] == variant]
    if not rows:
        return None
    finest = max(rows, key=lambda r: int(r["n_cells"]))
    value = finest[f"p_{norm}"]
    return float(value) if value else float("nan")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("csv", type=Path, help="convergence.csv")
    parser.add_argument("--variant", required=True)
    parser.add_argument("--norm", required=True, choices=["L1", "L2", "Linf"])
    parser.add_argument("--min", type=float, required=True, dest="lo")
    parser.add_argument("--max", type=float, required=True, dest="hi")
    args = parser.parse_args()

    order = finest_order(args.csv, args.variant, args.norm)
    if order is None:
        print(f"FAIL: variant '{args.variant}' not found in {args.csv}")
        return 1
    ok = args.lo <= order <= args.hi
    print(f"{'PASS' if ok else 'FAIL'}: {args.variant} {args.norm} order = {order:.3f} "
          f"(expected in [{args.lo}, {args.hi}])")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
