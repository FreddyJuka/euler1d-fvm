#!/usr/bin/env bash
# Regenerates the figures shown in the README into docs/figures/.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

scripts/run_convergence.sh
scripts/run_shock_tubes.sh

mkdir -p docs/figures
cp results/convergence/convergence_L1.png docs/figures/
cp results/shock_tubes/sod_comparison.png docs/figures/
cp results/shock_tubes/lax_comparison.png docs/figures/
