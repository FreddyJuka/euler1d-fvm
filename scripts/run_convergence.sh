#!/usr/bin/env bash
# Grid-convergence study: runs every scheme variant on a sequence of grids and
# summarizes the observed orders with python/convergence.py.
# Relative paths are resolved from the repository root.
set -euo pipefail

usage() {
    echo "Usage: $0 [-c CONFIG] [-n \"25 50 100 200 400 800\"] [-o OUTPUT_ROOT] [-b BUILD_DIR]" >&2
}

CONFIG="cases/density_wave.cfg"
GRIDS="25 50 100 200 400 800"
OUTPUT_ROOT="results/convergence"
BUILD_DIR="build"

while getopts "c:n:o:b:h" opt; do
    case "$opt" in
        c) CONFIG="$OPTARG" ;;
        n) GRIDS="$OPTARG" ;;
        o) OUTPUT_ROOT="$OPTARG" ;;
        b) BUILD_DIR="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 2 ;;
    esac
done

cd "$(dirname "${BASH_SOURCE[0]}")/.."

EXE="$BUILD_DIR/euler1d"
if [[ ! -x "$EXE" ]]; then
    echo "Building $EXE ..."
    cmake -S . -B "$BUILD_DIR" >/dev/null
    cmake --build "$BUILD_DIR" -j >/dev/null
fi

# name:flux:reconstruction:limiter ("-" = no limiter)
VARIANTS=(
    "hllc_first_order:hllc:first_order:-"
    "rusanov_first_order:rusanov:first_order:-"
    "hllc_muscl_minmod:hllc:muscl:minmod"
    "hllc_muscl_vanleer:hllc:muscl:vanleer"
    "rusanov_muscl_vanleer:rusanov:muscl:vanleer"
)

for spec in "${VARIANTS[@]}"; do
    IFS=: read -r name flux recon limiter <<<"$spec"
    for n in $GRIDS; do
        outdir="$OUTPUT_ROOT/$name/N$n"
        args=("$CONFIG" --quiet --n-cells "$n" --output-dir "$outdir"
              --flux "$flux" --reconstruction "$recon")
        [[ "$limiter" != "-" ]] && args+=(--limiter "$limiter")
        if "$EXE" "${args[@]}"; then
            printf '%-24s N=%-5s OK\n' "$name" "$n"
        else
            printf '%-24s N=%-5s FAILED\n' "$name" "$n"
            echo "error: run failed for variant '$name' with N=$n" >&2
            exit 1
        fi
    done
done

"${PYTHON:-python3}" python/convergence.py "$OUTPUT_ROOT"
