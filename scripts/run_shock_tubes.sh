#!/usr/bin/env bash
# Runs the Sod and Lax shock tubes with several schemes and plots them against the
# exact Riemann solution. Run from anywhere; paths are relative to the repository root.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

BUILD_DIR="build"
OUTPUT_ROOT="results/shock_tubes"
EXE="$BUILD_DIR/euler1d"

if [[ ! -x "$EXE" ]]; then
    echo "Building $EXE ..."
    cmake -S . -B "$BUILD_DIR" >/dev/null
    cmake --build "$BUILD_DIR" -j >/dev/null
fi

# name:flux:reconstruction:limiter:label ("-" = no limiter)
VARIANTS=(
    "hllc_first_order:hllc:first_order:-:HLLC first order"
    "rusanov_muscl_minmod:rusanov:muscl:minmod:Rusanov MUSCL minmod"
    "hllc_muscl_minmod:hllc:muscl:minmod:HLLC MUSCL minmod"
    "hllc_muscl_vanleer:hllc:muscl:vanleer:HLLC MUSCL van Leer"
)

for case_name in sod lax; do
    run_args=()
    for spec in "${VARIANTS[@]}"; do
        IFS=: read -r name flux recon limiter label <<<"$spec"
        outdir="$OUTPUT_ROOT/$case_name/$name"
        args=("cases/$case_name.cfg" --quiet --n-cells 200 --output-dir "$outdir"
              --flux "$flux" --reconstruction "$recon")
        [[ "$limiter" != "-" ]] && args+=(--limiter "$limiter")
        "$EXE" "${args[@]}"
        run_args+=(--run "$outdir:$label")
    done
    "${PYTHON:-python3}" python/plot_riemann.py --config "cases/$case_name.cfg" \
        "${run_args[@]}" --out "$OUTPUT_ROOT/${case_name}_comparison.png" \
        --title "${case_name^} shock tube, N = 200"
done
