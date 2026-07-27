#!/usr/bin/env bash
set -e

# Multivariate LAFP CLI Demonstration Script
# Fits a 3-variate time series using `lafp-fit` and exports posterior output files.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

BUILD_CLI="$REPO_DIR/build/lafp-fit"
Y_FILE="$SCRIPT_DIR/y_multivariate_demo.txt"
TOBS_FILE="$SCRIPT_DIR/tobs_demo.txt"
OUT_PREFIX="$SCRIPT_DIR/multivariate_demo_output"

if [ ! -f "$BUILD_CLI" ]; then
    echo "[ERROR] lafp-fit CLI binary not found at $BUILD_CLI."
    echo "Please build the project first: cmake -S . -B build && cmake --build build"
    exit 1
fi

echo "=========================================================================="
echo " Running lafp-fit CLI on 3-variate observation data (Durante et al. 2014)"
echo "=========================================================================="

"$BUILD_CLI" \
    --input-y "$Y_FILE" \
    --input-t "$TOBS_FILE" \
    --niter 500 \
    --n-factors 2 \
    --n-dict 2 \
    --out-prefix "$OUT_PREFIX"

echo ""
echo "Exported posterior output files:"
echo "  - Mean trajectories      : ${OUT_PREFIX}_Mu.txt"
echo "  - Dynamic Covariances    : ${OUT_PREFIX}_Sigma0.txt"
echo "  - Fitted values          : ${OUT_PREFIX}_yhat.txt"
echo "  - Dictionary process draws: ${OUT_PREFIX}_Ksi.txt"
echo "  - Factor process draws    : ${OUT_PREFIX}_Psi.txt"
echo "Done!"
