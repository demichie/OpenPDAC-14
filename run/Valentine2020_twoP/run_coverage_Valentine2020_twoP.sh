#!/bin/sh
set -e

# =============================================================================
# Coverage workflow for OpenPDAC-13/run/Valentine2020_twoP
# -----------------------------------------------------------------------------
# Purpose
#   This script runs a lightweight line-coverage workflow for the
#   Valentine2020_twoP case.
#
# Strategy
#   - Run the case fully in serial.
#   - Perform the standard initialization run using the .init dictionaries.
#   - Capture coverage after the initialization run.
#   - Reset counters.
#   - Perform a short main run using controlDict.coverage instead of
#     controlDict.run.
#   - Capture coverage after the short run.
#   - Merge the two coverage tracefiles and generate an HTML report.
#
# Notes
#   - This script intentionally skips any extra post-processing.
#   - Coverage capture is focused on applications/OpenPDAC.
# =============================================================================

cd "${0%/*}" || exit 1

COVERAGE_ROOT="${COVERAGE_ROOT:-$(pwd)/coverage}"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd ../.. && pwd)}"
COVERAGE_DIR="${COVERAGE_DIR:-$PROJECT_ROOT/applications/OpenPDAC}"

mkdir -p "$COVERAGE_ROOT"

echo "==> Project root:  $PROJECT_ROOT"
echo "==> Coverage root: $COVERAGE_ROOT"
echo "==> Coverage dir:  $COVERAGE_DIR"

print_coverage_debug() {
    label="$1"

    echo "============================================================"
    echo "Coverage debug: $label"
    echo "Coverage dir: $COVERAGE_DIR"
    echo "------------------------------------------------------------"

    echo "[gcno] First matches:"
    find "$COVERAGE_DIR" -name "*.gcno" | head -50 || true
    echo "[gcno] Count:"
    find "$COVERAGE_DIR" -name "*.gcno" | wc -l || true

    echo "------------------------------------------------------------"
    echo "[gcda] First matches:"
    find "$COVERAGE_DIR" -name "*.gcda" | head -50 || true
    echo "[gcda] Count:"
    find "$COVERAGE_DIR" -name "*.gcda" | wc -l || true

    echo "============================================================"
}

require_gcda_files() {
    label="$1"
    count=$(find "$COVERAGE_DIR" -name "*.gcda" | wc -l)

    if [ "$count" -eq 0 ]; then
        echo "ERROR: No .gcda files found in $COVERAGE_DIR after $label"
        echo "This usually means that either:"
        echo "  - the instrumented code was not actually executed, or"
        echo "  - runtime coverage files were written somewhere else."
        exit 1
    fi
}

. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

print_coverage_debug "before clean"

echo "==> Cleaning the case"
chmod +x ./Allclean
./Allclean

echo "==> Preparing controlDict for blockMesh"
cp ./system/controlDict.init ./system/controlDict

echo "==> Running blockMesh"
runApplication blockMesh

echo "==> Running checkMesh"
runApplication checkMesh -allTopology -allGeometry

echo "==> Applying changeDictionary"
runApplication changeDictionary

echo "==> Capturing baseline coverage"
lcov --capture --initial --directory "$COVERAGE_DIR" --output-file "$COVERAGE_ROOT/coverage_base.info"

print_coverage_debug "after baseline capture"

echo "==> Preparing initialization run"
cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution

rm -rf 0
cp -r org.0 0

echo "==> Setting fields for the '.init' run"
for field in alpha.air T.air U.air; do
    mv "0/${field}.init" "0/${field}"
done
for particle in particles1 particles2; do
    for field in alpha T U; do
        mv "0/${field}.${particle}.init" "0/${field}.${particle}"
    done
done

echo "==> Running initialization with $(getApplication)"
runApplication $(getApplication)
mv log.foamRun log.foamRun0 || true

print_coverage_debug "after initialization run"
require_gcda_files "initialization run"

echo "==> Capturing coverage after initialization run"
lcov --capture --directory "$COVERAGE_DIR" --output-file "$COVERAGE_ROOT/coverage_init.info"

echo "==> Zeroing counters before short main run"
lcov --zerocounters --directory "$COVERAGE_DIR"

print_coverage_debug "after zeroing counters"

echo "==> Preparing short main run"
cp ./system/controlDict.coverage ./system/controlDict
cp ./system/fvSolution.run ./system/fvSolution

echo "==> Setting fields for the '.coverage' run"
for field in alpha.air T.air U.air; do
    mv "0/${field}.run" "0/${field}"
done
for particle in particles1 particles2; do
    for field in alpha T U; do
        mv "0/${field}.${particle}.run" "0/${field}.${particle}"
    done
done

echo "==> Running short main simulation with $(getApplication)"
runApplication $(getApplication)

print_coverage_debug "after short main run"
require_gcda_files "short main run"

echo "==> Capturing coverage after short main run"
lcov --capture --directory "$COVERAGE_DIR" --output-file "$COVERAGE_ROOT/coverage_sim.info"

echo "==> Merging tracefiles"
lcov \
  -a "$COVERAGE_ROOT/coverage_base.info" \
  -a "$COVERAGE_ROOT/coverage_init.info" \
  -a "$COVERAGE_ROOT/coverage_sim.info" \
  -o "$COVERAGE_ROOT/coverage_total.info"

echo "==> Filtering external/OpenFOAM files"
lcov \
  --remove "$COVERAGE_ROOT/coverage_total.info" \
  '/opt/openfoam*' \
  '/usr/*' \
  '*/lnInclude/*' \
  --output-file "$COVERAGE_ROOT/coverage_clean.info"

echo "==> Generating HTML report"
genhtml "$COVERAGE_ROOT/coverage_clean.info" --output-directory "$COVERAGE_ROOT/html"

echo "==> Coverage run completed successfully"
echo "HTML report: $COVERAGE_ROOT/html/index.html"
