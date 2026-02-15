#!/usr/bin/env bash
# ===========================================================================
# Production Runner: Study 1 (Non-Gaussian) + Study 2 (Grid Refinement)
# ===========================================================================
#
# Usage:
#   bash ci-benchmark/run-production.sh [study1|study2|both]
#
# Defaults to "both" if no argument given.
#
# Study 1: 150 reps x 12 DGP cells = 1800 pffr fits
#   Estimated time: ~8-12h on 6 cores (73-157s/fit)
#
# Study 2: 120 reps x 3 DGPs x 3 grids = 1080 pffr fits
#   Estimated time: ~18-24h on 6 cores (210-810s/fit, paired design)
#
# Both studies support resume-on-restart via incremental saves.
# ===========================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

MODE="${1:-both}"

LOG_DIR="ci-benchmark/logs"
mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

run_study1() {
    echo "=========================================="
    echo " Study 1: Non-Gaussian Sandwich (full)"
    echo " 150 reps x 12 DGP cells"
    echo " Started: $(date)"
    echo "=========================================="

    Rscript ci-benchmark/sim-study-nongaussian-sandwich.R full \
        2>&1 | tee "${LOG_DIR}/study1_${TIMESTAMP}.log"

    echo ""
    echo "Study 1 complete: $(date)"
    echo "Results in: ci-benchmark/study1-nongaussian/"
}

run_study2() {
    echo "=========================================="
    echo " Study 2: Grid Refinement (main)"
    echo " 120 reps x 3 DGPs x 3 grids (coarse/medium/fine)"
    echo " Started: $(date)"
    echo "=========================================="

    # Verify pilot results exist (needed for grid selection)
    if [ ! -f "ci-benchmark/study2-grid-refinement/pilot_summary.rds" ]; then
        echo "ERROR: Pilot results not found. Run pilot first:"
        echo "  Rscript ci-benchmark/sim-study-grid-refinement.R pilot"
        exit 1
    fi

    Rscript ci-benchmark/sim-study-grid-refinement.R main \
        2>&1 | tee "${LOG_DIR}/study2_${TIMESTAMP}.log"

    echo ""
    echo "Study 2 complete: $(date)"
    echo "Results in: ci-benchmark/study2-grid-refinement/main/"
}

run_analysis() {
    echo "=========================================="
    echo " Combined Analysis"
    echo " Started: $(date)"
    echo "=========================================="

    Rscript ci-benchmark/study-followup-analysis.R \
        2>&1 | tee "${LOG_DIR}/analysis_${TIMESTAMP}.log"

    echo ""
    echo "Analysis complete: $(date)"
}

case "$MODE" in
    study1)
        run_study1
        ;;
    study2)
        run_study2
        ;;
    both)
        run_study1
        echo ""
        echo "=========================================="
        echo " Study 1 done. Starting Study 2..."
        echo "=========================================="
        echo ""
        run_study2
        echo ""
        echo "=========================================="
        echo " Both studies done. Running analysis..."
        echo "=========================================="
        echo ""
        run_analysis
        ;;
    analysis)
        run_analysis
        ;;
    *)
        echo "Usage: $0 [study1|study2|both|analysis]"
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo " All requested tasks complete."
echo " $(date)"
echo "=========================================="
