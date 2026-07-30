#!/bin/bash
# Runs the full synthetic-benchmark pipeline end-to-end: dyngen dataset generation,
# then Cyclum, scPrisma, and Oscope, then evaluation.
#
# Each step runs as its own `Rscript` process rather than being sourced inside one
# long-lived R session. This matters specifically for Cyclum/scPrisma: reticulate
# binds an R session to one conda env's Python for the lifetime of that session, so
# there's no way to switch environments mid-session - a fresh process per step avoids
# that entirely, since each one starts with no leftover Python state.
#
# Not included here (separate, optional workflows - run manually if needed):
#   - cyclum_gs.R / scPrisma_gs.R / oscope_gs.R + evaluate_gs.R: hyperparameter
#     grid search over a tuning dataset, not the main cells/genes benchmark sweep.
#   - plotting.R / plot_c1000g2007_1.R: one-off exploratory plots.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

run_step () {
  echo
  echo "========================================================================"
  echo "STEP: $1"
  echo "========================================================================"
  Rscript "$SCRIPT_DIR/$1"
}

run_step "generate_datasets.R"
run_step "run_cyclum.R"
run_step "run_scPrisma.R"
run_step "run_oscope.R"
run_step "evaluate.R"

echo
echo "Pipeline complete."
