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
# Not included here - run manually first if you want it:
#   - cyclum_gs.R / scPrisma_gs.R / oscope_gs.R: hyperparameter grid search, tuned
#     separately per (n_cells, n_genes) combination against the tuning-set replicates
#     generate_datasets.R produces (dyngen_new/gridsearch/c<n_cells>g<n_genes>/).
#     Optional but recommended - if run BEFORE run_cyclum.R/run_scPrisma.R/
#     run_oscope.R below, each combination's own best hyperparameters are picked up
#     automatically (from dyngen_new/gridsearch/best_hyperparams_<algorithm>.csv);
#     otherwise those steps fall back to fixed defaults. Not wired into this script
#     by default since the full grid search is substantially more compute than the
#     rest of the pipeline combined.
#
# Also not included:
#   - evaluate_gs.R: now superseded by tune_per_combo() (called from within
#     cyclum_gs.R/scPrisma_gs.R/oscope_gs.R), which computes AUC and picks the best
#     hyperparameters per combination automatically as part of the grid search
#     itself - evaluate_gs.R's separate load-and-rescan-the-log approach predates
#     that and no longer matches the per-combination directory layout.
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
