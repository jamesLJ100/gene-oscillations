#!/bin/bash
# Runs the full segmentation_clock pipeline end-to-end: preprocessing for all three
# datasets, then Cyclum and scPrisma, then evaluation.
#
# Each step runs as its own `Rscript` process rather than being sourced inside one
# long-lived R session. This matters specifically for Cyclum/scPrisma: reticulate
# binds an R session to one conda env's Python for the lifetime of that session, so
# there's no way to switch environments mid-session - a fresh process per step avoids
# that entirely, since each one starts with no leftover Python state.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../common/run_step.sh"

run_step "preprocess_mme95.R"
run_step "preprocess_mesc.R"
run_step "preprocess_hipsc.R"
run_step "run_cyclum.R"
run_step "run_scPrisma.R"
run_step "analyse.R"

echo
echo "Pipeline complete."
