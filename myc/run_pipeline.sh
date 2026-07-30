#!/bin/bash
# Runs the full myc (HCT116 MYC oscillations) pipeline end-to-end: preprocessing
# (auto-downloading the raw GEO data if missing), then Cyclum and scPrisma, then
# evaluation.
#
# Each step runs as its own `Rscript` process rather than being sourced inside one
# long-lived R session. This matters specifically for Cyclum/scPrisma: reticulate
# binds an R session to one conda env's Python for the lifetime of that session, so
# there's no way to switch environments mid-session - a fresh process per step avoids
# that entirely, since each one starts with no leftover Python state.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

run_step () {
  echo
  echo "========================================================================"
  echo "STEP: $1"
  echo "========================================================================"
  Rscript "$SCRIPT_DIR/$1"
}

run_step "preprocess_hct116.R"
run_step "run_cyclum.R"
run_step "run_scPrisma.R"
run_step "analyse.R"

echo
echo "Pipeline complete."
