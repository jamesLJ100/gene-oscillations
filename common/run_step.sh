# Shared helper for the per-dataset run_pipeline.sh scripts (myc/, nfkb/,
# segmentation_clock/). Runs a step script as its own `Rscript` process rather than
# sourcing it inside one long-lived R session - this matters specifically for
# Cyclum/scPrisma: reticulate binds an R session to one conda env's Python for the
# lifetime of that session, so there's no way to switch environments mid-session; a
# fresh process per step avoids that entirely.
#
# Callers must set SCRIPT_DIR (their own directory) before calling run_step.

run_step () {
  echo
  echo "========================================================================"
  echo "STEP: $1"
  echo "========================================================================"
  Rscript "$SCRIPT_DIR/$1"
}
