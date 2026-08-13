rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_scPrisma.R"))

use_condaenv("scPrisma_env", required = TRUE)

# All three datasets use the neighbor-purity-filtered cells (see common/cell_purity_filter.R).
datasets <- list(
  mme95 = file.path(proj_root, "segmentation_clock/data/mme95/counts/purity_filtered"),
  mESC  = file.path(proj_root, "segmentation_clock/data/mmESC/counts/purity_filtered"),
  hIPSC = file.path(proj_root, "segmentation_clock/data/hIPSC/counts/purity_filtered")
)

for (name in names(datasets)) {
  cat("\n========================================\n")
  cat("scPrisma:", name, "\n")
  cat("========================================\n")

  run_scPrisma(
    input_dir                = datasets[[name]],
    output_dir               = file.path(proj_root, "segmentation_clock/results/scPrisma", name),
    regularisation_strength  = 0.1,
    iternum                  = 100,
    # Raised from the 5000-gene default: checked directly that several segmentation-clock
    # marker genes (Hes7, Dkk1, Axin2) get dropped by HVG selection at 5000 in some clusters
    # but survive at 10000 - not a full fix (Lfng/Dkk1 still drop out in a few clusters even
    # at 10000), but a real improvement. Segmentation-clock-specific: other scPrisma callers
    # (myc/, nfkb/, synthetic/) are untouched and keep the 5000 default.
    n_top_genes              = 10000
  )
}
