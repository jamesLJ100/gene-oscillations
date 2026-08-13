rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_cyclum.R"))

use_condaenv("cyclum_env", required = TRUE)

# All three datasets use the neighbor-purity-filtered cells (see common/cell_purity_filter.R).
datasets <- list(
  mme95 = file.path(proj_root, "segmentation_clock/data/mme95/counts/purity_filtered"),
  mESC  = file.path(proj_root, "segmentation_clock/data/mmESC/counts/purity_filtered"),
  hIPSC = file.path(proj_root, "segmentation_clock/data/hIPSC/counts/purity_filtered")
)

for (name in names(datasets)) {
  cat("\n========================================\n")
  cat("Cyclum:", name, "\n")
  cat("========================================\n")

  run_cyclum(
    input_dir     = datasets[[name]],
    output_dir    = file.path(proj_root, "segmentation_clock/results/cyclum", name),
    encoder_width = c(30, 20),
    epochs        = 500,
    learning_rate = 2e-4
  )
}
