rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_scPrisma.R"))

use_condaenv("scPrisma_env", required = TRUE)

# mme95 uses the silhouette-filtered cells; mESC/hIPSC don't have a silhouette
# step yet, so they use the raw (QC'd) counts directly.
datasets <- list(
  mme95 = file.path(proj_root, "segmentation_clock/data/mme95/counts/silhouette"),
  mESC  = file.path(proj_root, "segmentation_clock/data/mmESC/counts/raw"),
  hIPSC = file.path(proj_root, "segmentation_clock/data/hIPSC/counts/raw")
)

for (name in names(datasets)) {
  cat("\n========================================\n")
  cat("scPrisma:", name, "\n")
  cat("========================================\n")

  run_scPrisma(
    input_dir                = datasets[[name]],
    output_dir               = file.path(proj_root, "segmentation_clock/results/scPrisma", name),
    regularisation_strength  = 0.1,
    iternum                  = 100
  )
}
