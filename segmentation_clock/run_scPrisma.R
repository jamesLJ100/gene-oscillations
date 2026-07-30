rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)

use_condaenv("scPrisma_env", required = TRUE)
python_exe <- py_config()$python

# mme95 uses the silhouette-filtered cells; mESC/hIPSC don't have a silhouette
# step yet, so they use the raw (QC'd) counts directly.
datasets <- list(
  mme95 = file.path(proj_root, "segmentation_clock/data/mme95/counts/silhouette"),
  mESC  = file.path(proj_root, "segmentation_clock/data/mmESC/counts/raw"),
  hIPSC = file.path(proj_root, "segmentation_clock/data/hIPSC/counts/raw")
)

for (name in names(datasets)) {
  input_dir  <- datasets[[name]]
  output_dir <- file.path(proj_root, "segmentation_clock/results/scPrisma", name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n========================================\n")
  cat("scPrisma:", name, "\n")
  cat("========================================\n")

  args <- c(
    "-u", "python/run_scPrisma.py",
    input_dir,
    "--output_dir", output_dir,
    "--regularisation_strength", "0.1",
    "--iternum", "100"
  )

  exit_code <- system2(
    command = python_exe,
    args    = args,
    stdout  = "",
    stderr  = ""
  )
}
