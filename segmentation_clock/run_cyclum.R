rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)

use_condaenv("cyclum_env", required = TRUE)
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
  output_dir <- file.path(proj_root, "segmentation_clock/results/cyclum", name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n========================================\n")
  cat("Cyclum:", name, "\n")
  cat("========================================\n")

  args <- c(
    "-u", "python/run_cyclum.py",
    input_dir,
    "--output_dir", output_dir,
    "--encoder_width", "30", "20",
    "--epochs", "500",
    "--learning_rate", "2e-4"
  )

  exit_code <- system2(
    command = python_exe,
    args    = args,
    stdout  = "",
    stderr  = ""
  )
}
