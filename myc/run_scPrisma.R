rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)

use_condaenv("scPrisma_env", required = TRUE)

input_dir  <- file.path(proj_root, "myc/data/GSM4286760")
output_dir <- file.path(proj_root, "myc/results/scPrisma")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

python_exe <- py_config()$python

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
