rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)

use_condaenv("scPrisma_env", required = TRUE)

dyngen_dir <- file.path(proj_root, "synthetic/data/dyngen_new")


python_exe <- py_config()$python

reg_str <- 0.1
iter <- 100

args <- c(
  "-u",
  "python/run_scPrisma.py",
  dyngen_dir,
  "--regularisation_strength", as.character(reg_str),
  "--iternum", as.character(iter)
)

exit_code <- system2(
  command = python_exe,
  args = args,
  stdout = "",
  stderr = ""
)

# Real datasets: see myc/run_scPrisma.R, nfkb/run_scPrisma.R, segmentation_clock/run_scPrisma.R

