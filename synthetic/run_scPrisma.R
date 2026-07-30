rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_scPrisma.R"))

use_condaenv("scPrisma_env", required = TRUE)

run_scPrisma(
  input_dir                = file.path(proj_root, "synthetic/data/dyngen_new"),
  regularisation_strength  = 0.1,
  iternum                  = 100
)

# Real datasets: see myc/run_scPrisma.R, nfkb/run_scPrisma.R, segmentation_clock/run_scPrisma.R
