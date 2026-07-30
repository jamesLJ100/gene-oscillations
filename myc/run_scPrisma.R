rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_scPrisma.R"))

use_condaenv("scPrisma_env", required = TRUE)

run_scPrisma(
  input_dir                = file.path(proj_root, "myc/data/GSM4286760"),
  output_dir               = file.path(proj_root, "myc/results/scPrisma"),
  regularisation_strength  = 0.1,
  iternum                  = 100
)

