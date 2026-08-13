rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_cyclum.R"))

use_condaenv("cyclum_env", required = TRUE)

run_cyclum(
  input_dir  = file.path(proj_root, "nfkb/data/GSE162992/processed"),
  output_dir = file.path(proj_root, "nfkb/results/cyclum"),
  encoder_width = c(30, 20),
  epochs        = 500,
  learning_rate = 2e-4
)
