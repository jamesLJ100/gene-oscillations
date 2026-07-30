rm(list=ls())
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_oscope.R"))

oscope_config <- list(
  KM_maxK           = 50,
  KM_quan           = 0.05,
  CalcMV_MeanCutLow = 0.1,
  FlagCluster_qt    = 0.9,
  FlagCluster_thre  = pi/4,
  Normalise         = TRUE
)

run_oscope(file.path(proj_root, "synthetic/data/dyngen_new"), oscope_config)
