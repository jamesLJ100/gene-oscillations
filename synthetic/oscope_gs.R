library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_oscope.R"))
source(file.path(proj_root, "common/grid_search.R"))

input_dir <- file.path(proj_root, "synthetic/data/dyngen/tuning")

# Fixed parameter (not part of grid search, but logged alongside the swept ones for context)
fixed_KM_maxK <- 50

# Define grid search parameters
param_grid <- expand.grid(
  KM_quan           = c(0.05),
  FlagCluster_qt    = c(0.90, 0.95, 0.85),
  FlagCluster_thre  = c(pi/5, pi/4, pi/3),
  CalcMV_MeanCutLow = c(0.1, 1, 10),
  stringsAsFactors = FALSE
)

run_fn <- function(params, run_number) {
  oscope_config <- list(
    KM_maxK           = fixed_KM_maxK,
    KM_quan           = params$KM_quan,
    CalcMV_MeanCutLow = params$CalcMV_MeanCutLow,
    FlagCluster_qt    = params$FlagCluster_qt,
    FlagCluster_thre  = params$FlagCluster_thre,
    Normalise         = TRUE
  )
  run_oscope(input_dir, oscope_config, run_number)
  list(KM_maxK = fixed_KM_maxK)
}

output_csv <- file.path(proj_root, "synthetic/data/dyngen/tuning/oscope/grid_search_log.csv")
run_grid_search(param_grid, run_fn, output_csv)
