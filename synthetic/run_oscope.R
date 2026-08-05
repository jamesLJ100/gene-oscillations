rm(list=ls())
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_oscope.R"))
source(file.path(proj_root, "synthetic/dyngen_utils.R"))

# Run against the held-out evaluation set (not the tuning set oscope_gs.R uses),
# with each (n_cells, n_genes) combination's own best hyperparameters if
# oscope_gs.R has been run for it, otherwise these defaults.
eval_root       <- file.path(proj_root, "synthetic/data/dyngen_new")
gridsearch_root <- file.path(eval_root, "gridsearch")
combos          <- list_combo_dirs(eval_root)

default_hyperparams <- list(
  KM_maxK           = 50,
  KM_quan           = 0.05,
  CalcMV_MeanCutLow = 0.1,
  FlagCluster_qt    = 0.9,
  FlagCluster_thre  = pi/4
)

for (i in seq_len(nrow(combos))) {
  cat("\n========================================\n")
  cat("Oscope:", combos$combo[i], "\n")
  cat("========================================\n")

  hp <- get_best_hyperparams("oscope", gridsearch_root, combos$n_cells[i], combos$n_genes[i],
                              defaults = default_hyperparams)

  oscope_config <- list(
    KM_maxK           = as.integer(hp$KM_maxK),
    KM_quan           = as.numeric(hp$KM_quan),
    CalcMV_MeanCutLow = as.numeric(hp$CalcMV_MeanCutLow),
    FlagCluster_qt    = as.numeric(hp$FlagCluster_qt),
    FlagCluster_thre  = as.numeric(hp$FlagCluster_thre)
  )

  run_oscope(combos$path[i], oscope_config)
}
