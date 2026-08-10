rm(list=ls())
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_oscope.R"))
source(file.path(proj_root, "synthetic/utils/dyngen_utils.R"))

# Run against the held-out evaluation set (not the tuning set gridsearch/oscope_gs.R uses),
# with each (n_cells, n_genes) combination's own best hyperparameters if
# gridsearch/oscope_gs.R has been run for it, otherwise these defaults.
eval_root       <- file.path(proj_root, "synthetic/data/dyngen")
gridsearch_root <- file.path(eval_root, "gridsearch")
combos          <- list_combo_dirs(eval_root)

default_hyperparams <- list(
  KM_maxK           = NULL,
  KM_quan           = 0.95,
  CalcMV_MeanCutLow = 100,
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
    # KM_maxK is NA in grid_search_log.csv when it was fixed to "auto" (NULL) rather
    # than tuned - as.integer(NA) would silently pass NA through to OscopeKM() instead
    # of letting it compute its own default, so map both NULL (no grid search run yet)
    # and NA (grid search ran, left as auto) back to NULL here.
    KM_maxK           = if (is.null(hp$KM_maxK) || is.na(hp$KM_maxK)) NULL else as.integer(hp$KM_maxK),
    KM_quan           = as.numeric(hp$KM_quan),
    CalcMV_MeanCutLow = as.numeric(hp$CalcMV_MeanCutLow),
    FlagCluster_qt    = as.numeric(hp$FlagCluster_qt),
    FlagCluster_thre  = as.numeric(hp$FlagCluster_thre)
  )

  run_oscope(combos$path[i], oscope_config, skip_if_exists = TRUE)
}
