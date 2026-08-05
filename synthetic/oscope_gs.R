library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_oscope.R"))
source(file.path(proj_root, "common/grid_search.R"))
source(file.path(proj_root, "common/utils.R"))
source(file.path(proj_root, "synthetic/dyngen_utils.R"))

# Find optimal hyperparameters separately for each (n_cells, n_genes) combination -
# see cyclum_gs.R for the full rationale.
tuning_root <- file.path(proj_root, "synthetic/data/dyngen_new/gridsearch")
combos      <- list_combo_dirs(tuning_root)

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

make_run_fn <- function(combo_dir, fnames) {
  force(combo_dir)
  force(fnames)

  function(params, run_number) {
    oscope_config <- list(
      KM_maxK           = fixed_KM_maxK,
      KM_quan           = params$KM_quan,
      CalcMV_MeanCutLow = params$CalcMV_MeanCutLow,
      FlagCluster_qt    = params$FlagCluster_qt,
      FlagCluster_thre  = params$FlagCluster_thre
    )
    run_oscope(combo_dir, oscope_config, run_number)

    aucs <- sapply(fnames, function(fname) {
      results_df <- get_model_scores("oscope", combo_dir, fname, run_id = run_number)
      score_against_ground_truth(fname, combo_dir, results_df)$auc
    })

    list(
      KM_maxK         = fixed_KM_maxK,
      mean_tuning_auc = mean(aucs, na.rm = TRUE),
      n_tuning_files  = sum(!is.na(aucs))
    )
  }
}

best_hyperparams <- tune_per_combo(combos, param_grid, make_run_fn, "oscope")
