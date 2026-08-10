library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_scPrisma.R"))
source(file.path(proj_root, "common/grid_search.R"))
source(file.path(proj_root, "common/utils.R"))
source(file.path(proj_root, "synthetic/utils/dyngen_utils.R"))

use_condaenv("scPrisma_env", required = TRUE)

# Find optimal hyperparameters separately for each (n_cells, n_genes) combination -
# see cyclum_gs.R for the full rationale.
tuning_root <- file.path(proj_root, "synthetic/data/dyngen/gridsearch")
combos      <- list_combo_dirs(tuning_root)

iternum <- 100

# Define grid search parameters
param_grid <- data.frame(
  regularisation_strength = c(0, 0.1, 0.2, 0.3, 0.4)
)

make_run_fn <- function(combo_dir, fnames) {
  force(combo_dir)
  force(fnames)

  function(params, run_number) {
    exit_code <- run_scPrisma(
      input_dir                = combo_dir,
      regularisation_strength  = params$regularisation_strength,
      iternum                  = iternum
    )
    if (exit_code != 0) {
      cat(sprintf("WARNING: Run %d exited with code %d\n", run_number, exit_code))
    }

    aucs <- sapply(fnames, function(fname) {
      results_df <- get_model_scores("scPrisma", combo_dir, fname)
      score_against_ground_truth(fname, combo_dir, results_df)$auc
    })

    list(
      iternum         = iternum,
      exit_code       = exit_code,
      mean_tuning_auc = mean(aucs, na.rm = TRUE),
      n_tuning_files  = sum(!is.na(aucs))
    )
  }
}

best_hyperparams <- tune_per_combo(combos, param_grid, make_run_fn, "scPrisma")
