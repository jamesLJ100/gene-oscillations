library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_cyclum.R"))
source(file.path(proj_root, "common/grid_search.R"))
source(file.path(proj_root, "common/utils.R"))
source(file.path(proj_root, "synthetic/dyngen_utils.R"))

use_condaenv("cyclum_env", required = TRUE)

#TODO nonlinear_reg & check for others.

# Find optimal hyperparameters separately for each (n_cells, n_genes) combination,
# using only that combination's tuning-set replicates (synthetic/data/dyngen_new/
# gridsearch/c<n_cells>g<n_genes>/, produced by generate_datasets.R with
# is_tuning = TRUE). Results are evaluated against the *other* held-out set
# (dyngen_new/c<n_cells>g<n_genes>/) in evaluate.R, using each combination's
# chosen best hyperparameters via run_cyclum.R - never against the same data the
# hyperparameters were chosen on.
tuning_root <- file.path(proj_root, "synthetic/data/dyngen_new/gridsearch")
combos      <- list_combo_dirs(tuning_root)

# Define grid search parameters. encoder_width is vector-valued, so it's stored as a
# comma-separated string in param_grid (kept a plain flat data frame, parsed back to
# an integer vector inside run_fn) rather than a list-column.
encoder_widths_list <- list(
  c(30, 20),
  c(40, 30),
  c(50, 50),
  c(50, 40, 30),
  c(60, 40, 20)
)
encoder_width_strs <- sapply(encoder_widths_list, paste, collapse = ", ")

param_grid <- expand.grid(
  encoder_width = encoder_width_strs,
  epochs        = c(300, 500, 700),
  learning_rate = c(1e-4, 2e-4, 5e-4, 1e-3),
  stringsAsFactors = FALSE
)

# Builds the run_fn tune_per_combo() calls for each (n_cells, n_genes) combination,
# capturing that combination's own tuning directory/file list for scoring.
make_run_fn <- function(combo_dir, fnames) {
  force(combo_dir)
  force(fnames)

  function(params, run_number) {
    encoder_width <- as.integer(trimws(strsplit(params$encoder_width, ",")[[1]]))

    exit_code <- run_cyclum(
      input_dir     = combo_dir,
      encoder_width = encoder_width,
      epochs        = params$epochs,
      learning_rate = params$learning_rate,
      run_number    = run_number
    )

    # Score against ground truth on each of this combination's tuning replicates;
    # the mean across them is what selects the best hyperparameters for this combo.
    aucs <- sapply(fnames, function(fname) {
      results_df <- get_model_scores("cyclum", combo_dir, fname, run_id = run_number)
      score_against_ground_truth(fname, combo_dir, results_df)$auc
    })

    list(
      exit_code       = exit_code,
      mean_tuning_auc = mean(aucs, na.rm = TRUE),
      n_tuning_files  = sum(!is.na(aucs))
    )
  }
}

best_hyperparams <- tune_per_combo(combos, param_grid, make_run_fn, "cyclum")
