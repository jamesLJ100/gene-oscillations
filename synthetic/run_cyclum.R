library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_cyclum.R"))
source(file.path(proj_root, "synthetic/utils/dyngen_utils.R"))

use_condaenv("cyclum_env", required = TRUE)

# Run against the held-out evaluation set (not the tuning set gridsearch/cyclum_gs.R uses),
# with each (n_cells, n_genes) combination's own best hyperparameters if
# gridsearch/cyclum_gs.R has been run for it, otherwise these defaults.
eval_root       <- file.path(proj_root, "synthetic/data/dyngen")
gridsearch_root <- file.path(eval_root, "gridsearch")
combos          <- list_combo_dirs(eval_root)

default_hyperparams <- list(encoder_width = c(30, 20), epochs = 500, learning_rate = 2e-4)

for (i in seq_len(nrow(combos))) {
  cat("\n========================================\n")
  cat("Cyclum:", combos$combo[i], "\n")
  cat("========================================\n")

  hp <- get_best_hyperparams("cyclum", gridsearch_root, combos$n_cells[i], combos$n_genes[i],
                              defaults = default_hyperparams)

  # encoder_width comes back as a "40, 30"-style string when sourced from a grid
  # search result (see gridsearch/cyclum_gs.R), or as-is (a numeric vector) from defaults.
  encoder_width <- if (is.character(hp$encoder_width)) {
    as.integer(trimws(strsplit(hp$encoder_width, ",")[[1]]))
  } else {
    as.integer(hp$encoder_width)
  }

  run_cyclum(
    input_dir     = combos$path[i],
    encoder_width = encoder_width,
    epochs        = as.integer(hp$epochs),
    learning_rate = as.numeric(hp$learning_rate),
    skip_if_exists = TRUE
  )
}
