rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_scPrisma.R"))
source(file.path(proj_root, "synthetic/dyngen_utils.R"))

use_condaenv("scPrisma_env", required = TRUE)

# Run against the held-out evaluation set (not the tuning set scPrisma_gs.R uses),
# with each (n_cells, n_genes) combination's own best hyperparameters if
# scPrisma_gs.R has been run for it, otherwise these defaults.
eval_root       <- file.path(proj_root, "synthetic/data/dyngen_new")
gridsearch_root <- file.path(eval_root, "gridsearch")
combos          <- list_combo_dirs(eval_root)

default_hyperparams <- list(regularisation_strength = 0.1, iternum = 100)

for (i in seq_len(nrow(combos))) {
  cat("\n========================================\n")
  cat("scPrisma:", combos$combo[i], "\n")
  cat("========================================\n")

  hp <- get_best_hyperparams("scPrisma", gridsearch_root, combos$n_cells[i], combos$n_genes[i],
                              defaults = default_hyperparams)

  run_scPrisma(
    input_dir                = combos$path[i],
    regularisation_strength  = as.numeric(hp$regularisation_strength),
    iternum                  = as.integer(hp$iternum)
  )
}

# Real datasets: see myc/run_scPrisma.R, nfkb/run_scPrisma.R, segmentation_clock/run_scPrisma.R
