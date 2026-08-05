#' Run a hyperparameter grid search, independent of which algorithm is being tuned
#'
#' Iterates over `param_grid` (one row per parameter combination), calling
#' `run_fn(params, run_number)` for each. `run_fn` is responsible for actually
#' invoking the algorithm (e.g. via run_cyclum()/run_scPrisma()/run_oscope() from
#' algorithms/) and returning a named list/single-row data frame of whatever's
#' worth logging beyond the parameters themselves (e.g. exit_code). Errors inside
#' `run_fn` are caught per-run - one bad combination doesn't abort the rest of the
#' grid search, and the error message is recorded in an `error` column instead.
#'
#' @param param_grid Data frame, one row per parameter combination to try (e.g. from
#'   expand.grid()). Every column becomes a column in the saved run log.
#' @param run_fn Function(params_row, run_number) -> named list / single-row data
#'   frame of extra columns to log (return list() if there's nothing to add)
#' @param output_csv Path to save the run log CSV to (parent directory created if missing)
#' @return The full run log as a data frame (also written to output_csv)
run_grid_search <- function(param_grid, run_fn, output_csv) {

  library(dplyr)

  total_runs <- nrow(param_grid)
  cat(sprintf("Starting grid search with %d total combinations\n", total_runs))

  run_log <- vector("list", total_runs)

  for (i in seq_len(total_runs)) {
    params <- param_grid[i, , drop = FALSE]

    cat(sprintf("\n========================================\n"))
    cat(sprintf("Run %d/%d\n", i, total_runs))
    print(params)
    cat(sprintf("========================================\n\n"))

    result <- tryCatch(
      run_fn(params, i),
      error = function(e) {
        cat(sprintf("ERROR in run %d: %s\n", i, conditionMessage(e)))
        list(error = conditionMessage(e))
      }
    )

    row <- cbind(data.frame(run_counter = i), params)
    if (length(result) > 0) row <- cbind(row, as.data.frame(result, stringsAsFactors = FALSE))
    rownames(row) <- NULL
    run_log[[i]] <- row
  }

  # bind_rows() (not rbind()) since a run that errors returns different columns
  # (just `error`) than one that succeeds (whatever run_fn's own result has) -
  # bind_rows() fills the gaps with NA instead of failing on the mismatch.
  run_log <- dplyr::bind_rows(run_log)

  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(run_log, output_csv, row.names = FALSE)

  cat("\n========================================\n")
  cat("Grid search completed!\n")
  cat(sprintf("Total runs: %d\n", total_runs))
  cat("========================================\n")
  cat(sprintf("\nGrid search log saved to: %s\n", output_csv))
  cat("\nParameter combinations tested:\n")
  print(run_log)

  run_log
}

#' Tune hyperparameters separately for each (n_cells, n_genes) combination
#'
#' Shared driver behind cyclum_gs.R/scPrisma_gs.R/oscope_gs.R: for each combination
#' in `combos`, runs the full `param_grid` sweep via run_grid_search() against that
#' combination's own tuning-set directory, picks the row with the highest
#' `mean_tuning_auc` (run_fn is expected to include that column in what it returns -
#' see synthetic/dyngen_utils.R's score_against_ground_truth()), and combines the
#' per-combination winners into one summary.
#'
#' @param combos Data frame from list_combo_dirs() (columns: combo, path, n_cells, n_genes)
#' @param param_grid Data frame, one row per hyperparameter combination to try
#'   (the same grid is used for every (n_cells, n_genes) combination)
#' @param make_run_fn Function(combo_dir, fnames) -> a run_fn(params, run_number)
#'   closure suitable for run_grid_search(), letting each algorithm capture its own
#'   combo_dir/fnames for scoring against that combination's tuning replicates
#' @param algorithm Algorithm name, used in log messages and the output filename
#' @return The combined best-hyperparameters data frame (one row per combination,
#'   also written to <parent of combos' directories>/best_hyperparams_<algorithm>.csv)
tune_per_combo <- function(combos, param_grid, make_run_fn, algorithm) {

  library(dplyr)

  if (nrow(combos) == 0) {
    stop("No tuning-set combinations found - run generate_datasets.R first.")
  }

  best_hyperparams <- list()

  for (i in seq_len(nrow(combos))) {
    combo_dir  <- combos$path[i]
    combo_name <- combos$combo[i]

    cat("\n\n########################################\n")
    cat(sprintf("### Tuning %s for %s (n_cells=%d, n_genes=%d)\n",
                algorithm, combo_name, combos$n_cells[i], combos$n_genes[i]))
    cat("########################################\n")

    tuning_files <- list.files(combo_dir, pattern = "\\.h5$", full.names = FALSE)
    tuning_files <- tuning_files[!grepl("_sim\\.h5$", tuning_files)]
    fnames <- tools::file_path_sans_ext(tuning_files)

    if (length(fnames) == 0) {
      cat("No tuning .h5 files found in", combo_dir, "- skipping.\n")
      next
    }

    run_fn     <- make_run_fn(combo_dir, fnames)
    output_csv <- file.path(combo_dir, "grid_search_log.csv")
    combo_log  <- run_grid_search(param_grid, run_fn, output_csv)

    if (!"mean_tuning_auc" %in% names(combo_log) || all(is.na(combo_log$mean_tuning_auc))) {
      cat("No valid mean_tuning_auc values for", combo_name, "- skipping (every run may have errored).\n")
      next
    }

    best_row <- combo_log[which.max(combo_log$mean_tuning_auc), ]
    best_hyperparams[[combo_name]] <- cbind(
      data.frame(n_cells = combos$n_cells[i], n_genes = combos$n_genes[i]),
      dplyr::select(best_row, -run_counter)
    )
  }

  best_hyperparams_df <- dplyr::bind_rows(best_hyperparams)

  gridsearch_root <- dirname(combos$path[1])
  best_csv <- file.path(gridsearch_root, paste0("best_hyperparams_", algorithm, ".csv"))
  write.csv(best_hyperparams_df, best_csv, row.names = FALSE)

  cat("\n\nBest hyperparameters per combination saved to:", best_csv, "\n")
  print(best_hyperparams_df)

  best_hyperparams_df
}
