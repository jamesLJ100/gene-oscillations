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
