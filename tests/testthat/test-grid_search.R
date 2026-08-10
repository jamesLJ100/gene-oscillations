source(file.path(proj_root, "common/grid_search.R"))
source(file.path(proj_root, "synthetic/utils/dyngen_utils.R"))  # list_combo_dirs()

test_that("run_grid_search logs one row per parameter combination", {
  out_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(out_csv))

  param_grid <- expand.grid(alpha = c(1, 2, 3), beta = c("x", "y"), stringsAsFactors = FALSE)
  run_fn <- function(params, run_number) {
    list(exit_code = 0, product = params$alpha * run_number)
  }

  log <- run_grid_search(param_grid, run_fn, out_csv)

  expect_equal(nrow(log), 6)
  expect_equal(log$run_counter, 1:6)
  expect_equal(log$product[1], 1 * 1)
  expect_true(file.exists(out_csv))
  expect_equal(nrow(read.csv(out_csv)), 6)
})

test_that("run_grid_search catches per-run errors and still logs every row (mismatched columns)", {
  out_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(out_csv))

  param_grid <- expand.grid(alpha = c(1, 2, 3), beta = c("x", "y"), stringsAsFactors = FALSE)
  run_fn <- function(params, run_number) {
    if (params$alpha == 2 && params$beta == "y") stop("simulated failure")
    list(exit_code = 0, product = params$alpha * run_number)
  }

  log <- suppressWarnings(run_grid_search(param_grid, run_fn, out_csv))
  failed_row <- log[log$alpha == 2 & log$beta == "y", ]

  expect_equal(nrow(log), 6)
  expect_equal(failed_row$error, "simulated failure")
  expect_true(is.na(failed_row$exit_code))
  expect_true(is.na(failed_row$product))
  expect_true(all(is.na(log$error[log$alpha != 2 | log$beta != "y"])))
})

test_that("run_grid_search handles run_fn returning an empty list()", {
  out_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(out_csv))

  param_grid <- data.frame(x = 1:3)
  log <- run_grid_search(param_grid, function(params, run_number) list(), out_csv)

  expect_equal(nrow(log), 3)
  expect_equal(names(log), c("run_counter", "x"))
})

test_that("run_grid_search resumes from an existing output_csv, skipping completed rows", {
  out_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(out_csv))

  param_grid <- expand.grid(alpha = c(1, 2, 3), beta = c("x", "y"), stringsAsFactors = FALSE)

  # Simulate a run interrupted after 3 of 6 combinations completed.
  partial <- data.frame(
    run_counter = 1:3, alpha = c(1, 2, 3), beta = c("x", "x", "x"),
    mean_tuning_auc = c(0.5, 0.6, 0.7)
  )
  write.csv(partial, out_csv, row.names = FALSE)

  called_with <- list()
  run_fn <- function(params, run_number) {
    called_with[[length(called_with) + 1]] <<- params
    list(mean_tuning_auc = 0.9)
  }

  log <- run_grid_search(param_grid, run_fn, out_csv)

  # Only the 3 missing combinations should have actually invoked run_fn.
  expect_equal(length(called_with), 3)
  called_betas <- vapply(called_with, function(p) p$beta, character(1))
  expect_true(all(called_betas == "y"))

  # The 3 resumed rows keep their original scores; the 3 fresh ones get the new value.
  expect_equal(nrow(log), 6)
  expect_equal(log$mean_tuning_auc[log$beta == "x"], c(0.5, 0.6, 0.7))
  expect_true(all(log$mean_tuning_auc[log$beta == "y"] == 0.9))
})

test_that("run_grid_search retries a previously-errored combination rather than treating it as done", {
  out_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(out_csv))

  param_grid <- data.frame(alpha = c(1, 2), stringsAsFactors = FALSE)

  # Row for alpha=1 previously errored (no mean_tuning_auc); row for alpha=2 succeeded.
  partial <- data.frame(
    run_counter = 1:2, alpha = c(1, 2),
    mean_tuning_auc = c(NA_real_, 0.8), error = c("boom", NA_character_)
  )
  write.csv(partial, out_csv, row.names = FALSE)

  called_with <- list()
  run_fn <- function(params, run_number) {
    called_with[[length(called_with) + 1]] <<- params
    list(mean_tuning_auc = 0.95)
  }

  log <- run_grid_search(param_grid, run_fn, out_csv)

  expect_equal(length(called_with), 1)
  expect_equal(called_with[[1]]$alpha, 1)
  expect_equal(log$mean_tuning_auc[log$alpha == 1], 0.95)
  expect_equal(log$mean_tuning_auc[log$alpha == 2], 0.8)
})

test_that("run_grid_search only resumes rows still present in a changed param_grid", {
  out_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(out_csv))

  # Previous run swept alpha = 1, 2; this run widens it to include 3.
  partial <- data.frame(run_counter = 1:2, alpha = c(1, 2), mean_tuning_auc = c(0.5, 0.6))
  write.csv(partial, out_csv, row.names = FALSE)

  param_grid <- data.frame(alpha = c(1, 2, 3), stringsAsFactors = FALSE)

  called_with <- list()
  run_fn <- function(params, run_number) {
    called_with[[length(called_with) + 1]] <<- params
    list(mean_tuning_auc = 0.99)
  }

  log <- run_grid_search(param_grid, run_fn, out_csv)

  expect_equal(length(called_with), 1)
  expect_equal(called_with[[1]]$alpha, 3)
  expect_equal(nrow(log), 3)
  expect_equal(log$mean_tuning_auc[log$alpha == 3], 0.99)
})

test_that("tune_per_combo picks the best row per combination independently", {
  base_dir <- tempfile()
  dir.create(file.path(base_dir, "c50g200"), recursive = TRUE)
  dir.create(file.path(base_dir, "c1000g200"), recursive = TRUE)
  on.exit(unlink(base_dir, recursive = TRUE))
  for (i in 1:3) {
    file.create(file.path(base_dir, "c50g200", sprintf("c50g200_%d.h5", i)))
    file.create(file.path(base_dir, "c1000g200", sprintf("c1000g200_%d.h5", i)))
  }

  combos <- list_combo_dirs(base_dir)
  param_grid <- expand.grid(alpha = c(1, 2, 3), stringsAsFactors = FALSE)

  make_run_fn <- function(combo_dir, fnames) {
    force(combo_dir); force(fnames)
    function(params, run_number) {
      is_best <- (grepl("c50g200", combo_dir) && params$alpha == 2) ||
                 (grepl("c1000g200", combo_dir) && params$alpha == 3)
      list(mean_tuning_auc = if (is_best) 0.95 else 0.5, n_tuning_files = length(fnames))
    }
  }

  best_df <- suppressWarnings(tune_per_combo(combos, param_grid, make_run_fn, "fakealgo"))

  expect_equal(nrow(best_df), 2)
  expect_equal(best_df$alpha[best_df$n_cells == 50], 2)
  expect_equal(best_df$alpha[best_df$n_cells == 1000], 3)
  expect_true(file.exists(file.path(base_dir, "best_hyperparams_fakealgo.csv")))
  expect_true(file.exists(file.path(base_dir, "c50g200", "grid_search_log_fakealgo.csv")))
})

test_that("tune_per_combo checkpoints best_hyperparams during a combo's sweep, not only after it fully completes", {
  base_dir <- tempfile()
  dir.create(file.path(base_dir, "c50g200"), recursive = TRUE)
  on.exit(unlink(base_dir, recursive = TRUE))
  file.create(file.path(base_dir, "c50g200", "c50g200_1.h5"))

  combos <- list_combo_dirs(base_dir)
  param_grid <- data.frame(alpha = c(1, 2, 3), stringsAsFactors = FALSE)
  best_csv <- file.path(base_dir, "best_hyperparams_fakealgo.csv")

  seen_after_row1 <- NULL
  make_run_fn <- function(combo_dir, fnames) {
    function(params, run_number) {
      if (params$alpha == 2) {
        # By the time row 2 runs, row 1's result should already be on disk - this is
        # what would let a wall-clock kill here still leave a usable best_hyperparams
        # file, instead of only ever writing one once the whole combo sweep finishes.
        seen_after_row1 <<- file.exists(best_csv)
      }
      list(mean_tuning_auc = params$alpha / 10, n_tuning_files = length(fnames))
    }
  }

  tune_per_combo(combos, param_grid, make_run_fn, "fakealgo")

  expect_true(isTRUE(seen_after_row1))
})

test_that("tune_per_combo namespaces grid_search_log by algorithm, so two algorithms tuning the same combo don't overwrite each other", {
  base_dir <- tempfile()
  dir.create(file.path(base_dir, "c50g200"), recursive = TRUE)
  on.exit(unlink(base_dir, recursive = TRUE))
  file.create(file.path(base_dir, "c50g200", "c50g200_1.h5"))

  combos <- list_combo_dirs(base_dir)

  make_run_fn <- function(param_name) {
    function(combo_dir, fnames) {
      force(combo_dir); force(fnames)
      function(params, run_number) list(mean_tuning_auc = 0.5, n_tuning_files = length(fnames))
    }
  }

  tune_per_combo(combos, data.frame(alpha = 1), make_run_fn("alpha"), "algo_one")
  tune_per_combo(combos, data.frame(beta = 1), make_run_fn("beta"), "algo_two")

  log_one <- read.csv(file.path(base_dir, "c50g200", "grid_search_log_algo_one.csv"))
  log_two <- read.csv(file.path(base_dir, "c50g200", "grid_search_log_algo_two.csv"))
  expect_true("alpha" %in% names(log_one))
  expect_true("beta" %in% names(log_two))
  expect_false("alpha" %in% names(log_two))
})

test_that("tune_per_combo errors clearly when there are no combinations", {
  empty_combos <- data.frame(combo = character(0), path = character(0),
                              n_cells = integer(0), n_genes = integer(0))
  expect_error(
    tune_per_combo(empty_combos, data.frame(alpha = 1), function(...) function(...) list(), "x"),
    "No tuning-set combinations found"
  )
})

test_that("tune_per_combo skips (not errors) a combination directory with no .h5 files", {
  base_dir <- tempfile()
  dir.create(file.path(base_dir, "c99g99"), recursive = TRUE)
  on.exit(unlink(base_dir, recursive = TRUE))

  combos <- list_combo_dirs(base_dir)
  make_run_fn <- function(combo_dir, fnames) function(params, run_number) list(mean_tuning_auc = 1)

  best_df <- tune_per_combo(combos, data.frame(alpha = 1), make_run_fn, "fakealgo")
  expect_equal(nrow(best_df), 0)
})
