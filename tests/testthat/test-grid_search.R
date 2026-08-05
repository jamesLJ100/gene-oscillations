source(file.path(proj_root, "common/grid_search.R"))
source(file.path(proj_root, "synthetic/dyngen_utils.R"))  # list_combo_dirs()

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
  expect_true(file.exists(file.path(base_dir, "c50g200", "grid_search_log.csv")))
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
