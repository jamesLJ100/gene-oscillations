library(reticulate)
source(file.path(proj_root, "algorithms/run_cyclum.R"))
source(file.path(proj_root, "algorithms/run_scPrisma.R"))

# Both run_cyclum()/run_scPrisma() call reticulate::py_config() (namespace-
# qualified) for the python executable and system2() (unqualified) to actually
# invoke it. py_config() needs testthat's local_mocked_bindings() since it's a
# namespaced call into another package's namespace; system2() is unqualified,
# so (as elsewhere in this suite) it has to be assigned directly into
# .GlobalEnv to be visible to code that was source()'d there - a plain
# `system2 <- function(...) ...` inside a test_that() block wouldn't be seen.
mock_system2_capturing_args <- function() {
  captured <- new.env()
  old <- if (exists("system2", envir = .GlobalEnv, inherits = FALSE)) get("system2", envir = .GlobalEnv) else NULL
  assign("system2", function(command, args, ...) {
    captured$command <- command
    captured$args <- args
    0L
  }, envir = .GlobalEnv)
  withr_like_restore <- function() {
    if (is.null(old)) rm("system2", envir = .GlobalEnv) else assign("system2", old, envir = .GlobalEnv)
  }
  list(captured = captured, restore = withr_like_restore)
}

test_that("run_cyclum builds the correct args with explicit hyperparameters and output_dir", {
  m <- mock_system2_capturing_args()
  on.exit(m$restore())
  local_mocked_bindings(py_config = function(...) list(python = "FAKE_PYTHON"), .package = "reticulate")

  out_dir <- tempfile()
  exit_code <- run_cyclum(
    input_dir = "IN", output_dir = out_dir,
    encoder_width = c(30, 20), epochs = 500, learning_rate = 2e-4,
    proj_root = "PROOT"
  )

  expect_equal(exit_code, 0L)
  expect_equal(m$captured$command, "FAKE_PYTHON")
  args <- m$captured$args
  expect_equal(args[1:3], c("-u", "PROOT/algorithms/run_cyclum.py", "IN"))
  expect_true("--encoder_width" %in% args)
  idx <- which(args == "--encoder_width")
  expect_equal(args[(idx + 1):(idx + 2)], c("30", "20"))
  expect_true(all(c("--epochs", "500", "--learning_rate", "2e-04", "--output_dir", out_dir) %in% args))
  expect_true(dir.exists(out_dir))
  unlink(out_dir, recursive = TRUE)
})

test_that("run_cyclum omits --output_dir when output_dir is NULL (defaults let Python decide)", {
  m <- mock_system2_capturing_args()
  on.exit(m$restore())
  local_mocked_bindings(py_config = function(...) list(python = "FAKE_PYTHON"), .package = "reticulate")

  run_cyclum(input_dir = "IN", proj_root = "PROOT")

  args <- m$captured$args
  expect_false("--output_dir" %in% args)
  # defaults mirror the Python CLI's own defaults
  expect_true(all(c("--epochs", "500", "--learning_rate", "2e-04") %in% args))
})

test_that("run_scPrisma builds the correct args and respects output_dir/n_top_genes", {
  m <- mock_system2_capturing_args()
  on.exit(m$restore())
  local_mocked_bindings(py_config = function(...) list(python = "FAKE_PYTHON"), .package = "reticulate")

  out_dir <- tempfile()
  exit_code <- run_scPrisma(
    input_dir = "IN", output_dir = out_dir,
    regularisation_strength = 0.2, iternum = 100, n_top_genes = 3000,
    proj_root = "PROOT"
  )

  expect_equal(exit_code, 0L)
  args <- m$captured$args
  expect_equal(args[1:3], c("-u", "PROOT/algorithms/run_scPrisma.py", "IN"))
  expect_true(all(c("--regularisation_strength", "0.2", "--iternum", "100",
                     "--output_dir", out_dir, "--n_top_genes", "3000") %in% args))
  unlink(out_dir, recursive = TRUE)
})

test_that("run_scPrisma omits --output_dir/--n_top_genes when not supplied", {
  m <- mock_system2_capturing_args()
  on.exit(m$restore())
  local_mocked_bindings(py_config = function(...) list(python = "FAKE_PYTHON"), .package = "reticulate")

  run_scPrisma(input_dir = "IN", regularisation_strength = 0.1, iternum = 100, proj_root = "PROOT")

  args <- m$captured$args
  expect_false("--output_dir" %in% args)
  expect_false("--n_top_genes" %in% args)
  expect_true(all(c("--regularisation_strength", "0.1", "--iternum", "100") %in% args))
})

test_that("run_cyclum/run_scPrisma include --skip_if_exists only when requested", {
  m <- mock_system2_capturing_args()
  on.exit(m$restore())
  local_mocked_bindings(py_config = function(...) list(python = "FAKE_PYTHON"), .package = "reticulate")

  run_cyclum(input_dir = "IN", proj_root = "PROOT")
  expect_false("--skip_if_exists" %in% m$captured$args)

  run_cyclum(input_dir = "IN", skip_if_exists = TRUE, proj_root = "PROOT")
  expect_true("--skip_if_exists" %in% m$captured$args)

  run_scPrisma(input_dir = "IN", proj_root = "PROOT")
  expect_false("--skip_if_exists" %in% m$captured$args)

  run_scPrisma(input_dir = "IN", skip_if_exists = TRUE, proj_root = "PROOT")
  expect_true("--skip_if_exists" %in% m$captured$args)
})
