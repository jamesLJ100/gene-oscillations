# algorithms/run_oscope.R's apply_oscope() (the actual sine-fitting/clustering
# statistics) needs the real Oscope/EBSeq packages, which aren't installed in
# every environment (see helper-setup.R) - that part is a real testing gap here,
# only exercisable where those packages are available. The file I/O and
# control-flow around it (read_data_file, process_data_file, run_oscope's
# directory handling) doesn't depend on Oscope/EBSeq at all and is tested here
# by mocking apply_oscope with a lightweight stand-in.
source(file.path(proj_root, "common/hdfrw.R"))

# algorithms/run_oscope.R itself library()s Oscope (and hdf5r, already attached
# via common/hdfrw.R above) - Oscope isn't installed in every environment (see
# helper-setup.R). Stub library() during the source() call, same pattern as
# test-gsea_functions.R.
real_library <- base::library
assign("library", function(pkg, ...) invisible(NULL), envir = .GlobalEnv)
source(file.path(proj_root, "algorithms/run_oscope.R"))
assign("library", real_library, envir = .GlobalEnv)
rm(real_library)

mock_apply_oscope <- function(oscillating_genes = character(0)) {
  function(gene_expr_matrix, oscope_config) {
    genes <- rownames(gene_expr_matrix)
    list(
      gene_classification = data.frame(
        gene = genes,
        is_oscillating = genes %in% oscillating_genes,
        stringsAsFactors = FALSE
      ),
      SineRes = NULL
    )
  }
}

# apply_oscope() was source()'d into .GlobalEnv, so (as with the other mocked
# functions in this suite) it has to be reassigned there directly - a plain
# `apply_oscope <- function(...) ...` inside a test_that() block wouldn't be
# visible to process_data_file()/run_oscope(), whose lexical scope is .GlobalEnv.
with_mocked_apply_oscope <- function(mock_fn, code) {
  old <- get("apply_oscope", envir = .GlobalEnv)
  assign("apply_oscope", mock_fn, envir = .GlobalEnv)
  on.exit(assign("apply_oscope", old, envir = .GlobalEnv))
  force(code)
}

test_that("read_data_file loads an .h5 count matrix", {
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f))
  m <- matrix(1:6, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  mat2hdf(m, f)

  loaded <- read_data_file(f)
  expect_equal(unname(loaded), unname(m))
})

test_that("process_data_file writes a symbol/score CSV sorted by oscillating status", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  genes <- c("A", "B", "C", "D")
  m <- matrix(1:16, nrow = 4, dimnames = list(genes, paste0("s", 1:4)))
  f <- file.path(d, "myfile.h5")
  mat2hdf(m, f)

  with_mocked_apply_oscope(mock_apply_oscope(oscillating_genes = c("B", "D")), {
    result <- process_data_file(f, oscope_config = list())
  })

  expect_equal(result$file, "myfile")
  expect_true(is.na(result$error))
  expect_equal(result$n_genes, 4)
  expect_equal(result$n_oscillating, 2)

  csv_path <- file.path(d, "oscope", "myfile.csv")
  expect_true(file.exists(csv_path))
  gene_df <- read.csv(csv_path)
  expect_equal(names(gene_df), c("symbol", "score"))
  expect_equal(gene_df$score, sort(gene_df$score, decreasing = TRUE))  # oscillating genes first
  expect_setequal(gene_df$symbol[gene_df$score == 1], c("B", "D"))
})

test_that("process_data_file overwrites a previous run's output for the same file", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  genes <- c("A", "B")
  m <- matrix(1:4, nrow = 2, dimnames = list(genes, c("s1", "s2")))
  f <- file.path(d, "myfile.h5")
  mat2hdf(m, f)

  with_mocked_apply_oscope(mock_apply_oscope(oscillating_genes = "A"), {
    process_data_file(f, oscope_config = list())
  })
  with_mocked_apply_oscope(mock_apply_oscope(oscillating_genes = "B"), {
    process_data_file(f, oscope_config = list())
  })

  csv_path <- file.path(d, "oscope", "myfile.csv")
  gene_df <- read.csv(csv_path)
  # Only the second run's result should remain - no myfile_r<n>.csv files left behind
  expect_setequal(gene_df$symbol[gene_df$score == 1], "B")
  expect_equal(list.files(file.path(d, "oscope"), pattern = "^myfile"), "myfile.csv")
})

test_that("process_data_file with skip_if_exists = TRUE skips a file whose output already exists", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  genes <- c("A", "B")
  m <- matrix(1:4, nrow = 2, dimnames = list(genes, c("s1", "s2")))
  f <- file.path(d, "myfile.h5")
  mat2hdf(m, f)

  with_mocked_apply_oscope(mock_apply_oscope(oscillating_genes = "A"), {
    process_data_file(f, oscope_config = list())
  })

  called <- FALSE
  with_mocked_apply_oscope(function(...) { called <<- TRUE; mock_apply_oscope("B")(...) }, {
    result <- process_data_file(f, oscope_config = list(), skip_if_exists = TRUE)
  })

  expect_false(called)
  expect_true(is.na(result$runtime_seconds))
  # Original output (from the first, non-skipped run) is untouched.
  gene_df <- read.csv(file.path(d, "oscope", "myfile.csv"))
  expect_setequal(gene_df$symbol[gene_df$score == 1], "A")
})

test_that("process_data_file with skip_if_exists = TRUE still runs when no output exists yet", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  genes <- c("A", "B")
  m <- matrix(1:4, nrow = 2, dimnames = list(genes, c("s1", "s2")))
  f <- file.path(d, "myfile.h5")
  mat2hdf(m, f)

  called <- FALSE
  with_mocked_apply_oscope(function(...) { called <<- TRUE; mock_apply_oscope("A")(...) }, {
    result <- process_data_file(f, oscope_config = list(), skip_if_exists = TRUE)
  })

  expect_true(called)
  expect_true(is.na(result$error))
  gene_df <- read.csv(file.path(d, "oscope", "myfile.csv"))
  expect_setequal(gene_df$symbol[gene_df$score == 1], "A")
})

test_that("process_data_file catches errors from apply_oscope and returns an error row instead of throwing", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  m <- matrix(1:4, nrow = 2, dimnames = list(c("A", "B"), c("s1", "s2")))
  f <- file.path(d, "myfile.h5")
  mat2hdf(m, f)

  with_mocked_apply_oscope(function(...) stop("simulated Oscope failure"), {
    result <- process_data_file(f, oscope_config = list())
  })

  expect_equal(result$error, "simulated Oscope failure")
  expect_true(is.na(result$runtime_seconds))
  expect_false(file.exists(file.path(d, "oscope", "myfile.csv")))
})

test_that("run_oscope processes every .h5 file except _sim.h5 and writes a combined runtimes.csv", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  for (fn in c("c50g20_1", "c50g20_2")) {
    m <- matrix(1:4, nrow = 2, dimnames = list(c("A", "B"), c("s1", "s2")))
    mat2hdf(m, file.path(d, paste0(fn, ".h5")))
    saveRDS(list(dummy = TRUE), file.path(d, paste0(fn, "_sim.rds")))
  }
  # a _sim.h5 file that must NOT be picked up as a data file
  mat2hdf(matrix(1:4, nrow = 2), file.path(d, "c50g20_1_sim.h5"))

  with_mocked_apply_oscope(mock_apply_oscope(oscillating_genes = "A"), {
    run_oscope(d, oscope_config = list())
  })

  runtimes <- read.csv(file.path(d, "oscope", "runtimes.csv"))
  expect_equal(nrow(runtimes), 2)  # not 3 - the _sim.h5 file was excluded
  expect_setequal(runtimes$file, c("c50g20_1", "c50g20_2"))
})

test_that("run_oscope with skip_if_exists = TRUE only recomputes files without existing output", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  for (fn in c("c50g20_1", "c50g20_2")) {
    m <- matrix(1:4, nrow = 2, dimnames = list(c("A", "B"), c("s1", "s2")))
    mat2hdf(m, file.path(d, paste0(fn, ".h5")))
  }

  # Pre-populate c50g20_1's output, as if from an earlier, interrupted run.
  dir.create(file.path(d, "oscope"), recursive = TRUE)
  write.csv(data.frame(symbol = c("A", "B"), score = c(1, 0)),
            file.path(d, "oscope", "c50g20_1.csv"), row.names = FALSE)

  processed <- character(0)
  with_mocked_apply_oscope(function(gene_expr_matrix, ...) {
    processed <<- c(processed, "called")
    mock_apply_oscope("B")(gene_expr_matrix)
  }, {
    run_oscope(d, oscope_config = list(), skip_if_exists = TRUE)
  })

  expect_equal(length(processed), 1)  # only c50g20_2 was actually recomputed

  # c50g20_1's pre-existing output is untouched (still A, not B).
  gene_df_1 <- read.csv(file.path(d, "oscope", "c50g20_1.csv"))
  expect_setequal(gene_df_1$symbol[gene_df_1$score == 1], "A")

  gene_df_2 <- read.csv(file.path(d, "oscope", "c50g20_2.csv"))
  expect_setequal(gene_df_2$symbol[gene_df_2$score == 1], "B")
})

test_that("run_oscope errors (doesn't quit the R session) when the directory is missing or empty", {
  expect_error(run_oscope(tempfile(), oscope_config = list()), "does not exist")

  empty_dir <- tempfile()
  dir.create(empty_dir)
  on.exit(unlink(empty_dir, recursive = TRUE))
  expect_error(run_oscope(empty_dir, oscope_config = list()), "No H5 files found")
})
