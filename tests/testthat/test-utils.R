source(file.path(proj_root, "common/utils.R"))

test_that("qc_cell_mask keeps cells passing both the count and mito-fraction thresholds", {
  genes <- c("MT-ND1", "MT-CO1", "GENE1", "GENE2")
  # cell1: 1000 total, 20% mito (fails max_mito_frac); cell2: 500 total (fails min_counts);
  # cell3: 1000 total, 10% mito (passes, boundary-inclusive); cell4: 2000 total, 5% mito (passes)
  counts <- matrix(
    c(100, 100, 400, 400,
       50,  50, 200, 200,
       50,  50, 450, 450,
       50,  50, 950, 950),
    nrow = 4, dimnames = list(genes, paste0("cell", 1:4))
  )

  mask <- qc_cell_mask(counts, mito_prefix = "^MT-", min_counts = 1000, max_mito_frac = 0.10)

  expect_equal(unname(mask), c(FALSE, FALSE, TRUE, TRUE))
})

test_that("qc_cell_mask handles no matching mitochondrial genes (mito_frac = 0)", {
  genes <- c("GENE1", "GENE2")
  # both cells sum to 1200 (>= min_counts); no gene matches mito_prefix, so mito_frac = 0
  counts <- matrix(c(600, 600, 600, 600), nrow = 2, dimnames = list(genes, c("cell1", "cell2")))

  mask <- qc_cell_mask(counts, mito_prefix = "^MT-", min_counts = 1000, max_mito_frac = 0.10)

  expect_equal(unname(mask), c(TRUE, TRUE))
})

test_that("get_cyclum_scores computes magnitude from the last two weight columns", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  genes <- c("geneA", "geneB", "geneC")
  expr <- matrix(1:9, nrow = 3, dimnames = list(genes, paste0("cell", 1:3)))
  # weight: 3 genes x 4 cols - only the last 2 columns (real, imag) matter
  weight <- matrix(
    c(0, 0, 0,      0, 0, 0,      3, 0, 5,      4, 0, 12),
    nrow = 3, dimnames = list(genes, NULL)
  )

  expr_file   <- file.path(d, "expr.h5")
  weight_file <- file.path(d, "weight.h5")
  mat2hdf(expr, expr_file)
  mat2hdf(weight, weight_file)

  scores <- get_cyclum_scores(expr_file, weight_file)

  expect_equal(scores$symbol, genes)
  # abs(3+4i) = 5, abs(0+0i) = 0, abs(5+12i) = 13
  expect_equal(scores$score, c(5, 0, 13))
})

test_that("get_model_scores loads cyclum results and returns NULL if the weight file is missing", {
  d <- tempfile()
  dir.create(file.path(d, "cyclum"), recursive = TRUE)
  on.exit(unlink(d, recursive = TRUE))

  genes <- c("g1", "g2")
  expr   <- matrix(1:4, nrow = 2, dimnames = list(genes, c("c1", "c2")))
  weight <- matrix(c(0, 0, 1, 1, 3, 0, 4, 0), nrow = 2, dimnames = list(genes, NULL))

  mat2hdf(expr, file.path(d, "myfile.h5"))
  mat2hdf(weight, file.path(d, "cyclum", "myfile.h5"))

  result <- get_model_scores("cyclum", d, "myfile")
  expect_equal(result$symbol, genes)
  expect_equal(unname(result$score), c(5, 0))

  expect_null(get_model_scores("cyclum", d, "no_such_file"))
})

test_that("get_model_scores loads scPrisma/oscope results (CSV) and returns NULL if missing", {
  d <- tempfile()
  dir.create(file.path(d, "scPrisma"), recursive = TRUE)
  dir.create(file.path(d, "oscope"), recursive = TRUE)
  on.exit(unlink(d, recursive = TRUE))

  write.csv(data.frame(symbol = c("a", "b"), score = c(0.1, 0.9)),
            file.path(d, "scPrisma", "myfile.csv"), row.names = FALSE)
  write.csv(data.frame(symbol = c("a", "b"), score = c(1, 0)),
            file.path(d, "oscope", "myfile.csv"), row.names = FALSE)

  sc <- get_model_scores("scPrisma", d, "myfile")
  expect_equal(sc$score, c(0.1, 0.9))

  os <- get_model_scores("oscope", d, "myfile")
  expect_equal(os$score, c(1, 0))

  expect_null(get_model_scores("scPrisma", d, "no_such_file"))
})

test_that("load_model_scores loads cyclum results from explicit results_dir/gene_expr_file", {
  # myc/nfkb-style layout: results live under a `results/<algorithm>` tree that's
  # separate from where the expression .h5 lives.
  results_dir <- tempfile()
  expr_dir    <- tempfile()
  dir.create(results_dir, recursive = TRUE)
  dir.create(expr_dir, recursive = TRUE)
  on.exit(unlink(c(results_dir, expr_dir), recursive = TRUE))

  genes  <- c("g1", "g2")
  expr   <- matrix(1:4, nrow = 2, dimnames = list(genes, c("c1", "c2")))
  weight <- matrix(c(0, 0, 1, 1, 3, 0, 4, 0), nrow = 2, dimnames = list(genes, NULL))

  expr_file <- file.path(expr_dir, "myfile.h5")
  mat2hdf(expr, expr_file)
  mat2hdf(weight, file.path(results_dir, "myfile.h5"))

  result <- load_model_scores("cyclum", results_dir, expr_file, "myfile")
  expect_equal(result$symbol, genes)
  expect_equal(unname(result$score), c(5, 0))

  expect_null(load_model_scores("cyclum", results_dir, expr_file, "no_such_file"))
})

test_that("load_model_scores loads scPrisma/oscope results (CSV) and returns NULL if missing", {
  results_dir <- tempfile()
  dir.create(results_dir, recursive = TRUE)
  on.exit(unlink(results_dir, recursive = TRUE))

  write.csv(data.frame(symbol = c("a", "b"), score = c(0.2, 0.8)),
            file.path(results_dir, "myfile.csv"), row.names = FALSE)

  result <- load_model_scores("scPrisma", results_dir, gene_expr_file = NULL, "myfile")
  expect_equal(result$score, c(0.2, 0.8))

  expect_null(load_model_scores("scPrisma", results_dir, gene_expr_file = NULL, "no_such_file"))
})

test_that("load_model_scores errors on an unknown algorithm", {
  expect_error(
    load_model_scores("not_a_real_algorithm", tempdir(), gene_expr_file = NULL, "myfile"),
    "Unknown algorithm"
  )
})

test_that("get_model_scores errors on an unknown algorithm", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  expect_error(get_model_scores("not_a_real_algorithm", d, "myfile"), "Unknown algorithm")
})
