source(file.path(proj_root, "common/hdfrw.R"))

test_that("mat2hdf/hdf2mat round-trips a matrix with both dimnames", {
  m <- matrix(1:12, nrow = 3, ncol = 4,
              dimnames = list(c("geneA", "geneB", "geneC"), paste0("cell", 1:4)))
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f))

  mat2hdf(m, f)
  m2 <- hdf2mat(f)

  expect_equal(dim(m2), dim(m))
  expect_equal(unname(m2), unname(m))
  expect_equal(rownames(m2), rownames(m))
  expect_equal(colnames(m2), colnames(m))
})

test_that("mat2hdf/hdf2mat round-trips a matrix with no dimnames", {
  m <- matrix(rnorm(20), nrow = 4, ncol = 5)
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f))

  mat2hdf(m, f)
  m2 <- hdf2mat(f)

  expect_equal(unname(m2), unname(m))
  expect_null(rownames(m2))
  expect_null(colnames(m2))
})

test_that("mat2hdf/hdf2mat round-trips a matrix with only rownames", {
  m <- matrix(1:6, nrow = 2, dimnames = list(c("x", "y"), NULL))
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f))

  mat2hdf(m, f)
  m2 <- hdf2mat(f)

  expect_equal(rownames(m2), c("x", "y"))
  expect_null(colnames(m2))
})

test_that("mat2hdf/hdf2mat preserves non-integer (decimal) values", {
  m <- matrix(c(0.1, 2.5, -3.333333, 100.001), nrow = 2,
              dimnames = list(c("a", "b"), c("s1", "s2")))
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f))

  mat2hdf(m, f)
  m2 <- hdf2mat(f)

  expect_equal(unname(m2), unname(m), tolerance = 1e-9)
})

test_that("hdf2mat errors clearly on a file with no 'matrix' dataset", {
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f))

  # Write an h5 file that has some other dataset, but not 'matrix'
  f5 <- hdf5r::H5File$new(f, mode = "w")
  f5[["not_a_matrix"]] <- 1:5
  hdf5r::h5close(f5)

  expect_error(hdf2mat(f), "no matrix found")
})
