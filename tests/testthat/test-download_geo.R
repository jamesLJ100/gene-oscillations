source(file.path(proj_root, "common/download_geo.R"))

test_that("download_geo_supp_if_missing skips the download when all required files exist", {
  skip_if_not_installed("GEOquery")

  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  f1 <- file.path(d, "a.csv.gz")
  f2 <- file.path(d, "b.csv.gz")
  file.create(f1, f2)

  # When all required_files already exist, GEOquery::getGEOSuppFiles() (network I/O)
  # must never run - the only observable side effect is the skip message.
  expect_output(
    download_geo_supp_if_missing("GSEXXXXXX", d, c(f1, f2), label = "Test dataset"),
    "Test dataset already present at.*skipping download"
  )
})
