test_that("every .R file in the repository parses without a syntax error", {
  # Catches regressions in the many orchestration scripts (preprocess_*.R,
  # run_*.R, analyse.R) that aren't unit-testable directly - they source heavy
  # Bioconductor packages and/or drive conda envs via reticulate - but should
  # at minimum always be syntactically valid.
  r_files <- list.files(proj_root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
  r_files <- r_files[!grepl("/\\.git/", r_files)]

  bad <- character(0)
  for (f in r_files) {
    ok <- tryCatch({ parse(f); TRUE }, error = function(e) FALSE)
    if (!ok) bad <- c(bad, f)
  }

  expect_length(bad, 0)
  if (length(bad) > 0) {
    cat("Files that failed to parse:\n")
    cat(paste0(" - ", bad, collapse = "\n"), "\n")
  }
})
