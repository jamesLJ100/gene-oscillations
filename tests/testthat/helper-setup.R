# Run once before any test file, via testthat::test_dir()'s automatic helper-file
# loading (this file matches "helper*.R"). Defines proj_root the same way every
# script in this repo does.

proj_root <- here::here()
setwd(proj_root)

# --- Why some test files stub out library() ---------------------------------
# common/gsea_functions.R and algorithms/run_oscope.R each library() several
# heavy Bioconductor/domain packages (clusterProfiler, msigdbr, dorothea,
# TFEA.ChIP, Oscope, EBSeq, ...) that aren't installed in every environment this
# suite might run in - they're needed for a real analysis session (see the
# README's environment setup), not for running tests. To still test the *logic*
# in those files (argument dispatch, data wrangling, control flow) without
# requiring the full package stack, the relevant test files temporarily
# reassign `library` to a no-op before sourcing the file under test:
#
#   library <- function(pkg, ...) invisible(NULL)
#   source(file.path(proj_root, "common/gsea_functions.R"))
#   rm(library)
#
# This only works because those files never call a function from the stubbed
# packages at *source* time - only inside function bodies, which the affected
# tests then call with their own mocked replacements for the handful of
# functions actually exercised (e.g. bitr(), GSEA()). Tests that need the real
# packages to be meaningful are wrapped in skip_if_not_installed() instead, so
# they skip cleanly here but still run in an environment with the full stack.
