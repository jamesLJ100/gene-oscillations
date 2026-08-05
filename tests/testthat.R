library(testthat)
library(here)

proj_root <- here::here()
setwd(proj_root)

test_dir(file.path(proj_root, "tests", "testthat"), reporter = "summary")
