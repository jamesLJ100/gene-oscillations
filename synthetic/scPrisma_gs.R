library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_scPrisma.R"))
source(file.path(proj_root, "common/grid_search.R"))

use_condaenv("scPrisma_env", required = TRUE)

input_dir <- file.path(proj_root, "synthetic/data/dyngen/gridsearch")
iternum   <- 100

# Define grid search parameters
param_grid <- data.frame(
  regularisation_strength = c(0, 0.1, 0.2, 0.3, 0.4)
)

run_fn <- function(params, run_number) {
  exit_code <- run_scPrisma(
    input_dir               = input_dir,
    regularisation_strength = params$regularisation_strength,
    iternum                 = iternum,
    run_number              = run_number
  )
  if (exit_code != 0) {
    cat(sprintf("WARNING: Run %d exited with code %d\n", run_number, exit_code))
  }
  list(exit_code = exit_code)
}

output_csv <- file.path(proj_root, "synthetic/data/dyngen/gridsearch/scPrisma/grid_search_log.csv")
run_grid_search(param_grid, run_fn, output_csv)
