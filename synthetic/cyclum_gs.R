library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_cyclum.R"))
source(file.path(proj_root, "common/grid_search.R"))

use_condaenv("cyclum_env", required = TRUE)

#TODO nonlinear_reg & check for others.

input_dir <- file.path(proj_root, "synthetic/data/dyngen/tuning")

# Define grid search parameters. encoder_width is vector-valued, so it's stored as a
# comma-separated string in param_grid (kept a plain flat data frame, parsed back to
# an integer vector inside run_fn) rather than a list-column.
encoder_widths_list <- list(
  c(30, 20),
  c(40, 30),
  c(50, 50),
  c(50, 40, 30),
  c(60, 40, 20)
)
encoder_width_strs <- sapply(encoder_widths_list, paste, collapse = ", ")

param_grid <- expand.grid(
  encoder_width = encoder_width_strs,
  epochs        = c(300, 500, 700),
  learning_rate = c(1e-4, 2e-4, 5e-4, 1e-3),
  stringsAsFactors = FALSE
)

run_fn <- function(params, run_number) {
  encoder_width <- as.integer(trimws(strsplit(params$encoder_width, ",")[[1]]))

  exit_code <- run_cyclum(
    input_dir     = input_dir,
    encoder_width = encoder_width,
    epochs        = params$epochs,
    learning_rate = params$learning_rate,
    run_number    = run_number
  )

  list(exit_code = exit_code)
}

output_csv <- file.path(proj_root, "synthetic/data/dyngen/tuning/cyclum/grid_search_log.csv")
run_grid_search(param_grid, run_fn, output_csv)
