library(reticulate)

#' Run Cyclum on a directory of raw-UMI-count .h5 files
#'
#' Thin R wrapper around algorithms/run_cyclum.py, invoked via reticulate's
#' currently-selected Python (call use_condaenv("cyclum_env", required = TRUE)
#' before this). Argument names/defaults mirror the Python CLI exactly.
#'
#' @param input_dir Directory containing input .h5 count files
#' @param output_dir Output directory for results (default: <input_dir>/cyclum, set by the
#'   Python script itself if NULL)
#' @param encoder_width Integer vector, encoder hidden-layer widths (default c(30, 20))
#' @param epochs Integer, training epochs (default 500)
#' @param learning_rate Numeric, training learning rate (default 2e-4)
#' @param skip_if_exists If TRUE, skip (rather than recompute and overwrite) any file
#'   whose output already exists - for resuming an interrupted evaluation run. Leave
#'   FALSE for grid search, which must always recompute since the same file's output
#'   path is reused across different hyperparameter combinations.
#' @param proj_root Project root directory (default here::here())
#' @return The exit code from the underlying Python process (invisibly)
run_cyclum <- function(input_dir, output_dir = NULL,
                        encoder_width = c(30L, 20L), epochs = 500L, learning_rate = 2e-4,
                        skip_if_exists = FALSE, proj_root = here::here()) {

  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  args <- c(
    "-u", file.path(proj_root, "algorithms/run_cyclum.py"),
    input_dir,
    "--encoder_width", as.character(encoder_width),
    "--epochs", as.character(epochs),
    "--learning_rate", sprintf("%.0e", learning_rate)
  )
  if (!is.null(output_dir)) args <- c(args, "--output_dir", output_dir)
  if (skip_if_exists) args <- c(args, "--skip_if_exists")

  exit_code <- system2(
    command = reticulate::py_config()$python,
    args    = args,
    stdout  = "",
    stderr  = ""
  )

  invisible(exit_code)
}
