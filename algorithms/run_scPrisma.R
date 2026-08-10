library(reticulate)

#' Run scPrisma on a directory of raw-UMI-count .h5 files
#'
#' Thin R wrapper around algorithms/run_scPrisma.py, invoked via reticulate's
#' currently-selected Python (call use_condaenv("scPrisma_env", required = TRUE)
#' before this). Argument names/defaults mirror the Python CLI exactly.
#'
#' @param input_dir Directory containing input .h5 count files
#' @param output_dir Output directory for results (default: <input_dir>/scPrisma, set by the
#'   Python script itself if NULL)
#' @param regularisation_strength Numeric, regularisation strength for cyclic gene
#'   filtering (default 0.1)
#' @param iternum Integer, number of iterations (default 100)
#' @param n_top_genes Optional integer, number of highly variable genes to select
#'   (default: Python script's own default of 5000 if NULL)
#' @param skip_if_exists If TRUE, skip (rather than recompute and overwrite) any file
#'   whose output already exists - for resuming an interrupted evaluation run. Leave
#'   FALSE for grid search, which must always recompute since the same file's output
#'   path is reused across different hyperparameter combinations.
#' @param proj_root Project root directory (default here::here())
#' @return The exit code from the underlying Python process (invisibly)
run_scPrisma <- function(input_dir, output_dir = NULL,
                          regularisation_strength = 0.1, iternum = 100L, n_top_genes = NULL,
                          skip_if_exists = FALSE, proj_root = here::here()) {

  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  args <- c(
    "-u", file.path(proj_root, "algorithms/run_scPrisma.py"),
    input_dir,
    "--regularisation_strength", as.character(regularisation_strength),
    "--iternum", as.character(iternum)
  )
  if (!is.null(output_dir))  args <- c(args, "--output_dir", output_dir)
  if (!is.null(n_top_genes)) args <- c(args, "--n_top_genes", as.character(n_top_genes))
  if (skip_if_exists)        args <- c(args, "--skip_if_exists")

  exit_code <- system2(
    command = reticulate::py_config()$python,
    args    = args,
    stdout  = "",
    stderr  = ""
  )

  invisible(exit_code)
}
