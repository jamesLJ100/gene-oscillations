source(file.path(proj_root, "common/hdfrw.R"))

# Cell QC mask: minimum UMI count + maximum mitochondrial fraction
qc_cell_mask <- function(counts, mito_prefix, min_counts = 1000, max_mito_frac = 0.10) {
  lib_sizes  <- colSums(counts)
  mito_genes <- grep(mito_prefix, rownames(counts), value = TRUE)
  mito_frac  <- colSums(counts[mito_genes, , drop = FALSE]) / lib_sizes

  cat("Mitochondrial genes found:", length(mito_genes), "\n")
  cat("Library size per cell:\n"); print(summary(lib_sizes))
  cat("Mitochondrial fraction per cell:\n"); print(summary(mito_frac))

  qc_flag <- lib_sizes >= min_counts & mito_frac <= max_mito_frac
  cat(sprintf(
    "QC filter (min_counts=%d, max_mito_frac=%.2f): keeping %d / %d cells\n",
    min_counts, max_mito_frac, sum(qc_flag), length(qc_flag)
  ))
  qc_flag
}

get_cyclum_scores <- function(gene_expr_file, model_weights_file) {
  expression <- hdf2mat(gene_expr_file)
  weight <- hdf2mat(model_weights_file)
  
  n <- ncol(weight)
  mag <- abs(weight[, n - 1] + 1i * weight[, n])
  cyclum_df <- data.frame(symbol = rownames(expression), score = mag)
  cyclum_df
}

# Results live under data_dir/algorithm/; the expression .h5 lives at
# data_dir/input_file_name.h5. Used by synthetic/*_gs.R and evaluate*.R, where both
# are nested under the same data_dir.
get_model_scores <- function(algorithm, data_dir, input_file_name) {
  load_model_scores(
    algorithm,
    results_dir     = file.path(data_dir, algorithm),
    gene_expr_file  = file.path(data_dir, paste0(input_file_name, ".h5")),
    input_file_name = input_file_name
  )
}

# Same loading logic as get_model_scores(), but with results_dir/gene_expr_file passed
# in directly rather than derived from a shared data_dir - for callers whose results and
# expression files don't live in that data_dir/algorithm layout (myc/nfkb/
# segmentation_clock's analyse.R, where results live in their own dataset/results/ tree).
load_model_scores <- function(algorithm, results_dir, gene_expr_file, input_file_name) {
  if (algorithm == "cyclum") {
    weight_file <- file.path(results_dir, paste0(input_file_name, ".h5"))
    if (!file.exists(weight_file)) return(NULL)
    get_cyclum_scores(gene_expr_file, weight_file)

  } else if (algorithm %in% c("scPrisma", "oscope")) {
    results_file <- file.path(results_dir, paste0(input_file_name, ".csv"))
    if (!file.exists(results_file)) return(NULL)
    read.csv(results_file)

  } else {
    stop(paste("Unknown algorithm:", algorithm))
  }
}