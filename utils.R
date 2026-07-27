source(file.path(proj_root, "hdfrw.R"))

get_cyclum_scores <- function(gene_expr_file, model_weights_file) {
  expression <- hdf2mat(gene_expr_file)
  weight <- hdf2mat(model_weights_file)
  
  n <- ncol(weight)
  mag <- abs(weight[, n - 1] + 1i * weight[, n])
  cyclum_df <- data.frame(symbol = rownames(expression), score = mag)
  cyclum_df
}

get_model_scores <- function(algorithm, data_dir, input_file_name, run_id=NULL) {
  if (is.null(run_id)) {
    input_file_name_with_run_id <- input_file_name
  } else {
    input_file_name_with_run_id <- paste0(input_file_name, "_r", run_id)
  }
  results_dir <- file.path(data_dir, algorithm)
  if (algorithm == "cyclum") {
    
    weight_file <- file.path(results_dir, paste0(input_file_name_with_run_id, ".h5"))
    if (!file.exists(weight_file)) return(NULL)
    gene_expr_file <- file.path(data_dir, paste0(input_file_name, ".h5"))
    get_cyclum_scores(gene_expr_file, weight_file)
    
  } else if (algorithm == "scPrisma") {
    results_file <- file.path(results_dir, paste0(input_file_name_with_run_id, ".csv"))
    if (!file.exists(results_file)) return(NULL)
    read.csv(results_file)
    
  } else if (algorithm == "oscope") {
    results_file <- file.path(results_dir, paste0(input_file_name_with_run_id, ".csv"))
    if (!file.exists(results_file)) return(NULL)
    read.csv(results_file)
    
  } else {
    stop(paste("Unknown algorithm:", algorithm))
  }
}