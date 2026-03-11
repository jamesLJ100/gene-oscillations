source(file.path(proj_root, "hdfrw.R"))

get_scores <- function(expr_file, weight_file) {
  expression <- hdf2mat(expr_file)
  weight <- hdf2mat(weight_file)
  
  n <- ncol(weight)
  mag <- abs(weight[, n - 1] + 1i * weight[, n])
  cyclum_df <- data.frame(symbol = rownames(expression), score = mag)
  cyclum_df
}

get_model_scores <- function(algorithm, file_name, expr_file) {
  if (algorithm == "Cyclum") {
    weight_file <- file.path(cyclum_dir, paste0(file_name, "_Cyclum.h5"))
    if (!file.exists(weight_file)) return(NULL)
    get_scores(expr_file, weight_file)
    
  } else if (algorithm == "scPrisma") {
    results_file <- file.path(scPrisma_dir, paste0(file_name, "_r0_scPrisma.csv"))
    if (!file.exists(results_file)) return(NULL)
    read.csv(results_file)
    
  } else if (algorithm == "Oscope") {
    results_file <- file.path(oscope_dir, paste0(file_name, "_oscope.csv"))
    if (!file.exists(results_file)) return(NULL)
    read.csv(results_file)
    
  } else {
    stop(paste("Unknown algorithm:", algorithm))
  }
}