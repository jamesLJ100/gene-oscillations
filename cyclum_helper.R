source(file.path(proj_root, "hdfrw.R"))

get_scores <- function(expr_file, weight_file) {
  expression <- hdf2mat(expr_file)
  weight <- hdf2mat(weight_file)
  
  n <- ncol(weight)
  mag <- abs(weight[, n - 1] + 1i * weight[, n])
  cyclum_df <- data.frame(symbol = rownames(expression), score = mag)
  cyclum_df
}