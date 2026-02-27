source(file.path(proj_root, "run_dyngen.R"))

generate_datasets <- function(n_cells, n_genes, n_replicates) {
  
  model_config <- make_config(n_cells, n_genes)
  
  for (i in seq_len(n_replicates)) {
    sim <- run_simulation(model_config)
    
    fname_base <- sprintf("c%dg%d_%d", ncol(sim$expression), nrow(sim$expression), i)
    
    h5_file <- here::here("data", "dyngen", paste0(fname_base, ".h5"))
    sim_file <- here::here("data", "dyngen", paste0(fname_base, "_sim.rds"))
    dir.create(dirname(h5_file), showWarnings = FALSE, recursive = TRUE)
    
    mat2hdf(sim$expression, h5_file)
    saveRDS(sim, sim_file)
    
    cat("Saved dataset", i, "to:", h5_file, "\n")
    cat("Dimensions:", nrow(sim$expression), "genes ×", ncol(sim$expression), "cells\n")
    cat("Row names (genes):", length(rownames(sim$expression)), "saved\n")
    cat("Column names (cells):", length(colnames(sim$expression)), "saved\n\n")
  }
}

# generate_all_data <- function() {
#   for (i in c(50, 100, 250, 500, 1000, 2500, 5000)) {
#     for (j in c(50, 200, 500, 2000, 5000)) {
#       generate_datasets(n_cells = i, n_genes = j, n_replicates = 10)
#     }
#   }
# }

generate_datasets(n_cells = 1000, n_genes = 200, n_replicates = 1)