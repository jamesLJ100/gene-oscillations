rm(list=ls())
proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "synthetic/run_dyngen.R"))

generate_datasets <- function(n_cells, n_genes, n_replicates, is_tuning = FALSE) {

  model_config <- make_config(n_cells, n_genes)

  # Determine subdirectory based on is_tuning flag
  subdir <- if (is_tuning) file.path("synthetic", "data", "dyngen_new", "gridsearch") else file.path("synthetic", "data", "dyngen_new")
  dir.create(here::here(subdir), showWarnings = FALSE, recursive = TRUE)
  error_log <- here::here(subdir, "generation_errors.log")

  for (i in seq_len(n_replicates)) {
    # Requested dimensions, used to identify this dataset if run_simulation() itself
    # fails (the actual simulated dimensions, used for the saved filename below, aren't
    # known until it succeeds).
    fname_base <- sprintf("c%dg%d_%d", n_cells, n_genes, i)

    tryCatch({
      sim <- run_simulation(model_config)

      fname_base <- sprintf("c%dg%d_%d", ncol(sim$counts), nrow(sim$counts), i)

      h5_file <- here::here(subdir, paste0(fname_base, ".h5"))
      sim_file <- here::here(subdir, paste0(fname_base, "_sim.rds"))

      mat2hdf(sim$counts, h5_file)
      saveRDS(sim, sim_file)

      cat("Saved dataset", i, "to:", h5_file, "\n")
      cat("Dimensions:", nrow(sim$counts), "genes ×", ncol(sim$counts), "cells\n")
      cat("Row names (genes):", length(rownames(sim$counts)), "saved\n")
      cat("Column names (cells):", length(colnames(sim$counts)), "saved\n\n")
    }, error = function(e) {
      msg <- conditionMessage(e)
      cat("ERROR generating", fname_base, "(n_cells =", n_cells, ", n_genes =", n_genes,
          ", replicate =", i, "):", msg, "\n\n")
      write(
        sprintf("%s\tdataset=%s\tn_cells=%d\tn_genes=%d\treplicate=%d\terror=%s",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"), fname_base, n_cells, n_genes, i,
                gsub("\n", " ", msg)),
        file = error_log, append = TRUE
      )
    })
  }
}

generate_data <- function() {
  for (i in c(50, 200, 500, 2000, 5000)) {
    generate_datasets(n_cells = 1000, n_genes = i, n_replicates = 5, is_tuning = FALSE)
  }
  for (i in c(50, 100, 250, 500, 2500, 5000)) {
    generate_datasets(n_cells = i, n_genes = 200, n_replicates = 5, is_tuning = FALSE)
  }
  #generate_datasets(n_cells = 5000, n_genes = 200, n_replicates = 5, is_tuning = FALSE)
}
generate_data()

# Generate separate tuning datasets
#generate_datasets(n_cells = 1000, n_genes = 200, n_replicates = 10, is_tuning = TRUE)