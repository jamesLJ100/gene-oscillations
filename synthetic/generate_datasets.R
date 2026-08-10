rm(list=ls())
proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "synthetic/run_dyngen.R"))

generate_datasets <- function(n_cells, n_genes, n_replicates, is_tuning = FALSE) {

  model_config <- make_config(n_cells, n_genes)

  # Each (n_cells, n_genes) combination gets its own subdirectory, for both the
  # tuning and evaluation sets. This isn't just organisational: run_cyclum.py/
  # run_scPrisma.py each glob every .h5 file in whatever directory they're pointed
  # at and apply one fixed hyperparameter setting to all of them in a single call -
  # so combinations that need different (tuned) hyperparameters must be physically
  # isolated in their own directory, not just distinguishable by filename.
  combo_dir <- sprintf("c%dg%d", n_cells, n_genes)
  subdir <- if (is_tuning) {
    file.path("synthetic", "data", "dyngen", "gridsearch", combo_dir)
  } else {
    file.path("synthetic", "data", "dyngen", combo_dir)
  }
  dir.create(here::here(subdir), showWarnings = FALSE, recursive = TRUE)
  error_log <- here::here(subdir, "generation_errors.log")

  for (i in seq_len(n_replicates)) {
    # Requested dimensions, used to identify this dataset if run_simulation() itself
    # fails (the actual simulated dimensions, used for the saved filename below, aren't
    # known until it succeeds).
    fname_base <- sprintf("c%dg%d_%d", n_cells, n_genes, i)

    # Skip already-generated replicates (resume support). Gene count in the filename
    # is the actual simulated count, not the requested n_genes, hence the wildcard.
    h5_pattern  <- sprintf("^c%dg[0-9]+_%d\\.h5$", n_cells, i)
    existing_h5 <- list.files(here::here(subdir), pattern = h5_pattern)

    already_done <- FALSE
    for (h5_name in existing_h5) {
      h5_path  <- here::here(subdir, h5_name)
      rds_path <- here::here(subdir, paste0(sub("\\.h5$", "", h5_name), "_sim.rds"))
      if (file.exists(rds_path) && file.size(h5_path) > 0 && file.size(rds_path) > 0) {
        already_done <- TRUE
        break
      }
    }

    if (already_done) {
      cat("Skipping", fname_base, "(n_cells =", n_cells, ", n_genes =", n_genes,
          ", replicate =", i, "): output already exists\n\n")
      next
    }

    tryCatch({
      sim <- run_simulation(model_config)

      fname_base <- sprintf("c%dg%d_%d", ncol(sim$counts), nrow(sim$counts), i)

      h5_file <- here::here(subdir, paste0(fname_base, ".h5"))
      sim_file <- here::here(subdir, paste0(fname_base, "_sim.rds"))

      # Write to temp names, then rename into place - avoids a truncated file at the
      # final path if the process is killed mid-write.
      h5_tmp  <- paste0(h5_file, ".tmp")
      sim_tmp <- paste0(sim_file, ".tmp")

      mat2hdf(sim$counts, h5_tmp)
      saveRDS(sim, sim_tmp)
      file.rename(h5_tmp, h5_file)
      file.rename(sim_tmp, sim_file)

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

# Two 1D sweeps sharing a common anchor (c1000g200): one varies genes with cells
# held fixed, the other varies cells with genes held fixed.
combos <- rbind(
  data.frame(n_cells = 1000, n_genes = c(50, 200, 500, 2000, 5000)),
  data.frame(n_cells = c(50, 100, 250, 500, 2500, 5000), n_genes = 200)
)

# For each (n_cells, n_genes) combination, generate two independent 5-replicate
# sets: a tuning set (for grid search - see gridsearch/cyclum_gs.R/scPrisma_gs.R/oscope_gs.R,
# which find the best hyperparameters per combination using only these) and a
# held-out evaluation set (scored in evaluate.R using those best hyperparameters,
# so the reported performance isn't inflated by evaluating on the same data the
# hyperparameters were chosen on).
generate_data <- function() {
  for (i in seq_len(nrow(combos))) {
    generate_datasets(n_cells = combos$n_cells[i], n_genes = combos$n_genes[i],
                       n_replicates = 5, is_tuning = FALSE)
    generate_datasets(n_cells = combos$n_cells[i], n_genes = combos$n_genes[i],
                       n_replicates = 5, is_tuning = TRUE)
  }

  # Errors are written per-combo (see generate_datasets() above), so scan for all of
  # them here rather than requiring a manual check of every combo's subdirectory.
  error_logs <- list.files(here::here("synthetic", "data", "dyngen"),
                            pattern = "^generation_errors\\.log$",
                            recursive = TRUE, full.names = TRUE)
  cat("\n========================================\n")
  if (length(error_logs) == 0) {
    cat("No generation errors.\n")
  } else {
    cat("Generation errors found in", length(error_logs), "combo(s):\n")
    for (log in error_logs) {
      cat(sprintf("  %s (%d error(s))\n", log, length(readLines(log))))
    }
  }
  cat("========================================\n")
}
generate_data()