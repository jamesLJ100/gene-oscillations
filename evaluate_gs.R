library(here)
library(dplyr)
library(ROCR)
library(tidyr)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "hdfrw.R"))
source(file.path(proj_root, "utils.R"))
source(file.path(proj_root, "dyngen_utils.R"))

direct_only <- FALSE

#algorithms  <- c("scPrisma", "cyclum", "oscope")
algorithms  <- c("scPrisma")

dyngen_dir  <- file.path(proj_root, "data/dyngen/gridsearch")
cyclum_dir  <- file.path(dyngen_dir, "cyclum")
scPrisma_dir <- file.path(dyngen_dir, "scPrisma")
oscope_dir <- file.path(dyngen_dir, "oscope")
expr_files  <- list.files(dyngen_dir, pattern = "\\.h5$", full.names = TRUE)

# Load hyperparameter grid search logs for each algorithm
algorithm_dirs <- list(
  "cyclum" = cyclum_dir,
  "scPrisma" = scPrisma_dir,
  "oscope" = oscope_dir
)

hyperparam_configs <- list()
for (algo in names(algorithm_dirs)) {
  grid_search_file <- file.path(algorithm_dirs[[algo]], "grid_search_log.csv")
  if (file.exists(grid_search_file)) {
    hyperparam_configs[[algo]] <- read.csv(grid_search_file, stringsAsFactors = FALSE)
    cat("Loaded", nrow(hyperparam_configs[[algo]]), "hyperparameter configurations for", algo, "\n")
    print(hyperparam_configs[[algo]])
  } else {
    cat("Warning: grid_search_log.csv not found for", algo, "at:", grid_search_file, "\n")
    hyperparam_configs[[algo]] <- data.frame(run_counter = integer(0))
  }
}

results <- list()

# Helper function to extract run_id from filename
# Expects format: c1000g204_1_r1_cyclum.h5
extract_run_id <- function(filename) {
  # Match pattern: _r<number>_
  match <- regmatches(filename, regexpr("_r(\\d+)_", filename))
  if (length(match) == 0) return(NA_integer_)
  as.integer(sub("_r(\\d+)_", "\\1", match))
}

# Helper function to extract dataset ID (e.g., "1" from c1000g204_1_r1_cyclum.h5)
extract_dataset_id <- function(filename) {
  # Match pattern: after g<digits>_<dataset_id>_r
  match <- regmatches(filename, regexpr("g\\d+_(\\d+)_r", filename))
  if (length(match) == 0) return(NA_character_)
  sub("g\\d+_(\\d+)_r", "\\1", match)
}

for (expr_file in expr_files) {
  fname <- tools::file_path_sans_ext(basename(expr_file))
  
  # Load simulation file once per expr_file
  sim_file <- file.path(dyngen_dir, paste0(fname, "_sim.rds"))
  if (!file.exists(sim_file)) {
    cat("Skipping", fname, "- no corresponding sim RDS found\n")
    next
  }
  sim <- readRDS(sim_file)
  gene_module_assignments <- propagate_module_assignments(sim, context = fname)
  
  cycling_genes <- gene_module_assignments %>%
    filter(module %in% c("B", "C", "D")) %>%
    { if (direct_only) filter(., hops <= 1) else . } %>%
    pull(gene)
  
  gene_symbols <- rownames(sim$expression)
  
  eval_symbols <- if (direct_only) {
    max_hop_genes <- gene_module_assignments %>% filter(hops <= 1) %>% pull(gene)
    gene_symbols[gene_symbols %in% max_hop_genes]
  } else {
    gene_symbols
  }
  
  label_df <- data.frame(symbol = eval_symbols, label = 0L)
  label_df$label[eval_symbols %in% cycling_genes] <- 1L
  
  # Parse metadata once
  parsed <- regmatches(fname, regexpr("c(\\d+)g(\\d+)", fname))
  n_cells <- as.integer(sub("c(\\d+)g.*",  "\\1", parsed))
  n_genes <- as.integer(sub(".*g(\\d+).*", "\\1", parsed))
  
  # Extract dataset ID from filename (e.g., "1" from c1000g204_1)
  dataset_id <- sub("c\\d+g\\d+_(\\d+)$", "\\1", fname)
  if (dataset_id == fname) dataset_id <- "1"  # fallback if no dataset ID
  
  # Now process each algorithm AND each run_id for this file
  for (algorithm in algorithms) {
    # Get hyperparameter configs for this algorithm
    algo_hyperparams <- hyperparam_configs[[algorithm]]
    if (is.null(algo_hyperparams) || nrow(algo_hyperparams) == 0) {
      cat("Skipping", algorithm, "- no hyperparameter configurations found\n")
      next
    }
    
    for (run_id in algo_hyperparams$run_counter) {
      cat("Processing:", fname, "with", algorithm, "run_id:", run_id, "\n")

      algorithm_results_df <- tryCatch(
        get_model_scores(algorithm, dyngen_dir, fname, run_id),
        error = function(e) { 
          cat("Error with", fname, algorithm, "run_id", run_id, "-", conditionMessage(e), "\n")
          NULL 
        }
      )
      
      if (is.null(algorithm_results_df)) {
        cat("Skipping", fname, algorithm, "run_id", run_id, "- no results file found\n")
        next
      }
      
      merge_df <- merge(label_df, algorithm_results_df, by = "symbol", all.x = TRUE)
      # If algorithm filtered out a gene, treat it like a non-osc prediction
      merge_df$score[is.na(merge_df$score)] <- 0
      
      pred_obj  <- prediction(merge_df$score, merge_df$label)
      auc_val   <- performance(pred_obj, "auc")@y.values[[1]]
      
      result_key <- paste(algorithm, fname, run_id, sep = "__")
      results[[result_key]] <- list(
        algorithm   = algorithm,
        file        = fname,
        dataset_id  = dataset_id,
        cells       = n_cells,
        genes       = n_genes,
        run_id      = run_id,
        auc         = auc_val
      )
      cat("AUC for", algorithm, fname, "run_id", run_id, ":", auc_val, "\n")
    }
  }
}

# Convert results to dataframe
results_df <- do.call(rbind, lapply(results, as.data.frame))
rownames(results_df) <- NULL

print("=== All Results (without hyperparameters) ===")
print(results_df)

# Process each algorithm separately
for (algo in algorithms) {
  cat("\n========================================\n")
  cat("Processing algorithm:", algo, "\n")
  cat("========================================\n")
  
  # Filter results for this algorithm
  algo_results <- results_df %>% filter(algorithm == algo)
  
  if (nrow(algo_results) == 0) {
    cat("No results found for", algo, "\n")
    next
  }
  
  # Get hyperparameters for this algorithm
  algo_hyperparams <- hyperparam_configs[[algo]]
  
  if (nrow(algo_hyperparams) == 0) {
    cat("No hyperparameter configurations found for", algo, "\n")
    next
  }
  
  # Merge with algorithm-specific hyperparameters
  algo_results_with_params <- algo_results %>%
    left_join(algo_hyperparams, by = c("run_id" = "run_counter"))
  
  print("=== Results with Hyperparameters ===")
  print(algo_results_with_params)
  
  # Save detailed results to algorithm directory
  algo_dir <- algorithm_dirs[[algo]]
  detailed_output <- file.path(algo_dir, paste0(algo, "_hyperparameter_results.csv"))
  write.csv(algo_results_with_params, detailed_output, row.names = FALSE)
  cat("\nDetailed results saved to:", detailed_output, "\n")
  
  # Get hyperparameter column names (excluding run_counter)
  hyperparam_cols <- setdiff(names(algo_hyperparams), "run_counter")
  
  # Create grouping columns dynamically
  group_cols <- c("algorithm", "run_id", hyperparam_cols)
  
  # Analyze performance by hyperparameters
  cat("\n=== Performance by Hyperparameter Configuration ===\n")
  hyperparam_summary <- algo_results_with_params %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      mean_auc = mean(auc, na.rm = TRUE),
      sd_auc   = sd(auc, na.rm = TRUE),
      se_auc   = sd(auc, na.rm = TRUE) / sqrt(n()),
      min_auc  = min(auc, na.rm = TRUE),
      max_auc  = max(auc, na.rm = TRUE),
      n_datasets = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_auc))
  
  print(hyperparam_summary)
  
  # Save hyperparameter summary
  summary_output <- file.path(algo_dir, paste0(algo, "_best_hyperparameters.csv"))
  write.csv(hyperparam_summary, summary_output, row.names = FALSE)
  cat("\nHyperparameter summary saved to:", summary_output, "\n")
  
  # Find best hyperparameters
  cat("\n=== Best Hyperparameter Configuration ===\n")
  best_config <- hyperparam_summary %>%
    slice_max(mean_auc, n = 1)
  print(best_config)
}
  