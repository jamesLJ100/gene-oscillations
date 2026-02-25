library(here)
library(dplyr)
library(ROCR)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "hdfrw.R"))
source(file.path(proj_root, "cyclum_helper.R"))
source(file.path(proj_root, "dyngen_utils.R"))

direct_only <- FALSE

algorithms  <- c("scPrisma", "Cyclum")
dyngen_dir  <- file.path(proj_root, "data/dyngen")
cyclum_dir  <- file.path(proj_root, "data/dyngen/cyclum")
scPrisma_dir <- file.path(proj_root, "data/dyngen/scPrisma")
expr_files  <- list.files(dyngen_dir, pattern = "\\.h5$", full.names = TRUE)
results     <- list()

# runtimes must be defined somewhere — placeholder here
# runtimes <- read.csv(file.path(proj_root, "runtimes.csv"))

get_model_scores <- function(algorithm, file_name, expr_file) {
  if (algorithm == "Cyclum") {
    weight_file <- file.path(cyclum_dir, paste0(file_name, "_Cyclum.h5"))
    if (!file.exists(weight_file)) return(NULL)
    get_scores(expr_file, weight_file)
    
  } else if (algorithm == "scPrisma") {
    results_file <- file.path(scPrisma_dir, paste0(file_name, "_scPrisma.csv"))
    if (!file.exists(results_file)) return(NULL)
    read.csv(results_file)
    
  } else {
    stop(paste("Unknown algorithm:", algorithm))
  }
}

for (expr_file in expr_files) {
  fname <- tools::file_path_sans_ext(basename(expr_file))
  
  # Load simulation file once per expr_file
  sim_file <- file.path(proj_root, "data/dyngen", paste0(fname, "_sim.rds"))
  if (!file.exists(sim_file)) {
    cat("Skipping", fname, "- no corresponding sim RDS found\n")
    next
  }
  sim <- readRDS(sim_file)
  
  # Extract ground truth labels once
  feature_net <- sim[["model"]][["feature_network"]]
  
  tf_modules <- feature_net %>%
    filter(!is.na(from_module)) %>%
    select(gene = from, module = from_module) %>%
    distinct()
  
  target_edges <- feature_net %>%
    filter(grepl("^Target", to)) %>%
    select(from, to)
  
  gene_module_assignments <- propagate_module_assignments(tf_modules, target_edges)
  
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
  
  # Now process each algorithm for this file
  for (algorithm in algorithms) {
    cat("Processing:", fname, "with", algorithm, "\n")
    
    algorithm_results_df <- tryCatch(
      get_model_scores(algorithm, fname, expr_file),
      error = function(e) { cat("Error with", fname, algorithm, "-", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(algorithm_results_df)) {
      cat("Skipping", fname, algorithm, "- no results file found\n")
      next
    }
    
    merge_df <- merge(label_df, algorithm_results_df, by = "symbol", all.x = TRUE)
    # If algorithm filtered out a gene, treat it like a non-osc prediction
    merge_df$score[is.na(merge_df$score)] <- 0
    
    pred_obj  <- prediction(merge_df$score, merge_df$label)
    auc_val   <- performance(pred_obj, "auc")@y.values[[1]]
    
    runtime <- if (exists("runtimes")) {
      rt <- runtimes$runtime_seconds[runtimes$file == fname & runtimes$algorithm == algorithm]
      if (length(rt) == 1) rt else NA_real_
    } else NA_real_
    
    result_key <- paste(algorithm, fname, sep = "__")
    results[[result_key]] <- list(
      algorithm = algorithm,
      file      = fname,
      cells     = n_cells,
      genes     = n_genes,
      auc       = auc_val,
      runtime   = runtime
    )
    cat("AUC for", algorithm, fname, ":", auc_val, "\n")
  }
}

results_df <- do.call(rbind, lapply(results, as.data.frame))
rownames(results_df) <- NULL
print(results_df)

results_summary <- results_df %>%
  group_by(algorithm, cells, genes) %>%
  summarise(
    mean_auc     = mean(auc,     na.rm = TRUE),
    se_auc       = sd(auc,       na.rm = TRUE) / sqrt(n()),
    mean_runtime = mean(runtime, na.rm = TRUE),
    se_runtime   = sd(runtime,   na.rm = TRUE) / sqrt(n()),
    n            = n(),
    .groups      = "drop"
  )
print(results_summary)