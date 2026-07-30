rm(list=ls())
library(here)
library(dplyr)
library(ROCR)
library(PRROC)
library(ggplot2)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/hdfrw.R"))
source(file.path(proj_root, "synthetic/dyngen_utils.R"))
source(file.path(proj_root, "common/utils.R"))


algorithms  <- c("scPrisma", "cyclum", "oscope")
dyngen_dir  <- file.path(proj_root, "synthetic/data/dyngen_new")
cyclum_dir  <- file.path(proj_root, "synthetic/data/dyngen_new/cyclum")
scPrisma_dir <- file.path(proj_root, "synthetic/data/dyngen_new/scPrisma")
oscope_dir <- file.path(proj_root, "synthetic/data/dyngen_new/oscope")
expr_files  <- list.files(dyngen_dir, pattern = "\\.h5$", full.names = TRUE)
results     <- list()

# Each algorithm writes a runtimes.csv (columns: file, runtime_seconds, ...) into
# its own output dir. Combine them into a single `runtimes` table tagged with the
# algorithm, averaging over any repeated runs (run_number) of the same file.
runtime_dirs <- c(scPrisma = scPrisma_dir, cyclum = cyclum_dir, oscope = oscope_dir)
runtimes <- do.call(rbind, lapply(names(runtime_dirs), function(algorithm) {
  rt_path <- file.path(runtime_dirs[[algorithm]], "runtimes.csv")
  if (!file.exists(rt_path)) {
    cat("No runtimes.csv for", algorithm, "at", rt_path, "\n")
    return(NULL)
  }
  rt <- read.csv(rt_path, stringsAsFactors = FALSE)
  rt <- rt[!is.na(rt$runtime_seconds), c("file", "runtime_seconds")]
  if (nrow(rt) == 0) return(NULL)
  # Average across repeated runs so each (file, algorithm) has one runtime
  rt <- aggregate(runtime_seconds ~ file, data = rt, FUN = mean)
  rt$algorithm <- algorithm
  rt
}))

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
  
  # Only evaluate on genes at most one hop from a TF (hops <= 1); this also
  # drops cycle genes with NA hops, since NA <= 1 is NA (filtered out).
  cycling_genes <- gene_module_assignments %>%
    filter(module %in% c("B", "C", "D"), hops <= 1) %>%
    pull(gene)

  gene_symbols <- rownames(sim$counts)

  eval_genes   <- gene_module_assignments %>% filter(hops <= 1) %>% pull(gene)
  eval_symbols <- gene_symbols[gene_symbols %in% eval_genes]

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
      get_model_scores(algorithm, dyngen_dir, fname),
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

    # AUPRC via PRROC (a proper precision-recall curve integration, not a naive
    # trapezoidal estimate over ROCR's cutoffs, which is biased under class imbalance).
    auprc_val <- tryCatch({
      pr <- pr.curve(
        scores.class0 = merge_df$score[merge_df$label == 1],
        scores.class1 = merge_df$score[merge_df$label == 0],
        curve = FALSE
      )
      pr$auc.integral
    }, error = function(e) NA_real_)

    # Threshold-dependent metrics (accuracy, F1, precision, recall) need a cutoff.
    # Algorithms produce scores on very different scales, so a fixed cutoff (e.g. 0.5)
    # wouldn't be comparable across them - instead use each (algorithm, file)'s own
    # F1-maximising cutoff, a standard choice for reporting a single operating point.
    f1_perf   <- performance(pred_obj, "f")
    acc_perf  <- performance(pred_obj, "acc")
    prec_perf <- performance(pred_obj, "prec")
    rec_perf  <- performance(pred_obj, "rec")

    f1_vals  <- f1_perf@y.values[[1]]
    best_idx <- which.max(f1_vals)

    f1_val   <- if (length(best_idx) > 0) f1_vals[best_idx]                    else NA_real_
    acc_val  <- if (length(best_idx) > 0) acc_perf@y.values[[1]][best_idx]     else NA_real_
    prec_val <- if (length(best_idx) > 0) prec_perf@y.values[[1]][best_idx]    else NA_real_
    rec_val  <- if (length(best_idx) > 0) rec_perf@y.values[[1]][best_idx]     else NA_real_

    runtime <- if (exists("runtimes")) {
      rt <- runtimes$runtime_seconds[runtimes$file == fname & runtimes$algorithm == algorithm]
      if (length(rt) == 1) rt else NA_real_
    } else NA_real_

    result_key <- paste(algorithm, fname, sep = "__")
    results[[result_key]] <- list(
      algorithm  = algorithm,
      file       = fname,
      cells      = n_cells,
      genes      = n_genes,
      eval_genes = nrow(merge_df),
      auc        = auc_val,
      auprc      = auprc_val,
      accuracy   = acc_val,
      f1         = f1_val,
      precision  = prec_val,
      recall     = rec_val,
      runtime    = runtime
    )
    cat("AUC for", algorithm, fname, ":", auc_val,
        "| AUPRC:", auprc_val, "| F1:", f1_val,
        "| Acc:", acc_val, "| Prec:", prec_val, "| Rec:", rec_val, "\n")
  }
}

results_df <- do.call(rbind, lapply(results, as.data.frame))
rownames(results_df) <- NULL
print(results_df)

results_summary <- results_df %>%
  group_by(algorithm, cells, genes) %>%
  summarise(
    mean_auc        = mean(auc,      na.rm = TRUE),
    sd_auc          = sd(auc,        na.rm = TRUE),
    se_auc          = sd(auc,        na.rm = TRUE) / sqrt(n()),
    mean_auprc      = mean(auprc,    na.rm = TRUE),
    sd_auprc        = sd(auprc,      na.rm = TRUE),
    se_auprc        = sd(auprc,      na.rm = TRUE) / sqrt(n()),
    mean_accuracy   = mean(accuracy, na.rm = TRUE),
    sd_accuracy     = sd(accuracy,   na.rm = TRUE),
    mean_f1         = mean(f1,       na.rm = TRUE),
    sd_f1           = sd(f1,         na.rm = TRUE),
    mean_precision  = mean(precision, na.rm = TRUE),
    sd_precision    = sd(precision,   na.rm = TRUE),
    mean_recall     = mean(recall,   na.rm = TRUE),
    sd_recall       = sd(recall,     na.rm = TRUE),
    mean_eval_genes = mean(eval_genes, na.rm = TRUE),
    mean_runtime    = mean(runtime, na.rm = TRUE),
    se_runtime      = sd(runtime,   na.rm = TRUE) / sqrt(n()),
    n               = n(),
    .groups         = "drop"
  )
print(results_summary)

# ---- AUROC sweeps: vs number of cells and vs number of genes ----------------
# The datasets form two sweeps sharing a common anchor (c1000g204): one varies
# cells with genes held fixed, the other varies genes with cells held fixed.
# Detect the held-fixed value automatically as the one covering the most points.
grid <- distinct(results_summary, cells, genes)
fixed_genes <- grid %>% count(genes) %>% arrange(desc(n)) %>% slice(1) %>% pull(genes)
fixed_cells <- grid %>% count(cells) %>% arrange(desc(n)) %>% slice(1) %>% pull(cells)

make_sweep_plot <- function(df, x, y, ymin, ymax, x_lab, y_lab, title,
                            log_y = FALSE) {
  p <- ggplot(df, aes(x = {{ x }}, y = {{ y }}, colour = algorithm, group = algorithm)) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = {{ ymin }}, ymax = {{ ymax }}),
                  width = 0.05) +
    scale_x_log10(breaks = sort(unique(pull(df, {{ x }})))) +
    labs(x = x_lab, y = y_lab, colour = "Algorithm",
         title = title) +
    theme_bw()
  if (log_y) p <- p + scale_y_log10()
  p
}

cells_df <- results_summary %>% filter(genes == fixed_genes)
genes_df <- results_summary %>% filter(cells == fixed_cells)

p_cells <- make_sweep_plot(cells_df, cells, mean_auc, mean_auc - sd_auc, mean_auc + sd_auc,
                           "Number of cells", "AUROC",
                           sprintf("AUROC vs number of cells (genes = %d)", fixed_genes))
p_genes <- make_sweep_plot(genes_df, genes, mean_auc, mean_auc - sd_auc, mean_auc + sd_auc,
                           "Number of genes", "AUROC",
                           sprintf("AUROC vs number of genes (cells = %d)", fixed_cells))

p_runtime_cells <- make_sweep_plot(cells_df, cells, mean_runtime,
                                   mean_runtime - se_runtime, mean_runtime + se_runtime,
                                   "Number of cells", "Runtime (seconds)",
                                   sprintf("Runtime vs number of cells (genes = %d)", fixed_genes),
                                   log_y = TRUE)
p_runtime_genes <- make_sweep_plot(genes_df, genes, mean_runtime,
                                   mean_runtime - se_runtime, mean_runtime + se_runtime,
                                   "Number of genes", "Runtime (seconds)",
                                   sprintf("Runtime vs number of genes (cells = %d)", fixed_cells),
                                   log_y = TRUE)

fig_dir <- file.path(proj_root, "synthetic/figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, "auroc_vs_cells.pdf"), p_cells, width = 7, height = 5, dpi = 150)
ggsave(file.path(fig_dir, "auroc_vs_genes.pdf"), p_genes, width = 7, height = 5, dpi = 150)
ggsave(file.path(fig_dir, "runtime_vs_cells.pdf"), p_runtime_cells, width = 7, height = 5, dpi = 150)
ggsave(file.path(fig_dir, "runtime_vs_genes.pdf"), p_runtime_genes, width = 7, height = 5, dpi = 150)
cat("Saved AUROC and runtime sweep plots to", fig_dir, "\n")