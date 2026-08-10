rm(list=ls())
library(here)
library(dplyr)
library(ROCR)
library(PRROC)
library(ggplot2)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/hdfrw.R"))
source(file.path(proj_root, "synthetic/utils/dyngen_utils.R"))
source(file.path(proj_root, "common/utils.R"))


algorithms <- c("scPrisma", "cyclum", "oscope")
dyngen_dir <- file.path(proj_root, "synthetic/data/dyngen")

# Each (n_cells, n_genes) combination has its own subdirectory (see
# generate_datasets.R) rather than one flat directory - "gridsearch" is excluded
# automatically since it doesn't match the c<n_cells>g<n_genes> naming pattern.
combos  <- list_combo_dirs(dyngen_dir)
results <- list()
gene_level_results <- list()

if (nrow(combos) == 0) {
  stop("No evaluation-set combinations found under ", dyngen_dir,
       " - run generate_datasets.R first.")
}

# Each algorithm writes a runtimes.csv (columns: file, runtime_seconds, ...) into
# its own per-combination output dir (dyngen/c<n_cells>g<n_genes>/<algorithm>/).
# Combine them all into a single `runtimes` table tagged with the algorithm,
# averaging over any repeated runs (run_number) of the same file.
runtimes <- do.call(rbind, lapply(seq_len(nrow(combos)), function(i) {
  do.call(rbind, lapply(algorithms, function(algorithm) {
    rt_path <- file.path(combos$path[i], algorithm, "runtimes.csv")
    if (!file.exists(rt_path)) return(NULL)
    rt <- read.csv(rt_path, stringsAsFactors = FALSE)
    rt <- rt[!is.na(rt$runtime_seconds), c("file", "runtime_seconds")]
    if (nrow(rt) == 0) return(NULL)
    # Average across repeated runs so each (file, algorithm) has one runtime
    rt <- aggregate(runtime_seconds ~ file, data = rt, FUN = mean)
    rt$algorithm <- algorithm
    rt
  }))
}))

for (combo_i in seq_len(nrow(combos))) {
  combo_dir <- combos$path[combo_i]
  n_cells   <- combos$n_cells[combo_i]
  n_genes   <- combos$n_genes[combo_i]

  expr_files <- list.files(combo_dir, pattern = "\\.h5$", full.names = TRUE)

  for (expr_file in expr_files) {
    fname <- tools::file_path_sans_ext(basename(expr_file))

    # Load simulation file once per expr_file
    sim_file <- file.path(combo_dir, paste0(fname, "_sim.rds"))
    if (!file.exists(sim_file)) {
      cat("Skipping", fname, "- no corresponding sim RDS found\n")
      next
    }
    sim <- readRDS(sim_file)

    gene_module_assignments <- propagate_module_assignments(sim, context = fname)

    # Only evaluate on genes at most one hop from a TF (hops <= 1); this also
    # drops cycle genes with NA hops, since NA <= 1 is NA (filtered out).
    cycling_genes <- gene_module_assignments %>%
      dplyr::filter(module %in% c("B", "C", "D"), hops <= 1) %>%
      pull(gene)

    gene_symbols <- rownames(sim$counts)

    eval_genes   <- gene_module_assignments %>% dplyr::filter(hops <= 1) %>% pull(gene)
    eval_symbols <- gene_symbols[gene_symbols %in% eval_genes]

    label_df <- data.frame(symbol = eval_symbols, label = 0L)
    label_df$label[eval_symbols %in% cycling_genes] <- 1L

    # dyngen's own simulated per-gene kinetics, for the mispredicted-gene analysis
    # near the end of this script (are algorithms more likely to miss genes with
    # certain production/degradation rates?).
    kinetics_df <- sim$model$feature_info %>%
      dplyr::select(symbol = feature_id, transcription_rate, translation_rate,
                     mrna_decay_rate, protein_decay_rate)

    # Now process each algorithm for this file
    for (algorithm in algorithms) {
      cat("Processing:", fname, "with", algorithm, "\n")

      algorithm_results_df <- tryCatch(
        get_model_scores(algorithm, combo_dir, fname),
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

      # Threshold-dependent metrics (accuracy, F1, precision, recall) need a cutoff -
      # see compute_threshold_metrics() for how binary vs continuous scores differ.
      threshold_metrics <- compute_threshold_metrics(merge_df$score, merge_df$label)
      acc_val  <- threshold_metrics$accuracy
      prec_val <- threshold_metrics$precision
      rec_val  <- threshold_metrics$recall
      f1_val   <- threshold_metrics$f1

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

      gene_level_results[[result_key]] <- data.frame(
        algorithm = algorithm,
        file      = fname,
        cells     = n_cells,
        genes     = n_genes,
        symbol    = merge_df$symbol,
        label     = merge_df$label,
        predicted = threshold_metrics$predicted,
        stringsAsFactors = FALSE
      ) %>% merge(kinetics_df, by = "symbol", all.x = TRUE)
    }
  }
}

results_df <- do.call(rbind, lapply(results, as.data.frame))
rownames(results_df) <- NULL
print(results_df)

results_summary <- results_df %>%
  dplyr::group_by(algorithm, cells, genes) %>%
  dplyr::summarise(
    mean_auc        = mean(auc,      na.rm = TRUE),
    sd_auc          = sd(auc,        na.rm = TRUE),
    se_auc          = sd(auc,        na.rm = TRUE) / sqrt(dplyr::n()),
    mean_auprc      = mean(auprc,    na.rm = TRUE),
    sd_auprc        = sd(auprc,      na.rm = TRUE),
    se_auprc        = sd(auprc,      na.rm = TRUE) / sqrt(dplyr::n()),
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
    se_runtime      = sd(runtime,   na.rm = TRUE) / sqrt(dplyr::n()),
    n               = dplyr::n(),
    .groups         = "drop"
  )
print(results_summary)

# ---- Metric sweeps: vs number of cells and vs number of genes ---------------
# The datasets form two sweeps sharing a common anchor (c1000g204): one varies
# cells with genes held fixed, the other varies genes with cells held fixed.
# Detect the held-fixed value automatically as the one covering the most points.
grid <- distinct(results_summary, cells, genes)
# dplyr::-qualified throughout this block: count/arrange/desc/slice all have
# same-named exports in other packages this codebase loads elsewhere in the same
# pipeline (matrixStats::count, IRanges::{slice,desc}, clusterProfiler re-exports
# several dplyr verbs) - any of them loaded earlier in the session masks dplyr's.
fixed_genes <- grid %>% dplyr::count(genes) %>% dplyr::arrange(dplyr::desc(n)) %>% dplyr::slice(1) %>% pull(genes)
fixed_cells <- grid %>% dplyr::count(cells) %>% dplyr::arrange(dplyr::desc(n)) %>% dplyr::slice(1) %>% pull(cells)

# x_col/mean_col/error_col are column names (strings) rather than embraced
# expressions so this can be driven from the metrics table below in a loop.
make_sweep_plot <- function(df, x_col, mean_col, error_col, x_lab, y_lab, title,
                            log_y = FALSE) {
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[mean_col]],
                       colour = algorithm, group = algorithm)) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = .data[[mean_col]] - .data[[error_col]],
                       ymax = .data[[mean_col]] + .data[[error_col]]),
                  width = 0.05) +
    scale_x_log10(breaks = sort(unique(df[[x_col]]))) +
    labs(x = x_lab, y = y_lab, colour = "Algorithm",
         title = title) +
    theme_bw()
  if (log_y) p <- p + scale_y_log10()
  p
}

cells_df <- results_summary %>% dplyr::filter(genes == fixed_genes)
genes_df <- results_summary %>% dplyr::filter(cells == fixed_cells)

# One sweep pair (vs cells, vs genes) per metric, all following the same
# mean/error-bar convention as the runtime plots below.
metrics <- list(
  list(file = "auroc",     mean = "mean_auc",       error = "sd_auc",       label = "AUROC"),
  list(file = "auprc",     mean = "mean_auprc",     error = "sd_auprc",     label = "AUPRC"),
  list(file = "accuracy",  mean = "mean_accuracy",  error = "sd_accuracy",  label = "Accuracy"),
  list(file = "f1",        mean = "mean_f1",        error = "sd_f1",        label = "F1"),
  list(file = "precision", mean = "mean_precision", error = "sd_precision", label = "Precision"),
  list(file = "recall",    mean = "mean_recall",    error = "sd_recall",    label = "Recall")
)

fig_dir <- file.path(proj_root, "synthetic/figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

for (m in metrics) {
  p_cells <- make_sweep_plot(cells_df, "cells", m$mean, m$error,
                             "Number of cells", m$label,
                             sprintf("%s vs number of cells (genes = %d)", m$label, fixed_genes))
  p_genes <- make_sweep_plot(genes_df, "genes", m$mean, m$error,
                             "Number of genes", m$label,
                             sprintf("%s vs number of genes (cells = %d)", m$label, fixed_cells))
  ggsave(file.path(fig_dir, sprintf("%s_vs_cells.pdf", m$file)), p_cells, width = 7, height = 5, dpi = 150)
  ggsave(file.path(fig_dir, sprintf("%s_vs_genes.pdf", m$file)), p_genes, width = 7, height = 5, dpi = 150)
}

p_runtime_cells <- make_sweep_plot(cells_df, "cells", "mean_runtime", "se_runtime",
                                   "Number of cells", "Runtime (seconds)",
                                   sprintf("Runtime vs number of cells (genes = %d)", fixed_genes),
                                   log_y = TRUE)
p_runtime_genes <- make_sweep_plot(genes_df, "genes", "mean_runtime", "se_runtime",
                                   "Number of genes", "Runtime (seconds)",
                                   sprintf("Runtime vs number of genes (cells = %d)", fixed_cells),
                                   log_y = TRUE)
ggsave(file.path(fig_dir, "runtime_vs_cells.pdf"), p_runtime_cells, width = 7, height = 5, dpi = 150)
ggsave(file.path(fig_dir, "runtime_vs_genes.pdf"), p_runtime_genes, width = 7, height = 5, dpi = 150)

cat("Saved metric and runtime sweep plots to", fig_dir, "\n")

# ---- Kinetics analysis: do mispredicted genes have distinctive production/
# degradation rates? -----------------------------------------------------------
# For each true cycling gene (label == 1), TP = correctly detected, FN = missed.
# For each true non-cycling gene (label == 0), TN = correctly rejected, FP =
# wrongly flagged. Comparing rate distributions across these tells us whether an
# algorithm's failures are concentrated among genes with particular kinetics
# (e.g. fast-degrading genes), rather than being uniform across the gene panel.
gene_level_df <- do.call(rbind, gene_level_results)
rownames(gene_level_df) <- NULL
gene_level_df$outcome <- dplyr::case_when(
  gene_level_df$label == 1 & gene_level_df$predicted == 1 ~ "TP",
  gene_level_df$label == 1 & gene_level_df$predicted == 0 ~ "FN",
  gene_level_df$label == 0 & gene_level_df$predicted == 1 ~ "FP",
  gene_level_df$label == 0 & gene_level_df$predicted == 0 ~ "TN",
  TRUE ~ NA_character_
)
gene_level_df$outcome <- factor(gene_level_df$outcome, levels = c("TP", "FN", "TN", "FP"))

rate_params <- list(
  list(col = "transcription_rate", label = "Transcription rate"),
  list(col = "translation_rate",   label = "Translation rate"),
  list(col = "mrna_decay_rate",    label = "mRNA decay rate"),
  list(col = "protein_decay_rate", label = "Protein decay rate")
)

# Wilcoxon rank-sum test per (algorithm, rate): does this rate differ between
# correctly- and incorrectly-classified genes, among the genes that share the
# same ground-truth label (TP vs FN among true positives, TN vs FP among true
# negatives)?
kinetics_tests <- do.call(rbind, lapply(algorithms, function(algo) {
  algo_df <- gene_level_df %>% dplyr::filter(algorithm == algo)
  do.call(rbind, lapply(rate_params, function(rp) {
    tp <- algo_df[[rp$col]][algo_df$outcome == "TP"]
    fn <- algo_df[[rp$col]][algo_df$outcome == "FN"]
    tn <- algo_df[[rp$col]][algo_df$outcome == "TN"]
    fp <- algo_df[[rp$col]][algo_df$outcome == "FP"]
    p_tp_fn <- tryCatch(wilcox.test(tp, fn)$p.value, error = function(e) NA_real_)
    p_tn_fp <- tryCatch(wilcox.test(tn, fp)$p.value, error = function(e) NA_real_)
    data.frame(algorithm = algo, rate = rp$label,
               p_TP_vs_FN = p_tp_fn, n_TP = length(tp), n_FN = length(fn),
               p_TN_vs_FP = p_tn_fp, n_TN = length(tn), n_FP = length(fp))
  }))
}))
cat("\nKinetics: rate differences between correctly- and incorrectly-classified genes\n")
cat("(Wilcoxon rank-sum p-values; small p means that rate differs by outcome)\n")
print(kinetics_tests)
write.csv(kinetics_tests, file.path(fig_dir, "kinetics_test_pvalues.csv"), row.names = FALSE)

make_kinetics_plot <- function(df, rate_col, y_lab) {
  ggplot(df, aes(x = outcome, y = .data[[rate_col]], fill = outcome)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~algorithm) +
    scale_y_log10() +
    labs(x = "Prediction outcome", y = y_lab,
         title = sprintf("%s by prediction outcome", y_lab)) +
    theme_bw() +
    theme(legend.position = "none")
}

for (rp in rate_params) {
  p <- make_kinetics_plot(gene_level_df, rp$col, rp$label)
  ggsave(file.path(fig_dir, sprintf("kinetics_%s.pdf", rp$col)), p, width = 8, height = 5, dpi = 150)
}

cat("Saved kinetics analysis (plots + p-values) to", fig_dir, "\n")