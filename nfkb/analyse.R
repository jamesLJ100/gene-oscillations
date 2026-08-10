rm(list=ls())

library(here)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(msigdbr)
library(BiocParallel)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dorothea)
library(TFEA.ChIP)

register(SerialParam())  # force single threaded - otherwise fgsea hangs...windows things

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/utils.R"))
source(file.path(proj_root, "common/gsea_functions.R"))
source(file.path(proj_root, "common/plotting_utils.R"))

nfkb_dir     <- file.path(proj_root, "nfkb/data/GSE162992/processed")
results_root <- file.path(proj_root, "nfkb/results")

# Get all .h5 files in directory
input_data_files <- list.files(nfkb_dir, pattern = "\\.h5$", full.names = TRUE)
cat("Found", length(input_data_files), "files to process\n")

# Built once and passed into every run_pathway_analyses() call below - species is
# fixed for this whole script, so these would otherwise be rebuilt on every call.
hallmark_geneset   <- build_hallmark_geneset("Mus musculus")
dorothea_term2gene <- build_dorothea_term2gene("mouse", org.Mm.eg.db)

# Run the full pipeline for each algorithm. Outputs go to a per-algorithm
# subdirectory (nfkb/figures/<algorithm>/) so the two runs never overwrite
# each other, and the algorithm is also carried in plot titles / filenames.
for (algorithm in c("cyclum", "scPrisma")) {

  cat("\n############################################\n")
  cat("### Algorithm:", algorithm, "\n")
  cat("############################################\n")

  # Create output directory if needed
  figures_dir <- file.path(proj_root, "nfkb/figures", algorithm)
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

  results_dir <- file.path(results_root, algorithm)

  # Collect all results for comparison
  all_results <- list()

  # Loop through all files (each is one WT/IkBaMM x stimulus condition)
  for (input_data_file in input_data_files) {

    fname <- tools::file_path_sans_ext(basename(input_data_file))
    cat("\n========================================\n")
    cat("Processing:", fname, "\n")
    cat("========================================\n")

    results_df <- load_model_scores(algorithm, results_dir, input_data_file, fname)
    if (is.null(results_df)) {
      stop("No ", algorithm, " results found for ", fname, " in ", results_dir)
    }

    geneList <- build_ranked_gene_list(results_df, org_db = org.Mm.eg.db)

    pathway_results <- run_pathway_analyses(
      geneList           = geneList,
      species            = "mouse",
      figures_dir        = figures_dir,
      fname              = fname,
      algorithm          = algorithm,
      hallmark           = hallmark_geneset,
      dorothea_term2gene = dorothea_term2gene
    )

    # Store scores and hallmark results separately (hallmark feeds the
    # cross-condition heatmap below; scores feed the violin comparison)
    all_results[[fname]] <- list(
      scores   = data.frame(
        cluster = fname,
        score   = results_df$score,
        symbol  = results_df$symbol
      ),
      hallmark = pathway_results$hallmark
    )
  }

  # ---- Hallmark Heatmap -------------------------------------------------------
  plot_hallmark_heatmap(
    hallmark_by_group = lapply(all_results, `[[`, "hallmark"),
    figures_dir       = figures_dir,
    algorithm         = algorithm
  )

  # ---- Pairwise Hallmark + Oscillation Score Comparison Plots ------------------
  # wt_tnf vs ss_tnf: same stimulus, does the TNF response differ by genotype?
  # wt_none vs wt_tnf: same genotype, what does TNF stimulation itself change?
  # A list rather than separate hardcoded calls so adding another pair (e.g. an
  # lps or pic comparison) is a one-line addition here, not a copy-pasted block.
  comparisons <- list(
    c(a = "wt_tnf",  b = "ss_tnf"),
    c(a = "wt_none", b = "wt_tnf")
  )

  for (cmp in comparisons) {
    # Pathway-level: one point per Hallmark pathway
    plot_hallmark_comparison(
      hallmark_a  = all_results[[cmp[["a"]]]]$hallmark,
      hallmark_b  = all_results[[cmp[["b"]]]]$hallmark,
      label_a     = cmp[["a"]],
      label_b     = cmp[["b"]],
      figures_dir = figures_dir,
      algorithm   = algorithm
    )
    # Gene-level analogue: one point per gene rather than per pathway
    plot_score_comparison(
      scores_a    = all_results[[cmp[["a"]]]]$scores,
      scores_b    = all_results[[cmp[["b"]]]]$scores,
      label_a     = cmp[["a"]],
      label_b     = cmp[["b"]],
      figures_dir = figures_dir,
      algorithm   = algorithm
    )
  }

  # ---- Score Distribution Comparison Plot --------------------------------------
  cat("\n=== Creating score distribution comparison plot ===\n")

  # Combine scores across all files
  df <- bind_rows(lapply(all_results, `[[`, "scores"))

  # Check what clusters we actually have
  cat("\nActual clusters found:\n")
  print(unique(df$cluster))

  cluster_order <- sort(unique(df$cluster))
  df$cluster <- factor(df$cluster, levels = cluster_order)

  # Find reference group (wt_tnf) - same condition name used for the pairwise
  # comparisons above
  if ("wt_tnf" %in% cluster_order) {
    reference_group <- "wt_tnf"
  } else {
    cat("\nWarning: Could not find 'wt_tnf' group. Using first group as reference.\n")
    reference_group <- cluster_order[1]
  }

  cat("\nUsing reference group:", reference_group, "\n")

  other_groups     <- setdiff(cluster_order, reference_group)
  comparisons_list <- lapply(other_groups, function(x) c(reference_group, x))

  color_values                  <- setNames(rep("#08007E", length(cluster_order)), cluster_order)
  color_values[reference_group] <- "#FF0187"

  p1 <- ggplot(df, aes(x = cluster, y = score, fill = cluster)) +
    geom_violin(alpha = 0.9, color = "white") +
    geom_boxplot(width = 0.15, outlier.colour = NA, fill = "white") +
    stat_compare_means(
      comparisons     = comparisons_list,
      label           = "p.signif",
      method          = "wilcox.test",
      p.adjust.method = "BH",
      size            = 5
    ) +
    scale_fill_manual(values = color_values) +
    labs(
      x = "",
      y = paste(algorithm, "Score")
    ) +
    theme_common(14, angle_x_labels = TRUE)

  ggsave(file.path(figures_dir, paste0(algorithm, "_score_comparison_statistical.pdf")),
         plot = p1, width = 8, height = 6, dpi = 300)

}  # end for (algorithm in ...)
