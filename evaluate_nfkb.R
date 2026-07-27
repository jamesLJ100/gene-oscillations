rm(list=ls())

library(here)
library(dplyr)
library(ROCR)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(readr)
library(msigdbr)
library(BiocParallel)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dorothea)
library(TFEA.ChIP)
library(GenomicRanges)
library(tidyr)

register(SerialParam())  # force single threaded - otherwise fgsea hangs...windows things

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "utils.R"))
source(file.path(proj_root, "gsea_functions.R"))

nfkb_dir <- file.path(proj_root, "data_old/GSE162992/processed")

# Get all .h5 files in directory
input_data_files <- list.files(nfkb_dir, pattern = "\\.h5$", full.names = TRUE)
cat("Found", length(input_data_files), "files to process\n")

# Run the full pipeline for each algorithm. Outputs go to a per-algorithm
# subdirectory (figures/nfkb/<algorithm>/) so the two runs never overwrite
# each other, and the algorithm is also carried in plot titles / filenames.
for (algorithm in c("cyclum", "scPrisma")) {

  cat("\n############################################\n")
  cat("### Algorithm:", algorithm, "\n")
  cat("############################################\n")

  # Create output directory if needed
  figures_dir <- file.path(proj_root, "figures", "nfkb", algorithm)
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

  # Collect all results for comparison
  all_results <- list()

# Loop through all files
for (input_data_file in input_data_files) {
  
  fname <- tools::file_path_sans_ext(basename(input_data_file))
  cat("\n========================================\n")
  cat("Processing:", fname, "\n")
  cat("========================================\n")
  
  # Get model scores
  results_df <- get_model_scores(algorithm, nfkb_dir, fname)
  
  # Convert gene symbols to entrez IDs (MOUSE)
  genes_entrez <- bitr(results_df$symbol,
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = "org.Mm.eg.db")
  
  results_df <- merge(results_df, genes_entrez, by.x = "symbol", by.y = "SYMBOL")
  
  # Build ranked gene list
  geneList <- results_df$score
  names(geneList) <- results_df$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  cat("=== Gene list diagnostics ===\n")
  cat("Genes in geneList:", length(geneList), "\n")
  cat("Score range:", min(geneList), "to", max(geneList), "\n")
  cat("Genes with score = 0:", sum(geneList == 0), "\n")
  
  #temp code-----------------------------------------
  # TEMPORARY: run ChIP-seq (TFEA.ChIP) only for the wt_tnf cohort.
  # Revert to `if (FALSE)` to disable, or remove the guard to run for all cohorts.
  if (fname == "wt_tnf") {
  # Don't re-run the (slow) analysis_from_table + ExperimentHub step; just load
  # the previously saved results and re-plot them.
  cat("\n=== Plotting saved TFEA.ChIP results (analysis skipped) ===\n")

  csv_path <- file.path(figures_dir, paste0(fname, "_tfea_chip_results.csv"))
  if (!file.exists(csv_path)) {
    cat("No saved TFEA.ChIP results at", csv_path, "- run the analysis first. Skipping plot.\n")
    tfea_df <- data.frame(Significant = logical(0))
  } else {
    cat("Loading saved TFEA.ChIP results from", csv_path, "\n")
    tfea_df <- read.csv(csv_path, stringsAsFactors = FALSE)
    tfea_df <- tfea_df[order(tfea_df$pval.adj), ]
    cat(sprintf("Loaded TFEA.ChIP results: %d TFs, %d significant (FDR < 0.05 & ES > 0.1)\n",
                nrow(tfea_df),
                sum(tfea_df$Significant, na.rm = TRUE)))
  }

  if (sum(tfea_df$Significant, na.rm = TRUE) > 0) {

    # The database holds many ChIP experiments per TF (e.g. CTCF x44), and the
    # permutation p-values are floored at 1/nperm so they can't rank the top
    # hits. Collapse to one row per distinct TF (median ES across its
    # experiments, robust to a single lucky cell line) and plot the top N by ES.
    top_n <- 25

    tf_summary <- tfea_df %>%
      dplyr::group_by(TF) %>%
      dplyr::summarise(
        ES_median = median(ES, na.rm = TRUE),
        n_expt    = dplyr::n(),
        n_sig     = sum(Significant, na.rm = TRUE),
        min_FDR   = min(FDR, na.rm = TRUE),
        .groups   = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(ES_median))

    # Save the collapsed per-TF table alongside the plot
    write.csv(tf_summary,
              file.path(figures_dir, paste0(fname, "_tfea_chip_byTF.csv")),
              row.names = FALSE)

    top_tf    <- head(tf_summary, top_n)
    top_tf$TF <- factor(top_tf$TF, levels = rev(top_tf$TF))  # highest ES on top

    p_top <- ggplot(top_tf, aes(x = ES_median, y = TF)) +
      geom_segment(aes(x = 0, xend = ES_median, y = TF, yend = TF),
                   colour = "grey75", linewidth = 0.6) +
      geom_point(aes(size = n_expt), colour = "#08007E") +
      scale_size_continuous(name = "ChIP\nexperiments", range = c(2, 6)) +
      labs(
        title    = paste0("TFEA.ChIP — top ", top_n, " enriched TFs — ", fname),
        subtitle = "one point per TF (median ES across its ChIP experiments)",
        x        = "Median enrichment score (ES)",
        y        = ""
      ) +
      theme_classic() +
      theme(
        axis.text.x   = element_text(size = 14),
        axis.text.y   = element_text(size = 13),
        axis.title    = element_text(size = 16),
        plot.title    = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 12, colour = "grey40"),
        legend.text   = element_text(size = 11),
        legend.title  = element_text(size = 12)
      )
    ggsave(file.path(figures_dir, paste0(fname, "_tfea_chip_topTF.png")),
           plot = p_top, width = 8, height = 9, dpi = 300)

  } else {
    cat("No significant TFs found for", fname, "\n")
  }
  }  # --- end TFEA.ChIP block (wt_tnf only, temporary) ---
  #temp code----------------------------------------

  
  #Run GSEA analyses
  hallmark_results <- run_hallmark_gsea(
    geneList    = geneList,
    species     = "Mus musculus",
    figures_dir = figures_dir,
    fname       = fname,
    algorithm   = algorithm
  )

  #
  # go_results <- run_go_gsea(
  #   geneList    = geneList,
  #   org_db      = org.Mm.eg.db,
  #   figures_dir = figures_dir,
  #   fname       = fname
  # )

  dorothea_results <- run_dorothea_gsea(
    geneList    = geneList,
    species     = "mouse",
    org_db      = org.Mm.eg.db,
    figures_dir = figures_dir,
    fname       = fname
  )
  
  # tfea_results <- run_tfea_chip(
  #   geneList    = geneList,
  #   figures_dir = figures_dir,
  #   fname       = fname,
  #   species = "mouse"
  # )
  
  # Store scores and hallmark results separately
  all_results[[fname]] <- list(
    scores   = data.frame(
      cluster = fname,
      score   = results_df$score,
      symbol  = results_df$symbol
    ),
    hallmark = if (exists("hallmark_results") && nrow(hallmark_results) > 0) hallmark_results else NULL
  )
}


# ---- Hallmark Heatmap -------------------------------------------------------
cat("\n=== Creating Hallmark heatmap ===\n")

hallmark_combined <- bind_rows(
  lapply(names(all_results), function(fname) {
    hr <- all_results[[fname]]$hallmark
    if (!is.null(hr) && nrow(hr) > 0) {
      hr %>%
        dplyr::select(ID, NES, p.adjust) %>%
        mutate(cluster = fname)
    }
  })
)

if (nrow(hallmark_combined) > 0) {
  
  # only keep pathways significant in at least one cluster
  sig_pathways <- hallmark_combined %>%
    filter(p.adjust < 0.05) %>%
    pull(ID) %>%
    unique()
  
  heatmap_df <- hallmark_combined %>%
    filter(ID %in% sig_pathways) %>%
    mutate(
      NES_plot = ifelse(p.adjust < 0.05, NES, NA),
      ID       = gsub("HALLMARK_", "", ID)
    )
  
  # NES values are all positive here (GSEA runs with scoreType = "pos"), so a
  # diverging scale centred at 0 pushes every tile into one colour. Centre the
  # diverging scale on the median NES so both tails get colour and contrast.
  nes_midpoint <- median(heatmap_df$NES_plot, na.rm = TRUE)
  p_heatmap <- ggplot(heatmap_df, aes(x = cluster, y = ID, fill = NES_plot)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(
      low      = "#08007E",
      mid      = "white",
      high     = "#FF0187",
      midpoint = nes_midpoint,
      na.value = "grey90",
      name     = "NES"
    ) +
    labs(
      title = paste(algorithm, "- Hallmark GSEA across clusters"),
      x     = "",
      y     = ""
    ) +
    theme_bw(base_size = 18) +
    theme(
      axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 14),
      plot.title   = element_text(size = 22, face = "bold"),
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 16),
      panel.grid   = element_blank()
    )
  
  ggsave(file.path(figures_dir, paste0(algorithm, "_hallmark_heatmap.png")),
         plot   = p_heatmap,
         width  = 4 + length(unique(heatmap_df$cluster)) * 1.2,
         height = 3 + length(sig_pathways) * 0.3,
         dpi    = 300)
  
  cat("Hallmark heatmap saved.\n")
  
} else {
  cat("No significant hallmark pathways found across any cluster.\n")
}


# ---- Score Distribution Comparison Plot -------------------------------------
cat("\n=== Creating score distribution comparison plot ===\n")

# Combine scores across all files
df <- bind_rows(lapply(all_results, `[[`, "scores"))

# Check what clusters we actually have
cat("\nActual clusters found:\n")
print(unique(df$cluster))

cluster_order <- sort(unique(df$cluster))
df$cluster <- factor(df$cluster, levels = cluster_order)

# --- Theme ---
theme_common <- function(base_size = 16) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.ticks.length = unit(-0.25, "cm"),
      legend.position   = "none",
      strip.background  = element_rect(fill = "transparent", colour = NA),
      axis.text.x       = element_text(angle = 45, hjust = 1)
    )
}

# Find reference group (wt_tnf)
reference_group <- grep("wt.*tnf|tnf.*wt", cluster_order, ignore.case = TRUE, value = TRUE)

if (length(reference_group) == 0) {
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
    title = paste(algorithm, "- Score comparison"),
    x     = "",
    y     = paste(algorithm, "Score")
  ) +
  theme_common(14)

ggsave(file.path(figures_dir, paste0(algorithm, "_score_comparison_statistical.png")),
       plot = p1, width = 8, height = 6, dpi = 300)

}  # end for (algorithm in ...)