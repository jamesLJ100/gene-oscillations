rm(list=ls())
library(here)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr)
library(BiocParallel)
library(dorothea)
library(TFEA.ChIP)
library(GenomicRanges)

register(SerialParam())

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/hdfrw.R"))
source(file.path(proj_root, "common/utils.R"))
source(file.path(proj_root, "common/gsea_functions.R"))
source(file.path(proj_root, "common/plotting_utils.R"))

# Algorithms to evaluate. Display name is used in plot titles / axis labels,
# and the key is used as a suffix on output filenames.
algo_display <- c(cyclum = "Cyclum", scPrisma = "scPrisma")

# ============================================================================
# Per-dataset configuration. `ref_cluster` is the "reference" cluster that
# Plot 2/3/4 compare everything else against (pPSM for mme95/mESC, since
# that's where the segmentation clock is actively oscillating; the closest
# analogue for hIPSC is aPSM, as it has no pPSM cluster). `cluster_order` uses
# the sanitised names that appear in the saved .h5/.csv filenames (spaces and
# hyphens become underscores).
# ============================================================================
dataset_configs <- list(
  mme95 = list(
    label         = "mmE9.5 PM",
    counts_dir    = file.path(proj_root, "segmentation_clock/data/mme95/counts/purity_filtered"),
    cyclum_dir    = file.path(proj_root, "segmentation_clock/results/cyclum/mme95"),
    scPrisma_dir  = file.path(proj_root, "segmentation_clock/results/scPrisma/mme95"),
    figures_dir   = file.path(proj_root, "segmentation_clock/figures/mme95"),
    file_prefix   = "mmE95_PM_",
    species       = "mouse",
    ref_cluster   = "pPSM",
    cluster_order = c("NMP", "MPC", "pPSM", "aPSM", "scSOM", "dmSOM")
  ),
  mESC = list(
    label         = "mESC",
    counts_dir    = file.path(proj_root, "segmentation_clock/data/mmESC/counts/purity_filtered"),
    cyclum_dir    = file.path(proj_root, "segmentation_clock/results/cyclum/mESC"),
    scPrisma_dir  = file.path(proj_root, "segmentation_clock/results/scPrisma/mESC"),
    figures_dir   = file.path(proj_root, "segmentation_clock/figures/mESC"),
    file_prefix   = "mmESC_",
    species       = "mouse",
    ref_cluster   = "d4_5_pPSM",
    cluster_order = c("d0_ESC", "d2_Epi", "d3_NMP", "d4_5_Col4", "d4_5_Krt", "d4_5_pPSM")
  ),
  hIPSC = list(
    label         = "hIPSC",
    counts_dir    = file.path(proj_root, "segmentation_clock/data/hIPSC/counts/purity_filtered"),
    cyclum_dir    = file.path(proj_root, "segmentation_clock/results/cyclum/hIPSC"),
    scPrisma_dir  = file.path(proj_root, "segmentation_clock/results/scPrisma/hIPSC"),
    figures_dir   = file.path(proj_root, "segmentation_clock/figures/hIPSC"),
    file_prefix   = "hIPSC_",
    species       = "human",
    ref_cluster   = "d3_4_aPSM",
    cluster_order = c("d0_iPSC", "d1_NMP", "d2_MPC", "d3_4_aPSM")
  )
)

reference_lists <- list(
  mouse = list(
    oscillating  = readRDS("segmentation_clock/data/oscillating_mouse.rds"),
    housekeeping = readRDS("segmentation_clock/data/housekeeping_mouse.rds"),
    cellcycle    = readRDS("segmentation_clock/data/cellcycle_mouse.rds")
  ),
  human = list(
    oscillating  = readRDS("segmentation_clock/data/oscillating_human.rds"),
    housekeeping = readRDS("segmentation_clock/data/housekeeping_human.rds"),
    cellcycle    = readRDS("segmentation_clock/data/cellcycle_human.rds")
  )
)

# ============================================================================
# Load per-cluster cycling scores for one algorithm within one dataset. Both
# run_cyclum.py and run_scPrisma.py write output using the same basename as
# the input count file, just in their own --output_dir:
#   cyclum   -> <name>.h5 weights (in cyclum_dir) + the count matrix
#   scPrisma -> <name>.csv (in scPrisma_dir)
# ============================================================================
load_scores <- function(algorithm, counts_files, file_prefix, cyclum_dir, scPrisma_dir) {
  results <- list()
  cluster_pattern <- paste0("^", file_prefix, "(.+)_counts$")
  results_dir <- if (algorithm == "cyclum") cyclum_dir else scPrisma_dir

  for (expr_file in counts_files) {
    file_name <- tools::file_path_sans_ext(basename(expr_file))
    cluster   <- gsub(cluster_pattern, "\\1", file_name)

    cat("Processing:", cluster, "with", algorithm, "\n")

    scores <- load_model_scores(algorithm, results_dir, expr_file, file_name)
    if (is.null(scores)) {
      cat("  no", algorithm, "results for:", cluster, "-", results_dir, "\n")
      next
    }

    scores$cluster <- cluster
    results[[cluster]] <- scores
  }

  if (length(results) == 0) return(data.frame())

  df <- bind_rows(results)
  # get_cyclum_scores returns "symbol" column — rename to gene
  if ("symbol" %in% colnames(df)) df <- dplyr::rename(df, gene = symbol)
  df
}

# ============================================================================
# Compute mean raw-count expression per gene, per cluster, from the same counts
# .h5 files load_scores() reads. Needed independently of load_scores() because
# scPrisma's own score loading just reads scPrisma's output CSV, not the input
# counts matrix, so the counts data isn't otherwise available for the
# score-vs-expression plots below.
# ============================================================================
compute_mean_expression <- function(counts_files, file_prefix) {
  results <- list()
  cluster_pattern <- paste0("^", file_prefix, "(.+)_counts$")

  for (expr_file in counts_files) {
    file_name <- tools::file_path_sans_ext(basename(expr_file))
    cluster   <- gsub(cluster_pattern, "\\1", file_name)

    counts_mat <- hdf2mat(expr_file)  # genes x cells

    results[[cluster]] <- data.frame(
      gene      = rownames(counts_mat),
      mean_expr = rowMeans(counts_mat),
      cluster   = cluster,
      row.names = NULL
    )
  }

  bind_rows(results)
}

# ============================================================================
# Run Hallmark/GO/DoRothEA/TFEA.ChIP pathway analysis separately for each cell
# cluster/type in `df` (not pooled across clusters - a segmentation-clock gene
# oscillating specifically in, say, pPSM cells is a different biological signal
# from it oscillating everywhere, and pooling would wash that out). Figures go
# to <config$figures_dir>/<algorithm>/ so the two algorithms don't collide, and
# a Hallmark heatmap comparing all of this dataset's clusters is saved there too.
# ============================================================================
run_pathway_analysis_per_cluster <- function(df, config, algorithm) {

  figures_dir <- file.path(config$figures_dir, algorithm)
  org_db      <- resolve_org_db(config$species)

  hallmark_by_cluster <- list()

  for (cluster_name in levels(droplevels(df$cluster))) {
    cluster_df <- df[df$cluster == cluster_name, c("gene", "score")]
    if (nrow(cluster_df) == 0) next

    cat("\n\n----------------------------------------\n")
    cat("Pathway analysis:", config$label, "-", algorithm, "-", cluster_name, "\n")
    cat("----------------------------------------\n")

    geneList <- build_ranked_gene_list(cluster_df, org_db = org_db, symbol_col = "gene")

    pathway_results <- run_pathway_analyses(
      geneList    = geneList,
      species     = config$species,
      figures_dir = figures_dir,
      fname       = cluster_name,
      algorithm   = algorithm
    )

    hallmark_by_cluster[[cluster_name]] <- pathway_results$hallmark
  }

  plot_hallmark_heatmap(
    hallmark_by_group = hallmark_by_cluster,
    figures_dir       = figures_dir,
    algorithm         = algorithm
  )

  invisible(hallmark_by_cluster)
}

# ============================================================================
# Run the full evaluation (load -> debug -> annotate -> plot -> save) for a
# single algorithm within a single dataset. `score_label` labels the y-axis;
# `suffix` is appended to output filenames so the two algorithms don't
# overwrite each other. Returns NULL (with a message) if no results exist yet
# for this dataset/algorithm combination.
# ============================================================================
run_evaluation <- function(algorithm, config, refs) {
  display     <- algo_display[[algorithm]]
  suffix      <- algorithm
  ref_cluster <- config$ref_cluster
  oscillating  <- refs$oscillating
  housekeeping <- refs$housekeeping
  cellcycle    <- refs$cellcycle

  counts_files <- list.files(config$counts_dir, pattern = "\\.h5$", full.names = TRUE)

  cat("\n\n########################################\n")
  cat("### ", config$label, " - ", display, " evaluation\n", sep = "")
  cat("########################################\n")

  df <- load_scores(algorithm, counts_files, config$file_prefix, config$cyclum_dir, config$scPrisma_dir)

  if (nrow(df) == 0) {
    cat("No", display, "results found for", config$label, "- skipping (run the pipeline for this dataset first).\n")
    return(NULL)
  }

  cat("Columns in df:", colnames(df), "\n")
  cat("Total rows:", nrow(df), "\n")

  # ========================================
  # === DEBUGGING SECTION STARTS HERE ===
  # ========================================

  cat("\n=============== DEBUGGING (", config$label, ", ", display, ") ===============\n", sep = "")

  # 1. Check unique genes overall
  unique_genes <- unique(df$gene)
  cat("\nTotal unique genes across all data:", length(unique_genes), "\n")

  # 2. Check genes per cluster
  cat("\nGenes per cluster:\n")
  genes_per_cluster <- df %>%
    group_by(cluster) %>%
    summarise(
      n_rows = n(),
      n_unique_genes = n_distinct(gene)
    )
  print(genes_per_cluster)

  # 3. Check if all clusters have the same genes
  cat("\nChecking if all clusters share the same gene set...\n")
  cluster_gene_sets <- df %>%
    group_by(cluster) %>%
    summarise(genes = list(sort(unique(gene))), .groups = 'drop')

  all_same <- TRUE
  if (nrow(cluster_gene_sets) > 1) {
    for (i in 2:nrow(cluster_gene_sets)) {
      if (!setequal(cluster_gene_sets$genes[[1]], cluster_gene_sets$genes[[i]])) {
        all_same <- FALSE
        cat("  Cluster", cluster_gene_sets$cluster[[i]], "has different genes!\n")
      }
    }
  }
  if (all_same) {
    cat("  ✓ All clusters have the same gene set\n")
  }

  # 4. Sample genes to check format
  cat("\nFirst 20 genes in your data:\n")
  print(head(unique_genes, 20))

  cat("\nFirst 20 oscillating genes in reference:\n")
  print(head(oscillating$SYMBOL, 20))

  # 5. Check reference list sizes
  cat("\nReference list sizes:\n")
  cat("  Oscillating genes:", nrow(oscillating), "\n")
  cat("  Housekeeping genes:", nrow(housekeeping), "\n")
  cat("  Cell cycle genes:", nrow(cellcycle), "\n")

  # 6. Check for duplicates within each cluster
  cat("\nChecking for duplicate genes within clusters:\n")
  duplicates <- df %>%
    group_by(cluster, gene) %>%
    summarise(n = n(), .groups = 'drop') %>%
    filter(n > 1)

  if (nrow(duplicates) > 0) {
    cat("  ⚠ Found", nrow(duplicates), "gene duplicates within clusters!\n")
    cat("  Examples:\n")
    print(head(duplicates, 10))
  } else {
    cat("  ✓ No duplicates within clusters\n")
  }

  # 7. Test case-insensitive matching on UNIQUE genes only
  cat("\nMatching UNIQUE genes to reference lists (case-insensitive):\n")
  unique_genes_lower <- tolower(unique_genes)

  n_osc_match <- sum(unique_genes_lower %in% tolower(oscillating$SYMBOL))
  n_hk_match  <- sum(unique_genes_lower %in% tolower(housekeeping$SYMBOL))
  n_cc_match  <- sum(unique_genes_lower %in% tolower(cellcycle$SYMBOL))

  cat("  Oscillating matches:", n_osc_match, "\n")
  cat("  Housekeeping matches:", n_hk_match, "\n")
  cat("  Cell cycle matches:", n_cc_match, "\n")
  cat("  Other (unmatched):", length(unique_genes) - n_osc_match - n_hk_match - n_cc_match, "\n")

  cat("\n=============== END DEBUGGING (", config$label, ", ", display, ") ===============\n\n", sep = "")

  # ========================================
  # === ANNOTATE + PLOT ===
  # ========================================

  # --- Case-insensitive gene category annotation ---
  df$gene_lower <- tolower(df$gene)

  df$dynamic <- "Other"
  df$dynamic[df$gene_lower %in% tolower(oscillating$SYMBOL)]  <- "Oscillating"
  df$dynamic[df$gene_lower %in% tolower(housekeeping$SYMBOL)] <- "Housekeeping"
  df$dynamic[df$gene_lower %in% tolower(cellcycle$SYMBOL)]    <- "Cell Cycle"
  df$dynamic <- factor(df$dynamic,
                       levels = c("Oscillating", "Housekeeping", "Cell Cycle", "Other"))

  cat("\nGene category counts (TOTAL ROWS, not unique genes):\n")
  print(table(df$dynamic))

  cat("\nGene category counts by cluster:\n")
  print(table(df$cluster, df$dynamic))

  # --- Cluster ordering ---
  cluster_order <- config$cluster_order
  df$cluster <- factor(df$cluster,
                       levels = cluster_order[cluster_order %in% unique(df$cluster)])

  # --- Score summary ---
  cat("\nScore summary by cluster:\n")
  print(df |> group_by(cluster) |> summarise(
    n_genes      = n(),
    mean_score   = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    n_nonzero    = sum(score > 0)
  ))

  score_label <- paste(display, "Cycling Score")

  # --- Plot 1: all genes, all clusters ---
  p1 <- ggplot(df, aes(x = cluster, y = score, fill = cluster)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "white") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA) +
    labs(x = "", y = score_label) +
    theme_common(14)

  # --- Plot 2: oscillating genes only, ref_cluster as reference ---
  plot_osc <- df[df$dynamic == "Oscillating", ]
  p2 <- NULL

  if (nrow(plot_osc) > 0 && ref_cluster %in% as.character(unique(plot_osc$cluster))) {
    other_clusters <- setdiff(levels(droplevels(plot_osc$cluster)), ref_cluster)

    if (length(other_clusters) > 0) {
      osc_stats <- compare_means(score ~ cluster, data = plot_osc, method = "wilcox.test") %>%
        filter(group1 == ref_cluster | group2 == ref_cluster)

      if (nrow(osc_stats) > 0) {
        osc_stats <- osc_stats %>%
          arrange(group1 == ref_cluster, group2) %>%
          mutate(y.position = seq(0.75, 0.98, length.out = n()))
      }

      osc_fill <- setNames(rep("skyblue", nlevels(droplevels(plot_osc$cluster))),
                            levels(droplevels(plot_osc$cluster)))
      osc_fill[ref_cluster] <- "#FF0187"

      p2 <- ggplot(plot_osc, aes(x = cluster, y = score, fill = cluster)) +
        geom_violin(alpha = 0.9, color = "white") +
        geom_boxplot(width = 0.15, outlier.colour = NA, fill = NA) +
        { if (nrow(osc_stats) > 0) stat_pvalue_manual(osc_stats, label = "p.signif", size = 5) } +
        ylim(0, 1) +
        scale_fill_manual(values = osc_fill) +
        labs(x = "", y = score_label) +
        theme_common(14)
    }
  }

  # --- Plot 3: ref_cluster cells only, score distribution by gene category ---
  plot_ref <- df[df$cluster == ref_cluster, ]
  p3 <- NULL

  if (nrow(plot_ref) > 0) {
    # stat_compare_means()'s automatic bracket placement (geom_signif) can silently drop
    # brackets when one group is much larger than the rest (tick-mark y-positions come back
    # NA). Compute p-values and bracket positions explicitly via stat_pvalue_manual() instead.
    ref_stats <- compare_means(score ~ dynamic, data = plot_ref, method = "wilcox.test") %>%
      filter(group1 == "Oscillating")

    if (nrow(ref_stats) > 0) {
      ref_stats <- ref_stats %>%
        arrange(group2) %>%
        mutate(y.position = seq(0.85, 0.99, length.out = n()))
    }

    p3 <- ggplot(plot_ref, aes(x = dynamic, y = score, fill = dynamic)) +
      geom_violin(alpha = 0.9, color = "white") +
      geom_boxplot(width = 0.15, outlier.colour = NA, fill = NA) +
      { if (nrow(ref_stats) > 0) stat_pvalue_manual(ref_stats, label = "p.signif", size = 5) } +
      ylim(0, 1) +
      scale_fill_manual(values = c(
        "Oscillating"  = "#FF0187",
        "Housekeeping" = "skyblue",
        "Cell Cycle"   = "skyblue",
        "Other"        = "skyblue"
      )) +
      labs(x = "", y = score_label) +
      theme_common(14)
  }

  # --- Plot 4: ref_cluster cells only, gene expression percentile distribution ---
  # Rank genes by cycling score and highlight the segmentation-clock targets,
  # annotating each with the percentile it falls at. Target list is mouse gene
  # symbols, but matching is case-insensitive so it also finds human orthologs
  # (e.g. HES7 in hIPSC data).
  target_genes <- c("Hes7", "Hes1", "Dll1", "Dkk1", "Lfng", "Axin2")
  p4 <- NULL

  if (nrow(plot_ref) > 0) {
    perc_df <- plot_ref %>%
      filter(!is.na(score)) %>%
      arrange(score) %>%
      mutate(percentile = row_number() / n() * 100)

    target_df <- perc_df %>%
      filter(tolower(gene) %in% tolower(target_genes))

    p4 <- ggplot(perc_df, aes(x = percentile, y = score)) +
      geom_line(color = "gray") +
      geom_point(data = target_df, color = "black", size = 3) +
      ggrepel::geom_text_repel(
        data = target_df,
        aes(label = paste0(gene, " (", round(percentile, 2), "%)")),
        color         = "#011c7e",
        box.padding   = 0.5,
        point.padding = 0.5,
        force         = 10,
        segment.color = "gray50",
        direction     = "both",
        max.overlaps  = Inf
      ) +
      labs(x = "Percentile", y = score_label) +
      theme_common(14)
  }

  # --- Save plots (algorithm suffix keeps the two runs from overwriting each other) ---
  if (!dir.exists(config$figures_dir)) dir.create(config$figures_dir, recursive = TRUE)

  if (!is.null(p3)) {
    ggsave(file.path(config$figures_dir, sprintf("%s_category_violin_%s.pdf", ref_cluster, suffix)),
           p3, width = 8, height = 5, dpi = 300)
  }
  if (!is.null(p2)) {
    ggsave(file.path(config$figures_dir, sprintf("oscillating_violin_%s.pdf", suffix)),
           p2, width = 8, height = 5, dpi = 300)
  }
  if (!is.null(p4)) {
    ggsave(file.path(config$figures_dir, sprintf("%s_percentile_%s.pdf", ref_cluster, suffix)),
           p4, width = 12, height = 8, dpi = 300)
  }

  # --- Pathway/TF-activity analysis, run separately per cluster ---
  run_pathway_analysis_per_cluster(df, config, algorithm)

  # --- Score vs mean expression, one plot per cluster -------------------------
  # Sanity check that the cycling score isn't just tracking expression level,
  # with reference oscillating genes highlighted to see if they stand out.
  expr_df       <- compute_mean_expression(counts_files, config$file_prefix)
  expr_figures_dir <- file.path(config$figures_dir, algorithm)

  score_expr_df <- df %>%
    dplyr::inner_join(expr_df, by = c("gene", "cluster")) %>%
    mutate(is_oscillating = dynamic == "Oscillating")

  for (cluster_name in levels(droplevels(score_expr_df$cluster))) {
    plot_score_vs_expression(
      cluster_df   = score_expr_df[score_expr_df$cluster == cluster_name, ],
      cluster_name = cluster_name,
      algorithm    = algorithm,
      score_label  = score_label,
      figures_dir  = expr_figures_dir
    )
  }

  invisible(list(df = df, p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}

# --- Run the evaluation for every algorithm, for every dataset ---
all_plots <- list()

for (dataset_name in names(dataset_configs)) {
  config <- dataset_configs[[dataset_name]]
  refs   <- reference_lists[[config$species]]

  cat("\n\n========================================================================\n")
  cat("DATASET: ", config$label, "\n", sep = "")
  cat("========================================================================\n")

  plots <- lapply(names(algo_display), run_evaluation, config = config, refs = refs)
  names(plots) <- names(algo_display)
  all_plots[[dataset_name]] <- plots

  # Display combined plots for this dataset (cyclum on top row, scPrisma on bottom),
  # only if both algorithms produced at least the category/percentile plots.
  have_cyclum   <- !is.null(plots$cyclum$p2) && !is.null(plots$cyclum$p3)
  have_scPrisma <- !is.null(plots$scPrisma$p2) && !is.null(plots$scPrisma$p3)

  if (have_cyclum && have_scPrisma) {
    print(
      (plots$cyclum$p2 | plots$cyclum$p3) /
      (plots$scPrisma$p2 | plots$scPrisma$p3)
    )
  } else if (have_cyclum) {
    print(plots$cyclum$p2 | plots$cyclum$p3)
  } else if (have_scPrisma) {
    print(plots$scPrisma$p2 | plots$scPrisma$p3)
  }
}
