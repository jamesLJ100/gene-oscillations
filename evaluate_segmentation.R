rm(list=ls())
library(here)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "hdfrw.R"))
source(file.path(proj_root, "cyclum_helper.R"))
source(file.path(proj_root, "utils.R"))

# Algorithms to evaluate. Display name is used in plot titles / axis labels,
# and the key is used as a suffix on output filenames.
algo_display <- c(cyclum = "Cyclum", scPrisma = "scPrisma")

# --- Directories ---
tpm_dir      <- file.path(proj_root, "data_old/mme95/tpm")
cyclum_dir   <- file.path(proj_root, "data_old/mme95/tpm/cyclum")
scPrisma_dir <- file.path(proj_root, "data_old/mme95/tpm/scPrisma")


oscillating  <- readRDS("data_old/oscillating_mouse.rds")
housekeeping <- readRDS("data_old/housekeeping_mouse.rds")
cellcycle    <- readRDS("data_old/cellcycle_mouse.rds")

# --- Find all TPM input files ---
tpm_files <- list.files(tpm_dir, pattern = "\\.h5$", full.names = TRUE)
cat("Found", length(tpm_files), "TPM files\n")

if (!dir.exists("figures")) dir.create("figures")

# --- Theme ---
theme_common <- function(base_size = 16) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(-0.25, "cm"),
      legend.position = "none",
      strip.background = element_rect(fill = "transparent", colour = NA)
    )
}

# ============================================================================
# Load per-cluster cycling scores for one algorithm.
# These data_old results use older naming conventions that get_model_scores()
# doesn't build, so read/derive them directly:
#   cyclum   -> <name>_cyclum.h5 weights + the expression matrix
#   scPrisma -> <name>_r0_scPrisma.csv
# ============================================================================
load_scores <- function(algorithm) {
  results <- list()
  for (expr_file in tpm_files) {
    file_name <- tools::file_path_sans_ext(basename(expr_file))
    cluster   <- gsub("mmE95_PM_(.+)_tpm", "\\1", file_name)

    cat("Processing:", cluster, "with", algorithm, "\n")

    if (algorithm == "cyclum") {
      weight_file <- file.path(cyclum_dir, paste0(file_name, "_cyclum.h5"))
      if (!file.exists(weight_file)) {
        cat("  no cyclum weights for:", cluster, "-", weight_file, "\n")
        next
      }
      scores <- get_cyclum_scores(expr_file, weight_file)
    } else if (algorithm == "scPrisma") {
      # run_scPrisma.py writes <name>.csv; older runs used <name>_r0_scPrisma.csv
      results_file <- file.path(scPrisma_dir, paste0(file_name, ".csv"))
      if (!file.exists(results_file)) {
        results_file <- file.path(scPrisma_dir, paste0(file_name, "_r0_scPrisma.csv"))
      }
      if (!file.exists(results_file)) {
        cat("  no scPrisma results for:", cluster, "-", results_file, "\n")
        next
      }
      scores <- read.csv(results_file)
    } else {
      stop("Unknown algorithm: ", algorithm)
    }

    scores$cluster <- cluster
    results[[cluster]] <- scores
  }

  df <- bind_rows(results)
  # get_scores returns "symbol" column — rename to gene
  if ("symbol" %in% colnames(df)) df <- rename(df, gene = symbol)
  df
}

# ============================================================================
# Run the full evaluation (load -> debug -> annotate -> plot -> save) for a
# single algorithm. `score_label` labels the y-axis; `suffix` is appended to
# output filenames so the two algorithms don't overwrite each other.
# ============================================================================
run_evaluation <- function(algorithm) {
  display <- algo_display[[algorithm]]
  suffix  <- algorithm

  cat("\n\n########################################\n")
  cat("### ", display, " evaluation\n", sep = "")
  cat("########################################\n")

  df <- load_scores(algorithm)

  cat("Columns in df:", colnames(df), "\n")
  cat("Total rows:", nrow(df), "\n")

  # ========================================
  # === DEBUGGING SECTION STARTS HERE ===
  # ========================================

  cat("\n=============== DEBUGGING (", display, ") ===============\n", sep = "")

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
  for (i in 2:nrow(cluster_gene_sets)) {
    if (!setequal(cluster_gene_sets$genes[[1]], cluster_gene_sets$genes[[i]])) {
      all_same <- FALSE
      cat("  Cluster", cluster_gene_sets$cluster[[i]], "has different genes!\n")
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

  cat("\n=============== END DEBUGGING (", display, ") ===============\n\n", sep = "")

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
  cluster_order <- c("NMP", "MPC", "pPSM", "aPSM", "scSOM", "dmSOM")
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
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    labs(title = paste0("All genes (", display, ")"), x = "", y = score_label) +
    theme_common(14)

  # --- Plot 2: oscillating genes only, pPSM as reference ---
  plot_osc <- df[df$dynamic == "Oscillating", ]

  p2 <- ggplot(plot_osc, aes(x = cluster, y = score, fill = cluster)) +
    geom_violin(alpha = 0.9, color = "white") +
    geom_boxplot(width = 0.15, outlier.colour = NA, fill = "white") +
    stat_compare_means(
      comparisons = list(
        c("pPSM", "NMP"),
        c("pPSM", "MPC"),
        c("pPSM", "aPSM"),
        c("pPSM", "scSOM"),
        c("pPSM", "dmSOM")
      ),
      label   = "p.signif",
      method  = "wilcox.test",
      label.y = c(0.92, 0.82, 0.80, 0.86, 0.92),
      size    = 5
    ) +
    ylim(0, 1) +
    scale_fill_manual(values = c(
      "NMP"   = "#08007E",
      "MPC"   = "#08007E",
      "pPSM"  = "#FF0187",
      "aPSM"  = "#08007E",
      "scSOM" = "#08007E",
      "dmSOM" = "#08007E"
    )) +
    labs(title = paste0("Oscillating genes only (", display, ")"),
         x = "", y = score_label) +
    theme_common(14)

  # --- Plot 3: pPSM cells only, score distribution by gene category ---
  plot_ppsm <- df[df$cluster == "pPSM", ]

  p3 <- ggplot(plot_ppsm, aes(x = dynamic, y = score, fill = dynamic)) +
    geom_violin(alpha = 0.9, color = "white") +
    geom_boxplot(width = 0.15, outlier.colour = NA, fill = "white") +
    stat_compare_means(
      comparisons = list(
        c("Oscillating", "Housekeeping"),
        c("Oscillating", "Cell Cycle"),
        c("Oscillating", "Other")
      ),
      label  = "p.signif",
      method = "wilcox.test",
      size   = 5
    ) +
    ylim(0, 1) +
    scale_fill_manual(values = c(
      "Oscillating"  = "#FF0187",
      "Housekeeping" = "#08007E",
      "Cell Cycle"   = "#08007E",
      "Other"        = "#08007E"
    )) +
    labs(title = paste0("pPSM cells — score by gene category (", display, ")"),
         x = "", y = score_label) +
    theme_common(14)

  # --- Plot 4: pPSM cells only, gene expression percentile distribution ---
  # Rank genes by cycling score and highlight the segmentation-clock targets,
  # annotating each with the percentile it falls at.
  target_genes <- c("Hes7", "Hes1", "Dll1", "Dkk1", "Lfng", "Axin2")

  perc_df <- plot_ppsm %>%
    filter(!is.na(score)) %>%
    arrange(desc(score)) %>%
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
    labs(title = paste0("pPSM gene percentile distribution (", display, ")"),
         x = "Percentile", y = score_label) +
    theme_common(14)

  # Save plots (algorithm suffix keeps the two runs from overwriting each other)
  ggsave(sprintf("figures/Figure_mmE95_pPSM_category_violin_%s.png", suffix),
         p3, width = 6, height = 5, dpi = 300)

  ggsave(sprintf("figures/Figure_mmE95_oscillating_violin_%s.png", suffix),
         p2, width = 6, height = 5, dpi = 300)

  ggsave(sprintf("figures/Figure_mmE95_pPSM_percentile_%s.png", suffix),
         p4, width = 12, height = 8, dpi = 300)

  invisible(list(df = df, p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}

# --- Run the evaluation for every algorithm ---
plots <- lapply(names(algo_display), run_evaluation)
names(plots) <- names(algo_display)

# Display combined plots (cyclum on top row, scPrisma on bottom)
(plots$cyclum$p2 | plots$cyclum$p3) /
  (plots$scPrisma$p2 | plots$scPrisma$p3)
