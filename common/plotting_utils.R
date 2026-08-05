library(ggplot2)
library(dplyr)

#' Shared minimal ggplot theme for score-comparison violin/boxplots
#'
#' @param base_size Base font size (default 16)
#' @param angle_x_labels Logical, rotate x-axis labels 45 degrees (useful when
#'   cluster/condition names are long, e.g. nfkb's "wt_tnf"/"ss_pic")
#' @return A ggplot2 theme object
theme_common <- function(base_size = 16, angle_x_labels = FALSE) {
  t <- theme_bw(base_size = base_size) +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.ticks.length = unit(-0.25, "cm"),
      legend.position   = "none",
      strip.background  = element_rect(fill = "transparent", colour = NA)
    )
  if (angle_x_labels) {
    t <- t + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  t
}

#' Scatter plot comparing per-gene oscillation scores between two conditions
#'
#' One point per gene, cycling score in condition A vs condition B - the gene-level
#' analogue of plot_hallmark_comparison() (which compares pathway-level enrichment
#' scores instead). Answers "which specific genes does the algorithm call cyclic in
#' one condition but not the other?" Genes present in only one condition's score
#' table are dropped (rather than 0-filled), since 0 would misleadingly imply
#' "present and scored as non-cyclic" when it's actually "not scored at all here."
#'
#' @param scores_a,scores_b Data frames with `symbol` and `score` columns, one per
#'   condition being compared
#' @param label_a,label_b Condition names, used as axis labels, title, and filename
#' @param figures_dir Directory to save the plot
#' @param algorithm Algorithm name, used in the title/axis labels and output filename
#' @param top_n Number of most score-divergent genes to label (default 20)
#' @return The merged comparison data frame (invisibly), or NULL if there are no
#'   genes shared between the two conditions' score tables
plot_score_comparison <- function(scores_a, scores_b, label_a, label_b,
                                   figures_dir, algorithm, top_n = 20) {

  cat(sprintf("\n=== Creating oscillation score comparison plot: %s vs %s ===\n", label_a, label_b))

  if (is.null(scores_a) || is.null(scores_b) || nrow(scores_a) == 0 || nrow(scores_b) == 0) {
    cat("Missing scores for", label_a, "or", label_b, "- skipping comparison plot.\n")
    return(invisible(NULL))
  }

  comparison_df <- dplyr::inner_join(
    dplyr::select(scores_a, symbol, score_a = score),
    dplyr::select(scores_b, symbol, score_b = score),
    by = "symbol"
  )

  if (nrow(comparison_df) == 0) {
    cat("No genes shared between", label_a, "and", label_b, "- skipping comparison plot.\n")
    return(invisible(NULL))
  }

  comparison_df$score_diff <- abs(comparison_df$score_a - comparison_df$score_b)
  top_genes <- comparison_df %>% dplyr::arrange(dplyr::desc(score_diff)) %>% head(top_n)

  p <- ggplot(comparison_df, aes(x = score_a, y = score_b)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(alpha = 0.35, size = 1, colour = "grey50") +
    geom_point(data = top_genes, colour = "#08007E", size = 1.8) +
    ggrepel::geom_text_repel(
      data         = top_genes,
      aes(label = symbol),
      size         = 3,
      max.overlaps = 20
    ) +
    labs(
      x = paste(label_a, algorithm, "score"),
      y = paste(label_b, algorithm, "score")
    ) +
    theme_common(14)

  out_file <- file.path(
    figures_dir,
    sprintf("score_comparison_%s_vs_%s.pdf",
            gsub("[^A-Za-z0-9]+", "_", label_a), gsub("[^A-Za-z0-9]+", "_", label_b))
  )
  ggsave(out_file, plot = p, width = 8, height = 7, dpi = 300)
  cat("Oscillation score comparison plot saved to:", out_file, "\n")

  invisible(comparison_df)
}

#' Scatter plot of oscillation score vs mean expression, one point per gene
#'
#' Standard sanity check: a cycling score should reflect periodicity, not just
#' expression magnitude, so a strong score-vs-expression trend here is a red flag.
#' Genes from a reference oscillating-gene list are highlighted (and drawn on top
#' of the background cloud) so you can see whether they stand out on score
#' regardless of where they sit on the expression axis.
#'
#' @param cluster_df Data frame with `score`, `mean_expr`, `is_oscillating`
#'   columns, already subset to one cell type/cluster
#' @param cluster_name Cell type/cluster name, used in the title and filename
#' @param algorithm Algorithm name, used in the title and filename
#' @param score_label Y-axis label (e.g. "Cyclum Cycling Score")
#' @param figures_dir Directory to save the plot
#' @return The plot (invisibly), or NULL if cluster_df is empty
plot_score_vs_expression <- function(cluster_df, cluster_name, algorithm, score_label, figures_dir) {

  if (is.null(cluster_df) || nrow(cluster_df) == 0) {
    cat("No data for", cluster_name, "-", algorithm, "- skipping score-vs-expression plot.\n")
    return(invisible(NULL))
  }

  # draw oscillating-gene points last (on top of the background cloud)
  cluster_df <- cluster_df[order(cluster_df$is_oscillating), ]

  p <- ggplot(cluster_df, aes(x = log1p(mean_expr), y = score, colour = is_oscillating)) +
    geom_point(aes(size = is_oscillating, alpha = is_oscillating)) +
    scale_colour_manual(
      values = c(`FALSE` = "grey70", `TRUE` = "#FF0187"),
      labels = c(`FALSE` = "Other", `TRUE` = "Oscillating (reference list)"),
      name   = ""
    ) +
    scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 1.8), guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.5, `TRUE` = 0.9), guide = "none") +
    labs(
      x = "Mean expression (log1p raw counts)",
      y = score_label
    ) +
    theme_common(14) +
    theme(legend.position = "right")

  out_file <- file.path(
    figures_dir,
    sprintf("%s_score_vs_expression_%s.pdf", cluster_name, algorithm)
  )
  ggsave(out_file, plot = p, width = 7, height = 6, dpi = 300)
  cat("Score-vs-expression plot saved to:", out_file, "\n")

  invisible(p)
}
