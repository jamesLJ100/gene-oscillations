# gsea_functions.R
# Functions for running GSEA analyses with Hallmark, Gene Ontology, DoRothEA, and TFEA.ChIP

library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(dorothea)
library(TFEA.ChIP)
library(ExperimentHub)

#' Build and save an NES ranking bar plot from a GSEA results data frame
#'
#' @param df GSEA results data frame with NES/p.adjust columns, sorted by desc(NES) (as
#'   returned by as.data.frame() on a gseaResult, or run_dorothea_gsea()'s return value)
#' @param figures_dir Directory to save the plot
#' @param fname File name prefix for saved figures
#' @param file_suffix Suffix appended after fname for the output filename
#' @param x_label Category axis label (e.g. "Pathway", "Transcription Factor")
#' @param label_col Column in df to use as the bar category (default "ID")
#' @param strip_hallmark_prefix Logical, strip "HALLMARK_" from labels for readability
#' @param top_n Cap on number of bars, keeping the top-N by NES (df is already
#'   desc(NES)-sorted). Defaults to 20 to match every dotplot() call's showCategory = 20,
#'   so the dot/bar plot pair for a given task always covers the same top-N; NULL keeps
#'   every row
#' @return The plot object (invisibly)
plot_nes_barplot <- function(df, figures_dir, fname, file_suffix, x_label = "Pathway",
                              label_col = "ID", strip_hallmark_prefix = FALSE, top_n = 20) {

  plot_df <- if (!is.null(top_n)) head(df, top_n) else df
  plot_df$.label <- plot_df[[label_col]]
  if (strip_hallmark_prefix) plot_df$.label <- gsub("HALLMARK_", "", plot_df$.label)

  p_nes <- ggplot(plot_df, aes(x = reorder(.label, NES), y = NES, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_c(name = "Adj. p-value") +
    coord_flip() +
    labs(
      x = x_label,
      y = "Normalised Enrichment Score"
    ) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      axis.title   = element_text(size = 16),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 14)
    )

  plot_height <- max(6, 2 + nrow(plot_df) * 0.25)
  plot_width  <- max(8, 3 + max(nchar(as.character(plot_df$.label))) * 0.15) + 2

  ggsave(file.path(figures_dir, paste0(fname, "_", file_suffix, ".pdf")),
         plot = p_nes, width = plot_width, height = plot_height, dpi = 300)

  invisible(p_nes)
}

#' Run Hallmark GSEA
#'
#' @param geneList Named numeric vector of gene scores (names = ENTREZ IDs)
#' @param species Character, either "Homo sapiens" or "Mus musculus"
#' @param figures_dir Directory to save output figures
#' @param fname File name prefix for saved figures
#' @param algorithm Algorithm name for plot titles
#' @return GSEA results data frame
run_hallmark_gsea <- function(geneList, species, figures_dir, fname, algorithm) {
  
  # Load Hallmark gene sets
  hallmark <- msigdbr(species = species, category = "H") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    mutate(entrez_gene = as.character(entrez_gene))
  
  cat("Hallmark gene sets loaded:", length(unique(hallmark$gs_name)), "\n")
  
  cat("\n=== Running Hallmark GSEA ===\n")
  gsea_hallmark <- GSEA(geneList,
                        TERM2GENE     = hallmark,
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        verbose       = TRUE,
                        eps           = 1e-10,
                        scoreType     = "pos")
  
  gsea_results_df <- as.data.frame(gsea_hallmark) %>%
    filter(!is.na(ID), !is.na(NES)) %>% arrange(desc(NES))
  cat("Significant hallmark terms found:", nrow(gsea_results_df), "\n")
  
  if (nrow(gsea_results_df) > 0) {
    print(gsea_results_df[, c("ID", "enrichmentScore", "NES", "pvalue", "p.adjust")])
    
    # Dotplot
    p_dot <- dotplot(gsea_hallmark, showCategory = 20)
    ggsave(file.path(figures_dir, paste0(fname, "_hallmark_gsea_dotplot.pdf")),
           plot = p_dot, width = 10, height = 10, dpi = 300)

    # NES bar plot
    plot_nes_barplot(gsea_results_df, figures_dir, fname, "hallmark_NES_barplot",
                      x_label = "Pathway", strip_hallmark_prefix = TRUE)
  } else {
    cat("No significant pathways found for", fname, "\n")
  }
  
  return(gsea_results_df)
}

#' Run Gene Ontology GSEA (Biological Process)
#'
#' @param geneList Named numeric vector of gene scores (names = ENTREZ IDs)
#' @param org_db OrgDb object (e.g., org.Hs.eg.db or org.Mm.eg.db)
#' @param figures_dir Directory to save output figures
#' @param fname File name prefix for saved figures
#' @return GO Biological Process GSEA results data frame
run_go_gsea <- function(geneList, org_db, figures_dir, fname) {

  cat("\n=== Running Gene Ontology GSEA (Biological Process) ===\n")

  gsea_go_bp <- gseGO(geneList,
                      OrgDb         = org_db,
                      ont           = "BP",
                      keyType       = "ENTREZID",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      verbose       = TRUE,
                      eps           = 1e-10,
                      scoreType     = "pos")

  go_bp_df <- as.data.frame(gsea_go_bp) %>%
    filter(!is.na(ID), !is.na(NES)) %>% arrange(desc(NES))
  cat("Significant GO BP terms found:", nrow(go_bp_df), "\n")

  if (nrow(go_bp_df) > 0) {
    print(head(go_bp_df[, c("Description", "NES", "pvalue", "p.adjust")], 10))

    p_dot_go_bp <- dotplot(gsea_go_bp, showCategory = 20)
    ggsave(file.path(figures_dir, paste0(fname, "_go_bp_dotplot.pdf")),
           plot = p_dot_go_bp, width = 12, height = 10, dpi = 300)

    plot_nes_barplot(go_bp_df, figures_dir, fname, "go_bp_NES_barplot",
                      x_label = "Biological Process", label_col = "Description")
  }

  go_bp_df
}

#' Run DoRothEA GSEA
#'
#' @param geneList Named numeric vector of gene scores (names = ENTREZ IDs)
#' @param species Character, either "human" or "mouse"
#' @param org_db OrgDb object (e.g., org.Hs.eg.db or org.Mm.eg.db)
#' @param figures_dir Directory to save output figures
#' @param fname File name prefix for saved figures
#' @return GSEA results data frame
run_dorothea_gsea <- function(geneList, species, org_db, figures_dir, fname) {
  
  # Load appropriate DoRothEA dataset
  if (species == "human") {
    data(dorothea_hs, package = "dorothea")
    dorothea_data <- dorothea_hs
  } else if (species == "mouse") {
    data(dorothea_mm, package = "dorothea")
    dorothea_data <- dorothea_mm
  } else {
    stop("Species must be 'human' or 'mouse'")
  }
  
  target_list <- dorothea_data %>%
    filter(confidence %in% c("A", "B")) %>%
    dplyr::select(tf, target) %>%
    dplyr::rename(Geneset = tf, SYMBOL = target)
  
  dorothea_entrez_map <- bitr(target_list$SYMBOL,
                              fromType = "SYMBOL",
                              toType   = "ENTREZID",
                              OrgDb    = org_db)
  
  dorothea_term2gene <- merge(target_list, dorothea_entrez_map, by = "SYMBOL") %>%
    dplyr::select(Geneset, ENTREZID)
  
  cat("DoRothEA regulons loaded:", length(unique(dorothea_term2gene$Geneset)), "TFs\n")
  
  cat("\n=== Running DoRothEA GSEA ===\n")
  gsea_dorothea <- GSEA(geneList,
                        TERM2GENE     = dorothea_term2gene,
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        verbose       = TRUE,
                        eps           = 1e-10,
                        scoreType     = "pos")
  
  dorothea_df <- as.data.frame(gsea_dorothea) %>%
    filter(!is.na(ID), !is.na(NES)) %>% arrange(desc(NES))
  cat("Significant DoRothEA TFs found:", nrow(dorothea_df), "\n")
  
  if (nrow(dorothea_df) > 0) {
    print(dorothea_df[, c("ID", "enrichmentScore", "NES", "pvalue", "p.adjust")])
    
    # Dotplot
    p_dot_doro <- dotplot(gsea_dorothea, x = "GeneRatio", showCategory = 20) +
      theme(
        axis.text.x  = element_text(size = 20),
        axis.text.y  = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        legend.text  = element_text(size = 20),
        legend.title = element_text(size = 24)
      )
    ggsave(file.path(figures_dir, paste0(fname, "_dorothea_dotplot.pdf")),
           plot = p_dot_doro, width = 10, height = 8, dpi = 300)

    # NES bar plot
    plot_nes_barplot(dorothea_df, figures_dir, fname, "dorothea_NES_barplot",
                      x_label = "Transcription Factor")
  } else {
    cat("No significant DoRothEA TFs found for", fname, "\n")
  }
  
  return(dorothea_df)
}

#' Run TFEA.ChIP Analysis
#'
#' @param geneList Named numeric vector of gene scores (names = ENTREZ IDs)
#' @param figures_dir Directory to save output figures
#' @param fname File name prefix for saved figures
#' @param species Character, either "human" or "mouse"
#' @param tf_filter Optional character vector of TF names to test (NULL = all TFs)
#' @param expressed Logical, whether to restrict to expressed TFs (default TRUE)
#' @return Data frame of TFEA.ChIP results with FDR correction
run_tfea_chip <- function(geneList, figures_dir, fname,
                          species = "human", tf_filter = NULL, expressed = TRUE) {
  
  cat("\n=== Running TFEA.ChIP Analysis ===\n")
  
  if (!species %in% c("human", "mouse")) {
    stop("species must be 'human' or 'mouse'")
  }
  
  gene_scores <- geneList[!is.na(geneList)]
  gene_scores <- sort(gene_scores, decreasing = TRUE)
  
  input_df <- data.frame(
    Genes          = names(gene_scores),
    log2FoldChange = as.numeric(gene_scores),
    pvalue         = 0,
    stringsAsFactors = FALSE
  )
  
  mode <- if (species == "mouse") "m2h" else "h2h"
  
  if (!exists("ChIPDB", envir = .GlobalEnv)) {
    cat("Loading full ChIP-seq database from ExperimentHub (EH9854)...\n")
    eh <- ExperimentHub::ExperimentHub()
    assign("ChIPDB", eh[["EH9854"]], envir = .GlobalEnv)
    cat("Full ChIP-seq database loaded and cached in global environment.\n")
  } else {
    cat("Using cached ChIPDB from global environment.\n")
  }
  
  cat(sprintf("Running TFEA.ChIP analysis (mode = '%s', method = gsea)...\n", mode))
  
  tfea_result <- TFEA.ChIP::analysis_from_table(
    input_df,
    mode      = mode,
    TFfilter  = tf_filter,
    expressed = expressed,
    method    = "gsea"
  )
  
  enrich_tbl <- tfea_result$result$Enrichment.table
  cat("Class of Enrichment.table:", class(enrich_tbl), "\n")
  cat("Dim of Enrichment.table:", nrow(enrich_tbl), "x", ncol(enrich_tbl), "\n")
  cat("Column names:", paste(colnames(enrich_tbl), collapse = ", "), "\n")
  cat("Column classes:\n")
  print(sapply(enrich_tbl, class))
  cat("Head of p.val:\n")
  print(head(enrich_tbl[["p.val"]]))
  cat("----------------------------------\n\n")
  
  tfea_df <- enrich_tbl
  
  # Use pval.adj from TFEA.ChIP directly (already BH-corrected across all datasets)
  tfea_df$FDR         <- tfea_df$pval.adj
  tfea_df$Significant <- tfea_df$FDR < 0.05 & tfea_df$ES > 0.1
  tfea_df             <- tfea_df[order(tfea_df$pval.adj), ]
  
  cat(sprintf("\nTFEA.ChIP Results: %d TFs tested, %d significant (FDR < 0.05 & ES > 0.1)\n",
              nrow(tfea_df),
              sum(tfea_df$Significant, na.rm = TRUE)))
  
  # Save full results to CSV
  csv_path <- file.path(figures_dir, paste0(fname, "_tfea_chip_results.csv"))
  write.csv(tfea_df, file = csv_path, row.names = FALSE)
  cat("Full TFEA.ChIP results saved to:", csv_path, "\n")
  
  if (sum(tfea_df$Significant, na.rm = TRUE) > 0) {
    
    # nperm must match the number of permutations used in analysis_from_table;
    # used as pseudocount floor: smallest resolvable p-value = 1/nperm
    nperm <- 1000
    
    sig_df           <- tfea_df[tfea_df$Significant & !is.na(tfea_df$Significant), ]
    plot_df          <- head(sig_df, 20)
    plot_df$log10FDR <- -log10(plot_df$FDR + 1/nperm)
    plot_df$TF_label <- paste0(plot_df$TF, " (", plot_df$Cell, ")")
    still_dup <- duplicated(plot_df$TF_label) | duplicated(plot_df$TF_label, fromLast = TRUE)
    accession_code <- sub("\\..*$", "", plot_df$Accession[still_dup])
    plot_df$TF_label[still_dup] <- paste0(
      plot_df$TF[still_dup], " (", plot_df$Cell[still_dup], ", ", accession_code, ")"
    )

    print(head(sig_df, 20))

    plot_height <- max(6, 2 + nrow(plot_df) * 0.25)
    plot_width  <- max(8, 3 + max(nchar(plot_df$TF_label)) * 0.15) + 2

    # Bar plot of top significant TFs. Bar height = ES (the effect-size statistic,
    # playing the role NES plays in the other bar plots), fill = -log10(FDR)
    p_bar <- ggplot(plot_df, aes(x = reorder(TF_label, ES), y = ES, fill = log10FDR)) +
      geom_bar(stat = "identity") +
      scale_fill_viridis_c(name = "-log10(FDR)", direction = -1) +
      coord_flip() +
      labs(
        x = "Transcription Factor (Cell)",
        y = "Enrichment Score"
      ) +
      theme_classic() +
      theme(
        axis.text.x  = element_text(size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title   = element_text(size = 16),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14)
      )
    ggsave(file.path(figures_dir, paste0(fname, "_tfea_chip_barplot.pdf")),
           plot = p_bar, width = plot_width, height = plot_height, dpi = 300)

    # Dot plot of top significant TFs - mirrors the dotplot() views the other GSEA
    # functions produce (x = effect-size statistic, colour = significance), using ES
    # (TFEA.ChIP's effect-size statistic, playing the role GeneRatio/NES play elsewhere)
    # since TFEA.ChIP results have no GeneRatio/setSize/Count equivalent to plot.
    p_dot_tfea <- ggplot(plot_df, aes(x = ES, y = reorder(TF_label, ES), colour = log10FDR)) +
      geom_point(size = 3) +
      scale_colour_viridis_c(name = "-log10(FDR)", direction = -1) +
      labs(
        x = "Enrichment Score",
        y = "Transcription Factor (Cell)"
      ) +
      theme_classic() +
      theme(
        axis.text.x  = element_text(size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title   = element_text(size = 16),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14)
      )
    ggsave(file.path(figures_dir, paste0(fname, "_tfea_chip_dotplot.pdf")),
           plot = p_dot_tfea, width = plot_width, height = plot_height, dpi = 300)

    # Volcano plot
    # Add pseudocount of 1/nperm before log transformation to avoid -log10(0) = Inf
    tfea_df$log10FDR <- -log10(tfea_df$FDR + 1/nperm)
    sig_df$log10FDR  <- -log10(sig_df$FDR  + 1/nperm)
    
    # Dynamic y-axis: pad 15% above highest value
    y_upper <- max(tfea_df$log10FDR, na.rm = TRUE) * 1.15
    
    p_scatter <- ggplot(tfea_df, aes(x = ES, y = log10FDR, colour = Significant)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#08007E"),
                          name   = "FDR < 0.05 &\nES > 0.1") +
      ggrepel::geom_text_repel(
        data         = head(sig_df, 20),
        aes(x = ES, y = log10FDR, label = TF),
        size         = 3,
        colour       = "#08007E",
        box.padding  = 0.4,
        max.overlaps = 20,
        ylim         = c(NA, Inf)
      ) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
      scale_y_continuous(
        limits = c(0, y_upper),
        expand = expansion(mult = c(0.02, 0.05))
      ) +
      labs(
        x = "Enrichment Score",
        y = paste0("-log10(FDR + ", 1/nperm, ")")
      ) +
      theme_classic() +
      theme(
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 16),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)
      )
    ggsave(file.path(figures_dir, paste0(fname, "_tfea_chip_volcano.pdf")),
           plot = p_scatter, width = 8, height = 6, dpi = 300)
    
  } else {
    cat("No significant TFs found for", fname, "\n")
  }

  return(tfea_df)
}

#' Build a ranked, Entrez-ID-keyed gene score vector for GSEA
#'
#' Converts gene symbols to Entrez IDs (dropping unmapped genes) and returns a
#' score vector sorted descending.
#'
#' @param results_df Data frame with a gene symbol column and a score column
#' @param org_db OrgDb object (e.g., org.Hs.eg.db or org.Mm.eg.db)
#' @param symbol_col Name of the gene symbol column (default "symbol")
#' @param score_col Name of the score column (default "score")
#' @return Named numeric vector (names = ENTREZ IDs), sorted descending
build_ranked_gene_list <- function(results_df, org_db, symbol_col = "symbol", score_col = "score") {

  genes_entrez <- bitr(results_df[[symbol_col]],
                        fromType = "SYMBOL",
                        toType   = "ENTREZID",
                        OrgDb    = org_db)

  merged <- merge(results_df, genes_entrez, by.x = symbol_col, by.y = "SYMBOL")

  geneList <- merged[[score_col]]
  names(geneList) <- merged$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)

  cat("=== Gene list diagnostics ===\n")
  cat("Genes in geneList:", length(geneList), "\n")
  cat("Score range:", min(geneList), "to", max(geneList), "\n")
  cat("Genes with score = 0:", sum(geneList == 0), "\n")

  geneList
}

#' Resolve a "human"/"mouse" species string to its OrgDb annotation object
#' @param species Character, "human" or "mouse"
#' @return An OrgDb object (org.Hs.eg.db or org.Mm.eg.db)
resolve_org_db <- function(species) {
  if (!species %in% c("human", "mouse")) stop("species must be 'human' or 'mouse'")
  if (species == "human") {
    requireNamespace("org.Hs.eg.db", quietly = TRUE)
    org.Hs.eg.db::org.Hs.eg.db
  } else {
    requireNamespace("org.Mm.eg.db", quietly = TRUE)
    org.Mm.eg.db::org.Mm.eg.db
  }
}

#' Run the full pathway/TF-activity analysis suite (Hallmark, GO, DoRothEA, TFEA.ChIP)
#'
#' @param geneList Named numeric vector of gene scores (names = ENTREZ IDs), as
#'   returned by build_ranked_gene_list()
#' @param species Character, "human" or "mouse"
#' @param figures_dir Directory to save output figures (created if missing)
#' @param fname File name prefix for saved figures
#' @param algorithm Algorithm name for plot titles (e.g. "cyclum", "scPrisma")
#' @param run_go Logical, whether to run the (slower) GO Biological Process GSEA (default TRUE)
#' @param run_tfea Logical, whether to run TFEA.ChIP (default TRUE)
#' @param tf_filter Optional character vector of TF names to test in TFEA.ChIP (NULL = all)
#' @return Named list with elements hallmark, go, dorothea, tfea - each a results data
#'   frame (go/tfea are NULL if skipped)
run_pathway_analyses <- function(geneList, species, figures_dir, fname, algorithm,
                                  run_go = TRUE, run_tfea = TRUE, tf_filter = NULL) {

  org_db <- resolve_org_db(species)
  msigdbr_species <- if (species == "human") "Homo sapiens" else "Mus musculus"

  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

  hallmark_df <- run_hallmark_gsea(geneList, species = msigdbr_species,
                                    figures_dir = figures_dir, fname = fname, algorithm = algorithm)

  go_results <- if (run_go) {
    run_go_gsea(geneList, org_db = org_db, figures_dir = figures_dir, fname = fname)
  } else {
    NULL
  }

  dorothea_df <- run_dorothea_gsea(geneList, species = species, org_db = org_db,
                                    figures_dir = figures_dir, fname = fname)

  tfea_df <- if (run_tfea) {
    run_tfea_chip(geneList, figures_dir = figures_dir, fname = fname,
                  species = species, tf_filter = tf_filter)
  } else {
    NULL
  }

  list(hallmark = hallmark_df, go = go_results, dorothea = dorothea_df, tfea = tfea_df)
}

#' Build and save a Hallmark-GSEA heatmap comparing NES across clusters/conditions
#'
#' @param hallmark_by_group Named list of per-group hallmark GSEA result data
#'   frames (as returned by run_hallmark_gsea(), or the `$hallmark` element of
#'   run_pathway_analyses()'s return value), one per cluster/condition
#' @param figures_dir Directory to save the heatmap
#' @param algorithm Algorithm name, used in the plot title and output filename
#' @return The combined long-format data frame used to build the heatmap
#'   (invisibly), or NULL if no group had any significant pathways
plot_hallmark_heatmap <- function(hallmark_by_group, figures_dir, algorithm) {

  cat("\n=== Creating Hallmark heatmap ===\n")

  hallmark_combined <- bind_rows(
    lapply(names(hallmark_by_group), function(group_name) {
      hr <- hallmark_by_group[[group_name]]
      if (!is.null(hr) && nrow(hr) > 0) {
        hr %>%
          dplyr::select(ID, NES, p.adjust) %>%
          mutate(cluster = group_name)
      }
    })
  )

  if (nrow(hallmark_combined) == 0) {
    cat("No significant hallmark pathways found across any cluster.\n")
    return(invisible(NULL))
  }

  # only keep pathways significant in at least one cluster
  sig_pathways <- hallmark_combined %>%
    filter(p.adjust < 0.05) %>%
    pull(ID) %>%
    unique()

  if (length(sig_pathways) == 0) {
    cat("No significant hallmark pathways found across any cluster.\n")
    return(invisible(NULL))
  }

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
      x = "",
      y = ""
    ) +
    theme_bw(base_size = 18) +
    theme(
      axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 14),
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 16),
      panel.grid   = element_blank()
    )

  ggsave(file.path(figures_dir, paste0(algorithm, "_hallmark_heatmap.pdf")),
         plot   = p_heatmap,
         width  = 4 + length(unique(heatmap_df$cluster)) * 1.2,
         height = 3 + length(sig_pathways) * 0.3,
         dpi    = 300)

  cat("Hallmark heatmap saved.\n")

  invisible(heatmap_df)
}

#' Scatter plot comparing Hallmark GSEA enrichment scores between two conditions
#'
#' One point per Hallmark pathway, NES in condition A vs NES in condition B - answers
#' "did the same pathways shift between these two conditions, and in the same direction?"
#' Points near the y=x diagonal are enriched similarly in both; points far from it are
#' specific to one condition. Since run_hallmark_gsea()'s result table only contains
#' pathways GSEA() found significant *for that condition* (pvalueCutoff = 0.05), a
#' pathway missing from one side's table is plotted at NES = 0 there - i.e. "not
#' significantly enriched in that condition", not "true NES was zero".
#'
#' @param hallmark_a,hallmark_b Hallmark GSEA result data frames (as returned by
#'   run_hallmark_gsea()/run_pathway_analyses()$hallmark), one per condition being
#'   compared. Either may be NULL/empty if that condition had no significant pathways.
#' @param label_a,label_b Condition names, used as axis labels, plot title, and filename
#' @param figures_dir Directory to save the plot
#' @param algorithm Algorithm name, used in the title and output filename
#' @return The merged comparison data frame (invisibly), or NULL if neither condition
#'   had any significant pathways
plot_hallmark_comparison <- function(hallmark_a, hallmark_b, label_a, label_b,
                                      figures_dir, algorithm) {

  cat(sprintf("\n=== Creating Hallmark comparison plot: %s vs %s ===\n", label_a, label_b))

  empty <- data.frame(ID = character(0), NES = numeric(0), p.adjust = numeric(0))
  a <- if (!is.null(hallmark_a) && nrow(hallmark_a) > 0) hallmark_a else empty
  b <- if (!is.null(hallmark_b) && nrow(hallmark_b) > 0) hallmark_b else empty

  if (nrow(a) == 0 && nrow(b) == 0) {
    cat("No significant hallmark pathways for either", label_a, "or", label_b,
        "- skipping comparison plot.\n")
    return(invisible(NULL))
  }

  comparison_df <- dplyr::full_join(
    dplyr::select(a, ID, NES_a = NES, padj_a = p.adjust),
    dplyr::select(b, ID, NES_b = NES, padj_b = p.adjust),
    by = "ID"
  ) %>%
    mutate(
      NES_a       = ifelse(is.na(NES_a), 0, NES_a),
      NES_b       = ifelse(is.na(NES_b), 0, NES_b),
      significant = (!is.na(padj_a) & padj_a < 0.05) | (!is.na(padj_b) & padj_b < 0.05),
      ID          = gsub("HALLMARK_", "", ID)
    )

  axis_lim <- max(abs(c(comparison_df$NES_a, comparison_df$NES_b)), na.rm = TRUE) * 1.15

  p <- ggplot(comparison_df, aes(x = NES_a, y = NES_b)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(aes(colour = significant), size = 2.5, alpha = 0.8) +
    ggrepel::geom_text_repel(
      data         = dplyr::filter(comparison_df, significant),
      aes(label = ID),
      size         = 3,
      max.overlaps = 20
    ) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#08007E"), guide = "none") +
    coord_equal(xlim = c(-axis_lim, axis_lim), ylim = c(-axis_lim, axis_lim)) +
    labs(
      x = paste(label_a, "NES"),
      y = paste(label_b, "NES")
    ) +
    theme_bw(base_size = 14)

  out_file <- file.path(
    figures_dir,
    sprintf("hallmark_comparison_%s_vs_%s.pdf",
            gsub("[^A-Za-z0-9]+", "_", label_a), gsub("[^A-Za-z0-9]+", "_", label_b))
  )
  ggsave(out_file, plot = p, width = 8, height = 7, dpi = 300)
  cat("Hallmark comparison plot saved to:", out_file, "\n")

  invisible(comparison_df)
}