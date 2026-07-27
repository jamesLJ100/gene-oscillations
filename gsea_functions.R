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
  
  gsea_results_df <- as.data.frame(gsea_hallmark) %>% arrange(desc(NES))
  cat("Significant hallmark terms found:", nrow(gsea_results_df), "\n")
  
  if (nrow(gsea_results_df) > 0) {
    print(gsea_results_df[, c("ID", "enrichmentScore", "NES", "pvalue", "p.adjust")])
    
    # Dotplot
    p_dot <- dotplot(gsea_hallmark, showCategory = 20,
                     title = paste("Hallmark GSEA -", algorithm, fname))
    ggsave(file.path(figures_dir, paste0(fname, "_hallmark_gsea_dotplot.png")),
           plot = p_dot, width = 10, height = 10, dpi = 300)
    
    # GSEA plot for top pathway
    p_gsea <- gseaplot2(gsea_hallmark, geneSetID = gsea_results_df$ID[1],
                        pvalue_table = TRUE, color = '#08007E', base_size = 20,
                        title = paste("GSEA -", gsea_results_df$ID[1], "-", fname))
    ggsave(file.path(figures_dir, paste0(fname, "_hallmark_gsea_top.png")),
           plot = p_gsea, width = 10, height = 8, dpi = 300)
  } else {
    cat("No significant pathways found for", fname, "\n")
  }
  
  return(gsea_results_df)
}

#' Run Gene Ontology GSEA
#'
#' @param geneList Named numeric vector of gene scores (names = ENTREZ IDs)
#' @param org_db OrgDb object (e.g., org.Hs.eg.db or org.Mm.eg.db)
#' @param figures_dir Directory to save output figures
#' @param fname File name prefix for saved figures
#' @return List of GSEA results for BP, MF, and CC
run_go_gsea <- function(geneList, org_db, figures_dir, fname) {
  
  cat("\n=== Running Gene Ontology GSEA ===\n")
  
  results_list <- list()
  
  # Biological Process
  cat("Running GO Biological Process...\n")
  gsea_go_bp <- gseGO(geneList,
                      OrgDb         = org_db,
                      ont           = "BP",
                      keyType       = "ENTREZID",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      verbose       = TRUE,
                      eps           = 1e-10,
                      scoreType     = "pos")
  
  go_bp_df <- as.data.frame(gsea_go_bp) %>% arrange(desc(NES))
  cat("Significant GO BP terms found:", nrow(go_bp_df), "\n")
  
  if (nrow(go_bp_df) > 0) {
    print(head(go_bp_df[, c("Description", "NES", "pvalue", "p.adjust")], 10))
    
    p_dot_go_bp <- dotplot(gsea_go_bp, showCategory = 20,
                           title = paste("GO Biological Process GSEA -", fname))
    ggsave(file.path(figures_dir, paste0(fname, "_go_bp_dotplot.png")),
           plot = p_dot_go_bp, width = 12, height = 10, dpi = 300)
    
    p_gsea_go_bp <- gseaplot2(gsea_go_bp, geneSetID = go_bp_df$ID[1],
                              pvalue_table = TRUE, color = '#08007E', base_size = 20,
                              title = paste("GO BP -", go_bp_df$Description[1]))
    ggsave(file.path(figures_dir, paste0(fname, "_go_bp_top_gsea.png")),
           plot = p_gsea_go_bp, width = 10, height = 8, dpi = 300)
  }
  
  results_list$BP <- go_bp_df
  
  # Molecular Function
  cat("Running GO Molecular Function...\n")
  gsea_go_mf <- gseGO(geneList,
                      OrgDb         = org_db,
                      ont           = "MF",
                      keyType       = "ENTREZID",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      verbose       = TRUE,
                      eps           = 1e-10,
                      scoreType     = "pos")
  
  go_mf_df <- as.data.frame(gsea_go_mf) %>% arrange(desc(NES))
  cat("Significant GO MF terms found:", nrow(go_mf_df), "\n")
  
  if (nrow(go_mf_df) > 0) {
    print(head(go_mf_df[, c("Description", "NES", "pvalue", "p.adjust")], 10))
    
    p_dot_go_mf <- dotplot(gsea_go_mf, showCategory = 20,
                           title = paste("GO Molecular Function GSEA -", fname))
    ggsave(file.path(figures_dir, paste0(fname, "_go_mf_dotplot.png")),
           plot = p_dot_go_mf, width = 12, height = 10, dpi = 300)
  }
  
  results_list$MF <- go_mf_df
  
  # Cellular Component
  cat("Running GO Cellular Component...\n")
  gsea_go_cc <- gseGO(geneList,
                      OrgDb         = org_db,
                      ont           = "CC",
                      keyType       = "ENTREZID",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      verbose       = TRUE,
                      eps           = 1e-10,
                      scoreType     = "pos")
  
  go_cc_df <- as.data.frame(gsea_go_cc) %>% arrange(desc(NES))
  cat("Significant GO CC terms found:", nrow(go_cc_df), "\n")
  
  if (nrow(go_cc_df) > 0) {
    print(head(go_cc_df[, c("Description", "NES", "pvalue", "p.adjust")], 10))
    
    p_dot_go_cc <- dotplot(gsea_go_cc, showCategory = 20,
                           title = paste("GO Cellular Component GSEA -", fname))
    ggsave(file.path(figures_dir, paste0(fname, "_go_cc_dotplot.png")),
           plot = p_dot_go_cc, width = 12, height = 10, dpi = 300)
  }
  
  results_list$CC <- go_cc_df
  
  return(results_list)
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
  
  dorothea_df <- as.data.frame(gsea_dorothea) %>% arrange(desc(NES))
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
    ggsave(file.path(figures_dir, paste0(fname, "_dorothea_dotplot.png")),
           plot = p_dot_doro, width = 10, height = 8, dpi = 300)
    
    # Top TF GSEA plot
    p_gsea_d <- gseaplot2(gsea_dorothea, geneSetID = dorothea_df$ID[1],
                          pvalue_table = TRUE, color = '#7E0008', base_size = 20,
                          title = paste("DoRothEA GSEA -", dorothea_df$ID[1], "-", fname))
    ggsave(file.path(figures_dir, paste0(fname, "_dorothea_top_gsea.png")),
           plot = p_gsea_d, width = 10, height = 8, dpi = 300)
    
    # NES barplot
    p_nes <- ggplot(dorothea_df, aes(x = reorder(ID, NES), y = NES, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "#08007E", high = "#a8b4e8", name = "Adj. p-value") +
      coord_flip() +
      labs(
        title = paste("DoRothEA TF Activity -", fname),
        x     = "Transcription Factor",
        y     = "Normalised Enrichment Score"
      ) +
      theme_classic() +
      theme(
        axis.text.x  = element_text(size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title   = element_text(size = 16),
        plot.title   = element_text(size = 18, face = "bold"),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14)
      )
    ggsave(file.path(figures_dir, paste0(fname, "_dorothea_NES_barplot.png")),
           plot = p_nes, width = 8, height = 6, dpi = 300)
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
  cat("This may take a few minutes...\n")
  
  tfea_result <- TFEA.ChIP::analysis_from_table(
    input_df,
    mode      = mode,
    TFfilter  = tf_filter,
    expressed = expressed,
    method    = "gsea"
  )
  
  # --- debug extraction ---
  cat("\n--- Debugging result structure ---\n")
  cat("Names at top level:", paste(names(tfea_result), collapse = ", "), "\n")
  cat("Names under $result:", paste(names(tfea_result$result), collapse = ", "), "\n")
  
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
    plot_df          <- head(sig_df, 30)
    plot_df$log10FDR <- -log10(plot_df$FDR + 1/nperm)
    
    print(head(sig_df, 20))
    
    # Bar plot of top significant TFs
    p_bar <- ggplot(plot_df, aes(x = reorder(TF, log10FDR), y = log10FDR, fill = log10FDR)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "#a8b4e8", high = "#08007E", name = "-log10(FDR)") +
      coord_flip() +
      labs(
        title = paste("TFEA.ChIP - Top TFs -", fname),
        x     = "Transcription Factor",
        y     = paste0("-log(p + ", 1/nperm, ")")
      ) +
      theme_classic() +
      theme(
        axis.text.x  = element_text(size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title   = element_text(size = 16),
        plot.title   = element_text(size = 18, face = "bold"),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14)
      )
    ggsave(file.path(figures_dir, paste0(fname, "_tfea_chip_barplot.png")),
           plot = p_bar, width = 8, height = 6, dpi = 300)
    
    # Volcano plot
    # Add pseudocount of 1/nperm before log transformation to avoid -log10(0) = Inf
    tfea_df$log10FDR <- -log10(tfea_df$FDR + 1/nperm)
    sig_df$log10FDR  <- -log10(sig_df$FDR  + 1/nperm)
    
    # Dynamic y-axis: pad 15% above highest value (all finite now due to pseudocount)
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
        ylim         = c(NA, Inf)   # allow labels to extend upward freely
      ) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
      scale_y_continuous(
        limits = c(0, y_upper),
        expand = expansion(mult = c(0.02, 0.05))
      ) +
      labs(
        title = paste("TFEA.ChIP Volcano -", fname),
        x     = "Enrichment Score",
        y     = paste0("-log(p + ", 1/nperm, ")")
      ) +
      theme_classic() +
      theme(
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 16),
        plot.title   = element_text(size = 18, face = "bold",
                                    margin = margin(b = 10)),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)
      )
    ggsave(file.path(figures_dir, paste0(fname, "_tfea_chip_volcano.png")),
           plot = p_scatter, width = 8, height = 6, dpi = 300)
    
  } else {
    cat("No significant TFs found for", fname, "\n")
  }
  
  return(tfea_df)
}