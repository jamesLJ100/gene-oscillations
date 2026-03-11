library(here)
library(dplyr)
library(ROCR)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(readr)
library(msigdbr)
library(BiocParallel)
library(ggplot2)
library(dorothea)
library(TFEA.ChIP)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

register(SerialParam())  # force single threaded - otherwise fgsea hangs...windows things
proj_root <- here::here()
setwd(proj_root)

source(file.path(proj_root, "utils.R"))

algorithm <- "scPrisma"

hct116_dir <- file.path(proj_root, "data/GSM4286760")
input_data_file  <- list.files(hct116_dir, pattern = "\\.h5$", full.names = TRUE)[1]


fname      <- tools::file_path_sans_ext(basename(input_data_file))


results_df <- get_model_scores(algorithm, hct116_dir, fname)

# convert gene symbols to entrez IDs
genes_entrez <- bitr(results_df$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
results_df <- merge(results_df, genes_entrez, by.x = "symbol", by.y = "SYMBOL")

# build ranked gene list
geneList <- results_df$score
names(geneList) <- results_df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

cat("=== Gene list diagnostics ===\n")
cat("Genes in geneList:", length(geneList), "\n")
cat("Score range:", min(geneList), "to", max(geneList), "\n")
cat("Genes with score = 0:", sum(geneList == 0), "\n")

if (!dir.exists("figures")) dir.create("figures")
png(file.path("figures", paste0(fname, "_score_distribution.png")))
hist(geneList, breaks = 50, main = paste(algorithm, "score distribution -", fname),
     xlab = paste(algorithm, "magnitude score"))
dev.off()

# # ---- Hallmark GSEA ----------------------------------------------------------
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
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

results_df <- as.data.frame(gsea_hallmark) %>% arrange(desc(NES))
cat("Significant hallmark terms found:", nrow(results_df), "\n")
print(results_df[, c("ID", "enrichmentScore", "NES", "pvalue", "p.adjust")])

p_dot <- dotplot(gsea_hallmark, showCategory = 20,
                 title = paste("Hallmark GSEA -", algorithm, fname))
ggsave(file.path("figures", paste0(fname, "_", algorithm, "_hallmark_dotplot.png")),
       plot = p_dot, width = 10, height = 10, dpi = 300)

p_gsea <- gseaplot2(gsea_hallmark, geneSetID = results_df$ID[1],
                    pvalue_table = TRUE, color = '#08007E', base_size = 20,
                    title = paste("GSEA -", results_df$ID[1], "-", fname))
ggsave(file.path("figures", paste0(fname, "_", algorithm, "_hallmark_top_gsea.png")),
       plot = p_gsea, width = 10, height = 8, dpi = 300)

# ---- Gene Ontology GSEA -----------------------------------------------------
cat("\n=== Running Gene Ontology GSEA ===\n")

# Biological Process
cat("Running GO Biological Process...\n")
gsea_go_bp <- gseGO(geneList,
                    OrgDb         = org.Hs.eg.db,
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
  ggsave(file.path("figures", paste0(fname, "_", algorithm, "_go_bp_dotplot.png")),
         plot = p_dot_go_bp, width = 12, height = 10, dpi = 300)
  
  p_gsea_go_bp <- gseaplot2(gsea_go_bp, geneSetID = go_bp_df$ID[1],
                            pvalue_table = TRUE, color = '#08007E', base_size = 20,
                            title = paste("GO BP -", go_bp_df$Description[1]))
  ggsave(file.path("figures", paste0(fname, "_", algorithm, "_go_bp_top_gsea.png")),
         plot = p_gsea_go_bp, width = 10, height = 8, dpi = 300)
}

# Molecular Function
cat("Running GO Molecular Function...\n")
gsea_go_mf <- gseGO(geneList,
                    OrgDb         = org.Hs.eg.db,
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
  ggsave(file.path("figures", paste0(fname, "_", algorithm, "_go_mf_dotplot.png")),
         plot = p_dot_go_mf, width = 12, height = 10, dpi = 300)
}

# Cellular Component
cat("Running GO Cellular Component...\n")
gsea_go_cc <- gseGO(geneList,
                    OrgDb         = org.Hs.eg.db,
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
  ggsave(file.path("figures", paste0(fname, "_", algorithm, "_go_cc_dotplot.png")),
         plot = p_dot_go_cc, width = 12, height = 10, dpi = 300)
}

# ---- DoRothEA GSEA ----------------------------------------------------------
data(dorothea_hs, package = "dorothea")
target_list <- dorothea_hs %>%
  filter(confidence %in% c("A", "B")) %>%
  dplyr::select(tf, target) %>%
  dplyr::rename(Geneset = tf, SYMBOL = target)

dorothea_entrez_map <- bitr(target_list$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
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
print(dorothea_df[, c("ID", "enrichmentScore", "NES", "pvalue", "p.adjust")])

p_dot_doro <- dotplot(gsea_dorothea, x = "GeneRatio", showCategory = 10) +
  theme(
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 24)
  )
ggsave(file.path("figures", paste0(fname, "_", algorithm, "_dorothea_dotplot.png")),
       plot = p_dot_doro, width = 10, height = 8, dpi = 300)

p_gsea_d <- gseaplot2(gsea_dorothea, geneSetID = dorothea_df$ID[1],
                      pvalue_table = TRUE, color = '#7E0008', base_size = 20,
                      title = paste("DoRothEA GSEA -", dorothea_df$ID[1], "-", fname))
ggsave(file.path("figures", paste0(fname, "_", algorithm, "_dorothea_top_gsea.png")),
       plot = p_gsea_d, width = 10, height = 8, dpi = 300)

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
ggsave(file.path("figures", paste0(fname, "_", algorithm, "_dorothea_NES_barplot.png")),
       plot = p_nes, width = 8, height = 6, dpi = 300)

# ---- TFEA.ChIP Analysis -----------------------------------------------------
cat("\n=== Running TFEA.ChIP Analysis ===\n")


library(TFEA.ChIP)



# Prepare input: TFEA.ChIP needs Entrez gene IDs and a ranking metric
# Create a named vector of scores with Entrez IDs as names
gene_scores <- results_df$score
names(gene_scores) <- results_df$ENTREZID

# Sort by score (descending) and remove NAs
gene_scores <- gene_scores[!is.na(gene_scores)]
gene_scores <- sort(gene_scores, decreasing = TRUE)

cat(sprintf("Prepared ranked gene list with %d genes\n", length(gene_scores)))
cat("Score range:", range(gene_scores), "\n")

# Create the input data frame format that TFEA.ChIP expects
# According to vignette: data frame with GeneID and log2FoldChange columns
tfea_input_df <- data.frame(
  Genes = names(gene_scores),
  log2FoldChange = as.numeric(gene_scores),
  stringsAsFactors = FALSE,
  pvalue =0
)

cat("Input data frame structure:\n")
str(tfea_input_df)
cat("First few rows:\n")
print(head(tfea_input_df))

cat("\nRunning TFEA.ChIP analysis using GSEA method...\n")
cat("This may take a few minutes...\n")

# Run TFEA.ChIP analysis
# method options: 'gsea' (Gene Set Enrichment Analysis) or 'mean' (average binding)
tfea_result <- TFEA.ChIP::analysis_from_table(
  tfea_input_df,
  TFfilter = NULL,        # NULL means test all TFs; can specify vector like c('MYC', 'E2F1')
  expressed = TRUE,       # Only consider TFs that are expressed
  method = 'gsea'        # Use GSEA for ranked list
)

cat("TFEA.ChIP analysis completed successfully\n")

# Extract and process results
tfea_df <- as.data.frame(tfea_result$Enrichment.table)

# Add FDR correction
tfea_df$FDR <- p.adjust(tfea_df$p.value, method = "BH")
tfea_df$Significant <- tfea_df$FDR < 0.05

# Sort by p-value
tfea_df <- tfea_df[order(tfea_df$p.value), ]

cat(sprintf("\nTFEA.ChIP Results: %d TFs tested, %d significant (FDR < 0.05)\n",
            nrow(tfea_df),
            sum(tfea_df$Significant, na.rm = TRUE)))