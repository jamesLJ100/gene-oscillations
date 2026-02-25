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
register(SerialParam())  # force single threaded - otherwise fgsea hangs...windows things
proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "cyclum_helper.R"))
hct116_dir <- file.path(proj_root, "data/GSM4286760")
cyclum_dir <- file.path(proj_root, "data/GSM4286760/cyclum")
expr_file  <- list.files(hct116_dir, pattern = "\\.h5$", full.names = TRUE)[1]
fname      <- tools::file_path_sans_ext(basename(expr_file))
weight_file <- file.path(cyclum_dir, paste0(fname, "_Cyclum.h5"))

# load hallmark gene sets
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))
cat("Hallmark gene sets loaded:", length(unique(hallmark$gs_name)), "\n")

# load DoRothEA regulons directly from the package
data(dorothea_hs, package = "dorothea")
target_list <- dorothea_hs %>%
  filter(confidence %in% c("A", "B")) %>%
  dplyr::select(tf, target) %>%
  dplyr::rename(Geneset = tf, SYMBOL = target)

dorothea_entrez_map <- bitr(target_list$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dorothea_term2gene <- merge(target_list, dorothea_entrez_map, by = "SYMBOL") %>%
  dplyr::select(Geneset, ENTREZID)
cat("DoRothEA regulons loaded:", length(unique(dorothea_term2gene$Geneset)), "TFs\n")

cat("Processing:", fname, "\n")
cyclum_df <- get_scores(expr_file, weight_file)

# convert gene symbols to entrez IDs
genes_entrez <- bitr(cyclum_df$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
cyclum_df <- merge(cyclum_df, genes_entrez, by.x = "symbol", by.y = "SYMBOL")

# build ranked gene list
geneList <- cyclum_df$score
names(geneList) <- cyclum_df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

cat("=== Gene list diagnostics ===\n")
cat("Genes in geneList:", length(geneList), "\n")
cat("Score range:", min(geneList), "to", max(geneList), "\n")
cat("Genes with score = 0:", sum(geneList == 0), "\n")

if (!dir.exists("figures")) dir.create("figures")
png(file.path("figures", paste0(fname, "_score_distribution.png")))
hist(geneList, breaks = 50, main = paste("Cyclum score distribution -", fname),
     xlab = "Cyclum magnitude score")
dev.off()

# ---- Hallmark GSEA ----------------------------------------------------------
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
                 title = paste("Hallmark GSEA - Cyclum", fname))
ggsave(file.path("figures", paste0(fname, "_cyclum_hallmark_dotplot.png")),
       plot = p_dot, width = 10, height = 10, dpi = 300)

p_gsea <- gseaplot2(gsea_hallmark, geneSetID = results_df$ID[1],
                    pvalue_table = TRUE, color = '#08007E', base_size = 20,
                    title = paste("GSEA -", results_df$ID[1], "-", fname))
ggsave(file.path("figures", paste0(fname, "_cyclum_hallmark_top_gsea.png")),
       plot = p_gsea, width = 10, height = 8, dpi = 300)

# ---- DoRothEA GSEA ----------------------------------------------------------
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
ggsave(file.path("figures", paste0(fname, "_cyclum_dorothea_dotplot.png")),
       plot = p_dot_doro, width = 10, height = 8, dpi = 300)

p_gsea_d <- gseaplot2(gsea_dorothea, geneSetID = dorothea_df$ID[1],
                      pvalue_table = TRUE, color = '#7E0008', base_size = 20,
                      title = paste("DoRothEA GSEA -", dorothea_df$ID[1], "-", fname))
ggsave(file.path("figures", paste0(fname, "_cyclum_dorothea_top_gsea.png")),
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
ggsave(file.path("figures", paste0(fname, "_cyclum_dorothea_NES_barplot.png")),
       plot = p_nes, width = 8, height = 6, dpi = 300)