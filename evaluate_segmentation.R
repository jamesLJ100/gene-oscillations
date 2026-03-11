library(here)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "hdfrw.R"))
source(file.path(proj_root, "cyclum_helper.R"))

#To evaluate
algorithm <- "scPrisma"

# --- Directories ---
tpm_dir    <- file.path(proj_root, "data/mme95/tpm")
cyclum_dir <- file.path(proj_root, "data/mme95/tpm/cyclum")
scPrisma_dir <- file.path(proj_root, "data/mme95/tpm/scPrisma")


oscillating  <- readRDS("data/oscillating_mouse.rds")
housekeeping <- readRDS("data/housekeeping_mouse.rds")
cellcycle    <- readRDS("data/cellcycle_mouse.rds")

# --- Find all TPM input files ---
tpm_files <- list.files(tpm_dir, pattern = "\\.h5$", full.names = TRUE)
cat("Found", length(tpm_files), "TPM files\n")

# --- Get Cyclum scores for each cluster ---
results <- list()

for (expr_file in tpm_files) {
  file_name   <- tools::file_path_sans_ext(basename(expr_file))
  cluster     <- gsub("mmE95_PM_(.+)_tpm", "\\1", file_name)
  #weight_file <- file.path(cyclum_dir, paste0(file_name, "_cyclum.h5"))
  #weight_file <- file.path(scPrisma_dir, paste0(file_name, "_r0_scPrisma.csv"))
  
  
  cat("Processing:", cluster, "\n")
  scores <- get_model_scores("scPrisma", file_name, expr_file)
  
  if (is.null(scores)) {
    cat("  get_scores returned NULL for:", cluster, "\n")
    next
  }
  
  scores$cluster <- cluster
  results[[cluster]] <- scores
}

# --- Combine and fix column name ---
df <- bind_rows(results)

# get_scores returns "symbol" column — rename to gene
if ("symbol" %in% colnames(df)) df <- rename(df, gene = symbol)

cat("Columns in df:", colnames(df), "\n")
cat("Total rows:", nrow(df), "\n")

# ========================================
# === DEBUGGING SECTION STARTS HERE ===
# ========================================

cat("\n=============== DEBUGGING ===============\n")

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

# 8. Show some example matches and non-matches
cat("\nExample oscillating gene matches:\n")
osc_matches <- unique_genes[unique_genes_lower %in% tolower(oscillating$SYMBOL)]
print(head(osc_matches, 10))

cat("\nExample genes that didn't match any category:\n")
matched_genes <- unique_genes_lower %in% c(
  tolower(oscillating$SYMBOL),
  tolower(housekeeping$SYMBOL),
  tolower(cellcycle$SYMBOL)
)
unmatched <- unique_genes[!matched_genes]
print(head(unmatched, 20))

cat("\n=============== END DEBUGGING ===============\n\n")

# ========================================
# === CONTINUE WITH ORIGINAL SCRIPT ===
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

# --- Plot 1: all genes, all clusters ---
p1 <- ggplot(df, aes(x = cluster, y = score, fill = cluster)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "white") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  labs(title = "All genes", x = "", y = "Cyclum Cycling Score") +
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
  labs(title = "Oscillating genes only", x = "", y = "Cyclum Cycling Score") +
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
  labs(title = "pPSM cells — score by gene category",
       x = "", y = "Cyclum Cycling Score") +
  theme_common(14)

if (!dir.exists("figures")) dir.create("figures")

# Save plots
ggsave("figures/Figure_mmE95_pPSM_category_violin.png",
       p3, width = 6, height = 5, dpi = 300)

ggsave("figures/Figure_mmE95_oscillating_violin.png",
       p2, width = 6, height = 5, dpi = 300)

# Display combined plot
 p2 / p3