library(dyngen)
library(tidyverse)
library(readr)
library(dyno)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)

# Select genes to plot
genes_to_plot <- c("C_TF1", "Target6", "Target250")

## define project root
proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "run_dyngen.R"))
source(file.path(proj_root, "dyngen_utils.R"))


## create directories
fig_dir <- "figures"
data_dir <- "data"
dir.create(fig_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)


set.seed(100)

## shared colours and theme
module_colors <- c(A = "#E75480", B = "#F89821", C = "#3AC6F3", D = "#4D4D4F")
extra_colors <- c(Target = "#02679A", HK = "mediumpurple")

base_theme <- function() {
  theme_bw(base_size = 16) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

model_config <- list(
  num_cells    = 1000L,
  num_tfs      = 4L,
  num_targets  = 100L,
  num_hks      = 100L,
  census_interval = 2,
  tau_seconds     = 30 / 3600,
  weight_bw       = 0.1,
  distance_metric = "euclidean",
  min_tfs_per_module = 1L,
  max_in_degree      = 1L,
  damping            = 0.01,
  target_resampling  = Inf
)

sim_file <- file.path(data_dir, "dyngen", "c1000g504_1_sim.rds")
if (file.exists(sim_file)) {
  sim <- readRDS(sim_file)
  cat(" Simulation loaded from disk (0 seconds)!\n")
} else {
  cat(" Running simulation (5 minutes)...\n")
  sim <- run_simulation(model_config)
  saveRDS(sim, sim_file)
  cat(" Simulation saved! Next runs instant.\n")
}


## OPTIONAL: Generate standard plots
# make_plots(sim$model, sim$expression, sim$sub_expression, fig_dir)

# FIXED MULTI-MODAL TRAJECTORIES PLOT (BUG CORRECTED)
sim_i <- 1
meta_sim <- sim$model$simulations$meta %>%
  filter(simulation_i == sim_i) %>%
  arrange(sim_time)

row_start <- min(which(meta_sim$sim_time == min(meta_sim$sim_time)))
n_rows <- nrow(meta_sim)
counts_sim <- as.matrix(sim$model$simulations$counts[row_start:(row_start + n_rows - 1), ])
times_sim <- meta_sim$sim_time

#  FINAL FIX: Reshape feature_info to have one row per molecule type
feature_info_long <- sim$model$feature_info %>%
  select(feature_id, mol_protein, mol_premrna, mol_mrna) %>%
  pivot_longer(
    cols = c(mol_protein, mol_premrna, mol_mrna),
    names_to = "moltype_col",
    values_to = "mol_feature_id"
  ) %>%
  mutate(
    moltype = case_when(
      moltype_col == "mol_protein" ~ "Protein",
      moltype_col == "mol_premrna" ~ "Pre-mRNA",
      moltype_col == "mol_mrna" ~ "mRNA"
    ),
    gene = feature_id # The base gene name
  ) %>%
  select(mol_feature_id, gene, moltype)

# Now join with the counts data
df <- counts_sim %>%
  as.data.frame() %>%
  tibble::rownames_to_column("time_ix") %>%
  tidyr::pivot_longer(-time_ix, names_to = "mol_feature_id", values_to = "count") %>%
  left_join(feature_info_long, by = "mol_feature_id") %>%
  mutate(
    time_ix = as.numeric(time_ix) - 1,
    time = times_sim[time_ix + 1]
  ) %>%
  filter(!is.na(time), !is.na(gene), !is.na(moltype))

cat("\n Unique genes:\n")
unique_genes <- unique(df$gene)
print(head(unique_genes, 20))



# Check if these genes exist, otherwise pick first 3
if (!all(genes_to_plot %in% unique_genes)) {
  cat("\n Some requested genes not found, using available genes\n")
  genes_to_plot <- head(unique_genes[grepl("Target|HK", unique_genes)], 3)
}

cat("\n Plotting genes:", genes_to_plot, "\n")

df_filtered <- df %>%
  filter(gene %in% genes_to_plot)

cat(" Rows in filtered data:", nrow(df_filtered), "\n")

# Create separate plots for each molecule type
# Define colors for genes
gene_colors <- setNames(c("#E75480", "#F89821", "#3AC6F3"), genes_to_plot)

# mRNA plot
p_mrna <- df_filtered %>%
  filter(moltype == "mRNA") %>%
  ggplot(aes(time, count, color = gene)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(x = "Simulation time", y = "mRNA Count",
       color = "Gene",
       title = "mRNA Expression Trajectories") +
  scale_color_manual(values = gene_colors)

# Pre-mRNA plot
p_premrna <- df_filtered %>%
  filter(moltype == "Pre-mRNA") %>%
  ggplot(aes(time, count, color = gene)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(x = "Simulation time", y = "Pre-mRNA Count",
       color = "Gene",
       title = "Pre-mRNA Expression Trajectories") +
  scale_color_manual(values = gene_colors)

# Protein plot
p_protein <- df_filtered %>%
  filter(moltype == "Protein") %>%
  ggplot(aes(time, count, color = gene)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(x = "Simulation time", y = "Protein Count",
       color = "Gene",
       title = "Protein Expression Trajectories") +
  scale_color_manual(values = gene_colors)

# Display plots
print(p_mrna)
print(p_premrna)
print(p_protein)

# Save individual plots
ggsave("figures/mrna_trajectories.png", p_mrna, width = 10, height = 6, dpi = 300)
ggsave("figures/premrna_trajectories.png", p_premrna, width = 10, height = 6, dpi = 300)
ggsave("figures/protein_trajectories.png", p_protein, width = 10, height = 6, dpi = 300)

# Optional: Create a combined plot with all three panels
library(patchwork)
p_combined <- p_premrna / p_mrna / p_protein
ggsave("figures/multi_modal_gene_trajectories.png", p_combined, width = 10, height = 14, dpi = 300)

cat(" All plots saved to 'figures/'\n")
cat("   - mrna_trajectories.png\n")
cat("   - premrna_trajectories.png\n")
cat("   - protein_trajectories.png\n")
cat("   - multi_modal_gene_trajectories.png (combined)\n")
cat(" Simulation saved to 'data/simulation_results.rds' (~200MB)\n")

## Heatmap of gene expression over simulation time, grouped by module

# Build module assignments
feature_net <- sim[["model"]][["feature_network"]]
feature_info <- sim[["model"]][["feature_info"]]

tf_modules <- feature_net %>%
  filter(!is.na(from_module)) %>%
  select(gene = from, module = from_module) %>%
  distinct()

target_edges <- feature_net %>%
  filter(grepl("^Target", to)) %>%
  select(from, to)

# gene_module <- setNames(tf_modules$module, tf_modules$gene)
# 
# #TODO fix: remove max hops
# for (i in seq_len(10)) {
#   unresolved <- target_edges$to[!(target_edges$to %in% names(gene_module))]
#   if (length(unresolved) == 0) break
#   newly_resolved <- target_edges %>%
#     filter(to %in% unresolved, from %in% names(gene_module)) %>%
#     mutate(module = gene_module[from])
#   if (nrow(newly_resolved) == 0) break
#   gene_module <- c(gene_module, setNames(newly_resolved$module, newly_resolved$to))
# }

#hk_genes <- feature_info$feature_id[grepl("^HK", feature_info$feature_id)]
#gene_module <- c(gene_module, setNames(rep("HK", length(hk_genes)), hk_genes))

feature_net <- sim[["model"]][["feature_network"]]


gene_module <- propagate_module_assignments(tf_modules, target_edges)

# Build wide matrix: genes x time for mRNA
df_mrna <- df %>%
  filter(moltype == "mRNA") %>%
  group_by(gene, time) %>%
  summarise(count = mean(count), .groups = "drop")

mat <- df_mrna %>%
  tidyr::pivot_wider(names_from = time, values_from = count) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# Order genes by module
genes_in_mat <- rownames(mat)
genes_with_module <- genes_in_mat[genes_in_mat %in% names(gene_module)]
mat <- mat[genes_with_module, ]
modules_ordered <- gene_module[genes_with_module]
mat <- mat[order(modules_ordered), ]
modules_ordered <- sort(modules_ordered)

# X-axis labels: find nearest column index for each desired label position
time_vals <- as.numeric(colnames(mat))
label_positions <- c(0, 250, 500, 750, 1000)
col_labels <- rep("", length(time_vals))
for (pos in label_positions) {
  nearest_idx <- which.min(abs(time_vals - pos))
  col_labels[nearest_idx] <- as.character(pos)
}

# Row annotation
ann_colors <- c(module_colors, HK = unname(extra_colors["HK"]))

row_ann <- rowAnnotation(
  Module = modules_ordered,
  col = list(Module = ann_colors),
  show_annotation_name = FALSE
)


# Raw counts heatmap
ht <- Heatmap(
  mat,
  name = "mRNA\nCount",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = col_labels,
  column_names_gp = gpar(fontsize = 9),
  left_annotation = row_ann,
  col = colorRamp2(c(0, max(mat) / 2, max(mat)), c("#02679A", "white", "#E75480")),
  row_split = modules_ordered,
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "Simulation Time",
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  use_raster = TRUE,
  row_gap = unit(0, "mm"),
  column_names_rot = 0
)

png("figures/expression_heatmap_time.png", width = 14, height = 10,
    units = "in", res = 300)
draw(ht)
dev.off()
draw(ht)
cat("  - expression_heatmap_time.png\n")

# Log1p heatmap
mat_log <- log1p(mat)

ht_log <- Heatmap(
  mat_log,
  name = "log1p\nmRNA",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = col_labels,
  column_names_gp = gpar(fontsize = 9),
  left_annotation = row_ann,
  col = colorRamp2(c(0, max(mat_log) / 2, max(mat_log)), c("#02679A", "white", "#E75480")),
  row_split = modules_ordered,
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "Simulation Time",
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  use_raster = TRUE,
  row_gap = unit(0, "mm"),
  column_names_rot = 0
)

png("figures/expression_heatmap_time_log1p.png", width = 14, height = 10,
    units = "in", res = 300)
draw(ht_log)
dev.off()
draw(ht_log)
cat("  - expression_heatmap_time_log1p.png\n")