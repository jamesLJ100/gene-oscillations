if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GEOquery")
#library(GEOquery)
#library(data.table)
#library(biomaRt)

# Create directory structure if not present
dirs <- c(
  "data/GSE116935"
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# download data -- hct116
#getGEOSuppFiles("GSM4286760", destdir = "data/GSM4286760")
#Mouse E9.5 Posterior
#getGEOSuppFiles("GSM3137206", destdir = "data/GSM3137206")
#getGEOSuppFiles("GSE114186", destdir = "data/GSE114186")
#Recapitulating the Human Segmentation Clock with Pluripotent Stem Cells
#mouse oscillating genes?
#gse <- getGEO("GSE116935", destdir = "data/GSE116935")

# mouse_eset <- getGEO(filename = "data/GSE116935/GSE116935-GPL19057_series_matrix.txt.gz")
# 
# # Expression matrix
# exprs_mat <- exprs(mouse_eset)
# 
# # Sample metadata — check time points, conditions
# pheno <- pData(mouse_eset)
# print(pheno[, c("title", "geo_accession", "characteristics_ch1")])

#Realnet
realnet_url <- "https://github.com/dynverse/dyngen/raw/data_files/regulatorycircuits_26_neuron-associated_cells_cancer.rds"
download.file(realnet_url, destfile = here::here("data", "realnet.rds"), mode = "wb")

#Realcount
realcount_url <- "https://github.com/dynverse/dyngen/raw/data_files/zenodo_1443566_real_gold_developing-dendritic-cells_schlitzer.rds"
download.file(realcount_url, destfile = here::here("data", "realcount.rds"), mode = "wb")
