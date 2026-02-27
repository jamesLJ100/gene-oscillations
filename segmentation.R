library(data.table)
library(Matrix)
library(Seurat)
library(here)

proj_root <- here::here()
setwd(proj_root)
data_dir <- file.path(proj_root, "data/GSE114186")

# Read cell and gene metadata first
celldata <- fread(file.path(data_dir, "GSE114186_mmE95_CellData.csv.gz"), header = TRUE) |> 
  as.data.frame()
rownames(celldata) <- make.unique(as.character(celldata[[1]]))
celldata <- celldata[, -1, drop = FALSE]

genedata <- fread(file.path(data_dir, "GSE114186_mmE95_GeneData.csv.gz"), header = TRUE) |> 
  as.data.frame()
rownames(genedata) <- make.unique(as.character(genedata[[1]]))
genedata <- genedata[, -1, drop = FALSE]

# Read X and assign row/col names from metadata
X_raw <- fread(file.path(data_dir, "GSE114186_mmE95_X.csv.gz"), header = FALSE) |> as.matrix()
rownames(X_raw) <- rownames(celldata)
colnames(X_raw) <- rownames(genedata)

dim(X_raw) # should be 4367 x 40523

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts    = Matrix(t(X_raw), sparse = TRUE),
  meta.data = celldata,
  project   = "mmE95"
)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Check what cluster columns the authors provided
colnames(celldata)

# Plot using author's assignments (replace "cluster" with actual column name)