proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/download_geo.R"))

geo_dir <- file.path(proj_root, "segmentation_clock/data/GSE114186")
required_files <- file.path(geo_dir, c(
  "GSE114186_mmE95_CellData.csv.gz",  "GSE114186_mmE95_GeneData.csv.gz",  "GSE114186_mmE95_X.csv.gz",
  "GSE114186_mmESC_CellData.csv.gz",  "GSE114186_mmESC_GeneData.csv.gz",  "GSE114186_mmESC_X.csv.gz",
  "GSE114186_hsIPSC_CellData.csv.gz", "GSE114186_hsIPSC_GeneData.csv.gz", "GSE114186_hsIPSC_X.csv.gz"
))

# GSE114186: mouse E9.5 PM, mESC, and hIPSC segmentation clock data
download_geo_supp_if_missing("GSE114186", file.path(proj_root, "segmentation_clock/data"),
                              required_files, label = "GSE114186 data")
