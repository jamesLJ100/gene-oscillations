proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/download_geo.R"))

raw_file <- file.path(proj_root, "myc/data/GSM4286760/GSM4286760_scRNA_HCT116_WT_10x_umiTable.txt.gz")

download_geo_supp_if_missing("GSM4286760", file.path(proj_root, "myc/data"), raw_file,
                              label = "HCT116 raw data")
