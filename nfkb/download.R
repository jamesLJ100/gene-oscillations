proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/download_geo.R"))

tar_file <- file.path(proj_root, "nfkb/data/GSE162992/GSE162992_RAW.tar")
raw_dir  <- file.path(proj_root, "nfkb/data/GSE162992/GSE162992_RAW")

download_geo_supp_if_missing("GSE162992", file.path(proj_root, "nfkb/data"), tar_file,
                              label = "NF-kB raw archive")

if (!dir.exists(raw_dir) || length(list.files(raw_dir)) == 0) {
  untar(tar_file, exdir = raw_dir)
} else {
  cat("NF-kB raw archive already extracted at", raw_dir, "- skipping untar.\n")
}
