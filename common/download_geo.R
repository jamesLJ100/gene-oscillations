# download_geo.R

ensure_geoquery <- function() {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("GEOquery")
  }
  library(GEOquery)
}

#' Download a GEO accession's supplementary files, unless already present
#'
#' @param accession GEO accession to fetch (e.g. "GSM4286760", "GSE162992")
#' @param base_dir Directory passed to GEOquery::getGEOSuppFiles() as baseDir
#' @param required_files Character vector of file paths whose presence (all of
#'   them, via file.exists()) means the download can be skipped
#' @param label Human-readable label for the skip/download message (e.g. "HCT116 raw data")
download_geo_supp_if_missing <- function(accession, base_dir, required_files, label) {
  ensure_geoquery()
  if (!all(file.exists(required_files))) {
    GEOquery::getGEOSuppFiles(accession, baseDir = base_dir)
  } else {
    cat(label, "already present at", required_files[1], "- skipping download.\n")
  }
}
