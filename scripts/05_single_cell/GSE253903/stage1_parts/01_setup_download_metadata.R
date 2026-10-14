###############################################################################
# GSE253903 STAGE 1
# QC, DONOR-AWARE CLUSTERING, AND BROAD CELL-TYPE ANNOTATION PREPARATION
#
# Goal:
#   Build an independent human carotid-plaque scRNA-seq reference from
#   GSE253903 (12 subjects), reproduce the GEO QC thresholds, cluster cells,
#   and generate the marker evidence needed for defensible broad cell-type
#   annotation.
#
# IMPORTANT:
#   This Stage 1 script DOES NOT yet test the brown/blue WGCNA modules.
#   We first lock broad cell-type labels from the actual marker evidence.
#   Stage 2 will then test:
#      Brown inflammatory module -> macrophage / myeloid localization
#      Blue fibrous module        -> SMC / fibroblast localization
#      Locked 65-gene inflammasome set -> cell-type localization
#
# Source:
#   GEO GSE253903
#   12 scRNA-seq carotid plaque subjects (6 symptomatic)
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260923)

###############################################################################
# 0. USER SETTINGS
###############################################################################

GSE <- "GSE253903"
BASE_DIR <- "GSE253903_STAGE1"
DOWNLOAD_DIR <- file.path(BASE_DIR, "00_download")
RAW_DIR <- file.path(BASE_DIR, "01_raw_10x")
OUTDIR <- file.path(BASE_DIR, "02_stage1_outputs")

dir.create(DOWNLOAD_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

AUTO_DOWNLOAD <- TRUE
MIN_FEATURES <- 200L
MAX_FEATURES <- 6000L
MAX_COUNTS <- 40000
MAX_PERCENT_MT <- 30
N_VARIABLE_FEATURES <- 3000L
N_PCS <- 40L
DIMS_USE <- 1:30
CLUSTER_RESOLUTION <- 0.6

###############################################################################
# 1. PACKAGES
###############################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cran_pkgs <- c("Seurat","SeuratObject","Matrix","data.table","ggplot2","patchwork","harmony")
bioc_pkgs <- c("GEOquery")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

missing_pkgs <- c(cran_pkgs, bioc_pkgs)[
  !vapply(c(cran_pkgs, bioc_pkgs), requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(GEOquery)
})

###############################################################################
# 2. HELPERS
###############################################################################

strip_characteristic_prefix <- function(x) {
  x <- as.character(x)
  trimws(sub("^.*?:\\s*", "", x))
}

find_meta_characteristic <- function(pd, keyword) {
  cols <- grep(keyword, colnames(pd), ignore.case = TRUE, value = TRUE)
  if (length(cols) > 0L) return(strip_characteristic_prefix(pd[[cols[1]]]))
  char_cols <- grep("characteristics", colnames(pd), ignore.case = TRUE, value = TRUE)
  if (length(char_cols) == 0L) return(rep(NA_character_, nrow(pd)))
  ans <- rep(NA_character_, nrow(pd))
  for (cc in char_cols) {
    z <- as.character(pd[[cc]])
    hit <- grepl(paste0("^", keyword, "\\s*:"), z, ignore.case = TRUE)
    ans[hit] <- strip_characteristic_prefix(z[hit])
  }
  ans
}

read_one_10x_triplet <- function(matrix_file) {
  prefix <- sub("_matrix\\.mtx\\.gz$", "", matrix_file, ignore.case = TRUE)
  barcode_file <- paste0(prefix, "_barcodes.tsv.gz")
  feature_candidates <- c(
    paste0(prefix, "_features.tsv.gz"),
    paste0(prefix, "_features.tsv.dup_marked.tsv.gz")
  )
  feature_file <- feature_candidates[file.exists(feature_candidates)][1]
  if (!file.exists(barcode_file)) stop("Barcode file missing for: ", matrix_file)
  if (is.na(feature_file) || !file.exists(feature_file)) stop("Feature file missing for: ", matrix_file)
  counts <- Matrix::readMM(gzfile(matrix_file, open = "rt"))
  features <- data.table::fread(feature_file, header = FALSE, data.table = FALSE)
  barcodes <- data.table::fread(barcode_file, header = FALSE, data.table = FALSE)
  if (nrow(features) != nrow(counts)) stop("Feature/count row mismatch for: ", matrix_file)
  if (nrow(barcodes) != ncol(counts)) stop("Barcode/count column mismatch for: ", matrix_file)
  gene_symbols <- if (ncol(features) >= 2L) as.character(features[[2]]) else as.character(features[[1]])
  gene_symbols <- make.unique(gene_symbols)
  rownames(counts) <- gene_symbols
  colnames(counts) <- as.character(barcodes[[1]])
  counts
}

save_plot_both <- function(plot_object, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), plot_object, width = width, height = height)
  ggsave(paste0(stem, ".png"), plot_object, width = width, height = height, dpi = 300)
}

###############################################################################
# 3. DOWNLOAD GEO SUPPLEMENTARY DATA
###############################################################################

raw_tar_expected <- file.path(DOWNLOAD_DIR, GSE, paste0(GSE, "_RAW.tar"))
matrix_files_existing <- list.files(
  RAW_DIR, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE,
  recursive = TRUE, ignore.case = TRUE
)

if (length(matrix_files_existing) == 0L && AUTO_DOWNLOAD) {
  message("Downloading GEO supplementary files for ", GSE, " ...")
  GEOquery::getGEOSuppFiles(GSE, makeDirectory = TRUE, baseDir = DOWNLOAD_DIR)
  tar_candidates <- list.files(
    file.path(DOWNLOAD_DIR, GSE), pattern = "RAW\\.tar$",
    full.names = TRUE, ignore.case = TRUE
  )
  if (length(tar_candidates) == 0L) stop("GEO RAW.tar was not found after download.")
  raw_tar <- tar_candidates[1]
  message("Extracting: ", raw_tar)
  utils::untar(raw_tar, exdir = RAW_DIR)
}

###############################################################################
# 4. DISCOVER 12 RAW 10X SAMPLE MATRICES
###############################################################################

matrix_files <- list.files(
  RAW_DIR, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE,
  recursive = TRUE, ignore.case = TRUE
)
if (length(matrix_files) == 0L) {
  stop("No *_matrix.mtx.gz files found under: ", normalizePath(RAW_DIR, winslash = "/", mustWork = FALSE))
}
message("10x matrices detected: ", length(matrix_files))
if (length(matrix_files) != 12L) {
  warning("Expected 12 GSE253903 scRNA-seq samples, but detected ", length(matrix_files), ". Inspect the extracted files before manuscript use.")
}

gsm_from_file <- basename(matrix_files)
gsm_from_file <- sub("_.*$", "", gsm_from_file)

###############################################################################
# 5. DOWNLOAD GEO SAMPLE METADATA
###############################################################################

message("Loading GEO sample metadata...")
gse_obj <- GEOquery::getGEO(GSE, GSEMatrix = TRUE, getGPL = FALSE)
if (inherits(gse_obj, "ExpressionSet")) eset <- gse_obj else eset <- gse_obj[[1]]
pd <- Biobase::pData(eset)

sample_metadata <- data.frame(
  GSM = rownames(pd),
  Title = if ("title" %in% colnames(pd)) as.character(pd$title) else NA_character_,
  Age = find_meta_characteristic(pd, "age"),
  Sex = find_meta_characteristic(pd, "sex"),
  Symptom = find_meta_characteristic(pd, "symptom presence"),
  stringsAsFactors = FALSE
)

write.csv(sample_metadata, file.path(OUTDIR, "01_GEO_sample_metadata.csv"), row.names = FALSE)
