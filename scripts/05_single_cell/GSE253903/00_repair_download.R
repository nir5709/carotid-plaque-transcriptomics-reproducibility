###############################################################################
# GSE253903 REPAIR DOWNLOAD
# Download all 12 scRNA-seq samples individually from GEO and validate them.
#
# After this script completes successfully, run:
#   source("scripts/05_single_cell/GSE253903/01_stage1_QC_clustering.R")
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

BASE_DIR <- "GSE253903_STAGE1"
RAW_DIR <- file.path(BASE_DIR, "01_raw_10x")
GSM_IDS <- sprintf("GSM%d", 8029962:8029973)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Matrix)
})

if (dir.exists(RAW_DIR)) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup_dir <- paste0(RAW_DIR, "_PARTIAL_BACKUP_", stamp)
  message("Renaming existing raw folder to: ", backup_dir)
  ok <- file.rename(RAW_DIR, backup_dir)
  if (!ok) {
    stop("Could not rename the existing RAW_DIR. Close Explorer/R connections using it and retry.")
  }
}
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

find_triplet <- function(sample_dir) {
  list(
    matrix = list.files(sample_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE,
                        recursive = TRUE, ignore.case = TRUE),
    barcodes = list.files(sample_dir, pattern = "_barcodes\\.tsv\\.gz$", full.names = TRUE,
                          recursive = TRUE, ignore.case = TRUE),
    features = list.files(sample_dir, pattern = "_features.*\\.tsv\\.gz$", full.names = TRUE,
                          recursive = TRUE, ignore.case = TRUE)
  )
}

count_gz_lines <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  n <- 0L
  repeat {
    x <- readLines(con, n = 100000L, warn = FALSE)
    if (length(x) == 0L) break
    n <- n + length(x)
  }
  n
}

validate_gsm <- function(gsm) {
  sample_dir <- file.path(RAW_DIR, gsm)
  trip <- find_triplet(sample_dir)
  if (length(trip$matrix) != 1L || length(trip$barcodes) != 1L || length(trip$features) != 1L) {
    stop(gsm, ": expected one matrix, barcode and feature file. Found ",
         length(trip$matrix), "/", length(trip$barcodes), "/", length(trip$features), ".")
  }
  message("  Reading and validating Matrix Market file...")
  con <- gzfile(trip$matrix, open = "rt")
  mat <- tryCatch(Matrix::readMM(con), finally = close(con))
  n_barcodes <- count_gz_lines(trip$barcodes)
  n_features <- count_gz_lines(trip$features)
  if (nrow(mat) != n_features) stop(gsm, ": matrix rows != feature rows.")
  if (ncol(mat) != n_barcodes) stop(gsm, ": matrix columns != barcode rows.")
  out <- data.frame(
    GSM = gsm,
    Matrix_file = basename(trip$matrix),
    Matrix_MB = round(file.info(trip$matrix)$size / 1024^2, 2),
    Genes = nrow(mat),
    Barcodes = ncol(mat),
    Nonzero_entries = length(mat@x),
    Valid = TRUE,
    stringsAsFactors = FALSE
  )
  rm(mat); gc()
  out
}

download_one <- function(gsm, max_attempts = 3L) {
  for (attempt in seq_len(max_attempts)) {
    message("")
    message("============================================================")
    message(gsm, " — attempt ", attempt, "/", max_attempts)
    message("============================================================")
    sample_dir <- file.path(RAW_DIR, gsm)
    if (dir.exists(sample_dir)) unlink(sample_dir, recursive = TRUE, force = TRUE)
    ok <- tryCatch({
      GEOquery::getGEOSuppFiles(gsm, makeDirectory = TRUE, baseDir = RAW_DIR)
      TRUE
    }, error = function(e) {
      message("Download error: ", conditionMessage(e)); FALSE
    })
    if (!ok) next
    val <- tryCatch(validate_gsm(gsm), error = function(e) {
      message("Validation error: ", conditionMessage(e)); NULL
    })
    if (!is.null(val)) {
      message(gsm, " validated successfully.")
      return(val)
    }
    message(gsm, " failed validation; redownloading.")
  }
  stop(gsm, " failed after ", max_attempts, " attempts.")
}

results <- vector("list", length(GSM_IDS))
for (i in seq_along(GSM_IDS)) results[[i]] <- download_one(GSM_IDS[i])
validation_table <- do.call(rbind, results)

write.csv(validation_table,
          file.path(BASE_DIR, "GSE253903_12sample_download_validation.csv"),
          row.names = FALSE)

matrix_files <- list.files(RAW_DIR, pattern = "_matrix\\.mtx\\.gz$",
                           full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

cat("\n============================================================\n")
cat("GSE253903 DOWNLOAD/VALIDATION COMPLETE\n")
cat("============================================================\n")
cat("Expected samples: 12\n")
cat("Validated samples:", nrow(validation_table), "\n")
cat("Matrix files detected:", length(matrix_files), "\n\n")
print(validation_table)

if (nrow(validation_table) != 12L || length(matrix_files) != 12L) {
  stop("Validation did not finish with all 12 samples.")
}

cat("\nAll 12 samples are valid.\n")
cat("Now run:\n")
cat('source("scripts/05_single_cell/GSE253903/01_stage1_QC_clustering.R")\n')
