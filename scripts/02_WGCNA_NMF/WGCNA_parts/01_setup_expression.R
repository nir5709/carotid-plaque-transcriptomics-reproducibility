###############################################################################
# IDEA 2: WHOLE-TRANSCRIPTOME UNSUPERVISED WGCNA
# Dataset: GSE125771 carotid plaque transcriptomes
#
# Discover co-expression modules without using phenotype/signature information;
# only after module discovery are module eigengenes related to the predefined
# transcriptional scores.
###############################################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(timeout = 1200)

SEED <- 20260919
set.seed(SEED)

EXPR_FILE <- "data/frozen/GSE125771_gene_level_matrix_frozen.csv"
EXPR_GENE_COLUMN <- 1
NMF_FILE <- "NMF_gene_modules.csv"
PROJECT_NAME <- "IDEA2_GSE125771_WHOLE_TRANSCRIPTOME_WGCNA"
OUTDIR <- PROJECT_NAME

NETWORK_TYPE <- "signed"
TOM_TYPE <- "signed"
COR_TYPE <- "bicor"
MAX_P_OUTLIERS <- 0.05
MIN_MODULE_SIZE <- 30
DEEP_SPLIT <- 2
MERGE_CUT_HEIGHT <- 0.25
MAX_BLOCK_SIZE <- 5000
POWERS <- c(1:10, seq(12, 30, by = 2))
SFT_R2_TARGET <- 0.80
FALLBACK_SIGNED_POWER_30_40 <- 14
TOP_HUBS_PER_MODULE <- 20
GO_FDR_CUTOFF <- 0.05

###############################################################################
# 1. INSTALL REQUIRED PACKAGES
###############################################################################

cran_packages <- c("WGCNA", "data.table", "ggplot2")
bioc_packages <- c(
  "limma", "GSVA", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
  "impute", "preprocessCore", "GO.db"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

###############################################################################
# 2. LOAD PACKAGES
###############################################################################

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(ggplot2)
  library(limma)
  library(GSVA)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

try(WGCNA::allowWGCNAThreads(), silent = TRUE)

###############################################################################
# 3. OUTPUT FOLDERS
###############################################################################

folders <- c(
  "00_Input_QC", "01_Trait_scores", "02_Soft_threshold",
  "03_WGCNA_modules", "04_Module_trait", "05_Hub_genes", "06_GO_BP",
  "07_NMF_overlap", "08_Reproducibility"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
for (f in folders) {
  dir.create(file.path(OUTDIR, f), recursive = TRUE, showWarnings = FALSE)
}

###############################################################################
# 4. HELPER FUNCTIONS
###############################################################################

clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "-")] <- NA_character_
  toupper(x)
}

write_csv_safe <- function(x, path) {
  if (is.null(x)) return(invisible(NULL))
  if (is.matrix(x)) {
    x <- data.frame(Row = rownames(x), x, check.names = FALSE)
  }
  data.table::fwrite(as.data.frame(x), path, na = "")
  invisible(path)
}

safe_file <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

run_ssgsea <- function(expression_matrix, gene_sets) {
  valid_sets <- gene_sets[lengths(gene_sets) >= 2L]
  parameters <- GSVA::ssgseaParam(
    exprData = as.matrix(expression_matrix),
    geneSets = valid_sets,
    minSize = 2,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE
  )
  GSVA::gsva(parameters, verbose = FALSE)
}

matrix_to_long <- function(value_matrix, value_name) {
  out <- expand.grid(
    Row = rownames(value_matrix),
    Column = colnames(value_matrix),
    stringsAsFactors = FALSE
  )
  out[[value_name]] <- as.vector(value_matrix)
  out
}

###############################################################################
# 5. READ THE LOCKED GENE-LEVEL EXPRESSION MATRIX
###############################################################################

if (!file.exists(EXPR_FILE)) {
  stop(
    paste0(
      "\nExpression matrix was not found:\n", EXPR_FILE,
      "\n\nPlace this script/repository in the same working directory as ",
      EXPR_FILE, " or change EXPR_FILE.\n"
    )
  )
}

expression_table <- data.table::fread(
  EXPR_FILE,
  data.table = FALSE,
  check.names = FALSE
)

gene_column_name <- if (is.numeric(EXPR_GENE_COLUMN)) {
  names(expression_table)[EXPR_GENE_COLUMN]
} else {
  EXPR_GENE_COLUMN
}

if (!gene_column_name %in% names(expression_table)) {
  stop("Gene column not found: ", gene_column_name)
}

genes_raw <- clean_gene(expression_table[[gene_column_name]])
expression_table[[gene_column_name]] <- NULL
expr_full <- as.matrix(expression_table)
storage.mode(expr_full) <- "numeric"
rownames(expr_full) <- genes_raw

###############################################################################
# 6. CLEAN THE MATRIX
###############################################################################

valid_rows <- !is.na(rownames(expr_full)) &
  rownames(expr_full) != "" &
  rowSums(is.na(expr_full)) < ncol(expr_full)

expr_full <- expr_full[valid_rows, , drop = FALSE]
expr_full <- limma::avereps(expr_full, ID = rownames(expr_full))

if (anyNA(expr_full)) {
  row_medians <- apply(expr_full, 1, median, na.rm = TRUE)
  for (i in seq_len(nrow(expr_full))) {
    missing_i <- is.na(expr_full[i, ])
    if (any(missing_i)) expr_full[i, missing_i] <- row_medians[i]
  }
}

gene_sd <- apply(expr_full, 1, sd)
keep_gene <- is.finite(gene_sd) & gene_sd > 0
expr_full <- expr_full[keep_gene, , drop = FALSE]

if (any(!is.finite(expr_full))) {
  stop("Non-finite values remain in the expression matrix.")
}

cat("\nExpression matrix after cleaning:\n")
print(dim(expr_full))

if (ncol(expr_full) != 40) {
  warning(paste("Expected 40 GSE125771 plaques but found", ncol(expr_full)))
}

if (nrow(expr_full) != 31568) {
  warning(paste(
    "The locked workflow previously contained 31,568 genes. This run contains",
    nrow(expr_full), "usable genes after cleaning."
  ))
}

input_summary <- data.frame(
  Input_file = EXPR_FILE,
  Genes = nrow(expr_full),
  Samples = ncol(expr_full),
  Minimum_expression = min(expr_full),
  Maximum_expression = max(expr_full),
  Mean_expression = mean(expr_full)
)

write_csv_safe(
  input_summary,
  file.path(OUTDIR, "00_Input_QC", "Input_summary.csv")
)
