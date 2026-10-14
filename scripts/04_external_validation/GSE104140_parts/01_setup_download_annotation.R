###############################################################################
# GSE104140 — FROZEN MANUSCRIPT-REPRODUCTION INPUTS
#
# Canonical reproduction uses:
#   1) the historical TMM/log2-CPM gene-level expression matrix;
#   2) the historical 32-sample pathology metadata.
#
# These frozen inputs avoid numerical drift caused by changes in
# Entrez-to-HGNC annotation resources across Bioconductor releases.
#
# The frozen expression matrix was originally generated from public
# GSE104140 NCBI GEO RNA-seq counts using the manuscript preprocessing
# workflow. This script does not re-query gene annotations.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260926)

GSE_ID <- "GSE104140"

MODULE_FILE <- "Whole_transcriptome_module_assignments.csv"

FROZEN_EXPR_FILE <- file.path(
  "data",
  "frozen",
  "GSE104140_TMM_logCPM_gene_level_frozen.csv"
)

FROZEN_META_FILE <- file.path(
  "data",
  "frozen",
  "GSE104140_sample_metadata_frozen.csv"
)

OUTDIR <- "GSE104140_FROZEN_CALCIFICATION_VALIDATION_v4"

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################################
# 1. PACKAGES
###############################################################################

cran_pkgs <- c(
  "dplyr",
  "tidyr",
  "ggplot2"
)

missing_pkgs <- cran_pkgs[
  !vapply(
    cran_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs)) {
  install.packages(
    missing_pkgs,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

###############################################################################
# 2. LOCKED SIGNATURES
###############################################################################

INFLAMMASOME_65 <- c(
  "ABHD17A","ABHD8","ADCY7","ATAT1","BRCC3","BTK","CARD8","CASR",
  "CCDC92","CD36","CPTP","CSNK1A1","DDX3X","DHX33","EIF2AK2","ELP6",
  "FBXL2","GBP2","GBP5","GKN2","HSPA8","IFI16","IRAK1","IRAK2",
  "IRAK4","IRGM","KCNJ8","KCNK13","KCNK6","LAMP2","LATS1","LATS2",
  "MAPK8","MARCHF5","MARK4","MAVS","MYD88","NEK7","NLRC3",
  "NLRP1","NLRP3","NLRP6","OGT","P2RX7","PLCG2","PPP2CA","PRKD1",
  "PTPN22","PYCARD","PYDC1","PYDC2","SIRT2","STMP1","TLR4","TLR6",
  "TLR9","TREM2","TRIM11","TRIM31","TRIM65","USP50","ZDHHC1",
  "ZDHHC12","ZDHHC5","MEFV"
)

FIBROUS_CHECK <- length(unique(INFLAMMASOME_65))

if (FIBROUS_CHECK != 65L) {
  stop(
    "INFLAMMASOME_65 must contain exactly 65 unique genes."
  )
}

RUNX2_SPP1 <- c(
  "RUNX2",
  "SPP1"
)

###############################################################################
# 3. LOAD FROZEN GSE125771 WGCNA MODULE ASSIGNMENTS
###############################################################################

if (!file.exists(MODULE_FILE)) {

  hits <- list.files(
    ".",
    pattern = "^Whole_transcriptome_module_assignments\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(hits) == 1L) {

    MODULE_FILE <- hits[1]

  } else {

    stop(
      "Cannot find a unique Whole_transcriptome_module_assignments.csv. ",
      "Generate it with the WGCNA workflow or set MODULE_FILE explicitly."
    )
  }
}

modules <- read.csv(
  MODULE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("Gene", "Module") %in% names(modules))) {
  stop(
    "MODULE_FILE must contain columns named Gene and Module."
  )
}

modules$Gene <- toupper(
  trimws(
    as.character(modules$Gene)
  )
)

modules$Module <- tolower(
  trimws(
    as.character(modules$Module)
  )
)

BROWN <- unique(
  modules$Gene[
    modules$Module == "brown"
  ]
)

BLUE <- unique(
  modules$Gene[
    modules$Module == "blue"
  ]
)

if (length(BROWN) != 4262L) {
  stop(
    "Expected 4262 Brown genes; detected ",
    length(BROWN),
    "."
  )
}

if (length(BLUE) != 5197L) {
  stop(
    "Expected 5197 Blue genes; detected ",
    length(BLUE),
    "."
  )
}

GENE_SETS <- list(
  Brown_WGCNA = BROWN,
  Blue_WGCNA = BLUE,
  Inflammasome_65 = INFLAMMASOME_65,
  RUNX2_SPP1 = RUNX2_SPP1
)

###############################################################################
# 4. LOAD FROZEN GSE104140 EXPRESSION MATRIX
###############################################################################

if (!file.exists(FROZEN_EXPR_FILE)) {
  stop(
    "Missing frozen GSE104140 expression matrix: ",
    FROZEN_EXPR_FILE
  )
}

expr_df <- read.csv(
  FROZEN_EXPR_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  nrow(expr_df) != 20861L ||
  ncol(expr_df) != 33L
) {
  stop(
    "Unexpected frozen GSE104140 matrix dimensions: ",
    nrow(expr_df),
    " x ",
    ncol(expr_df),
    "; expected 20861 x 33."
  )
}

if (!identical(names(expr_df)[1], "Gene")) {
  stop(
    "First column of frozen GSE104140 matrix must be Gene."
  )
}

genes <- toupper(
  trimws(
    as.character(expr_df$Gene)
  )
)

if (
  any(is.na(genes)) ||
  any(!nzchar(genes))
) {
  stop(
    "Frozen GSE104140 matrix contains missing/blank gene names."
  )
}

if (anyDuplicated(genes)) {
  stop(
    "Frozen GSE104140 gene-level matrix contains duplicated gene names."
  )
}

expr_gene <- as.matrix(
  expr_df[
    ,
    -1,
    drop = FALSE
  ]
)

storage.mode(expr_gene) <- "numeric"
rownames(expr_gene) <- genes

if (!all(is.finite(expr_gene))) {
  stop(
    "Frozen GSE104140 expression matrix contains non-finite values."
  )
}

cat(
  "\nFrozen GSE104140 expression matrix loaded: ",
  nrow(expr_gene),
  " genes x ",
  ncol(expr_gene),
  " samples\n",
  sep = ""
)

###############################################################################
# 5. LOAD FROZEN GSE104140 SAMPLE METADATA
###############################################################################

if (!file.exists(FROZEN_META_FILE)) {
  stop(
    "Missing frozen GSE104140 metadata: ",
    FROZEN_META_FILE
  )
}

meta <- read.csv(
  FROZEN_META_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_meta <- c(
  "GSM",
  "Title",
  "Pathology_group",
  "Metadata_text"
)

missing_meta <- setdiff(
  required_meta,
  names(meta)
)

if (length(missing_meta)) {
  stop(
    "Frozen metadata is missing columns: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

if (nrow(meta) != 32L) {
  stop(
    "Expected 32 frozen GSE104140 metadata rows; found ",
    nrow(meta),
    "."
  )
}

meta$GSM <- trimws(
  as.character(meta$GSM)
)

if (anyDuplicated(meta$GSM)) {
  stop(
    "Duplicated GSM accessions in frozen GSE104140 metadata."
  )
}

expected_groups <- c(
  Calcified_FA = 11L,
  Noncalcified_FA = 8L,
  DIT = 8L,
  Other_exclude = 5L
)

observed_groups <- table(
  factor(
    meta$Pathology_group,
    levels = names(expected_groups)
  )
)

if (
  !all(
    as.integer(observed_groups) ==
      expected_groups
  )
) {
  stop(
    "Frozen GSE104140 pathology-group counts do not match ",
    "11 calcified, 8 non-calcified, 8 DIT and 5 excluded."
  )
}

###############################################################################
# 6. VERIFY SAMPLE ALIGNMENT
###############################################################################

if (
  !setequal(
    colnames(expr_gene),
    meta$GSM
  )
) {
  stop(
    "Frozen expression matrix and frozen metadata contain ",
    "different GSM sample sets."
  )
}

meta <- meta[
  match(
    colnames(expr_gene),
    meta$GSM
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(expr_gene),
    meta$GSM
  )
)

cat(
  "\nFrozen pathology groups:\n"
)

print(
  table(meta$Pathology_group)
)

###############################################################################
# 7. RECREATE STANDARD OUTPUT INPUT FILES
###############################################################################

write.csv(
  meta,
  file.path(
    OUTDIR,
    "01_GSE104140_sample_metadata.csv"
  ),
  row.names = FALSE
)

write.csv(
  meta[
    meta$Pathology_group %in%
      c(
        "Calcified_FA",
        "Noncalcified_FA"
      ),
    ,
    drop = FALSE
  ],
  file.path(
    OUTDIR,
    "01B_PRIMARY_19_advanced_fibroatheroma_samples.csv"
  ),
  row.names = FALSE
)

write.csv(
  meta[
    meta$Pathology_group == "Other_exclude",
    ,
    drop = FALSE
  ],
  file.path(
    OUTDIR,
    "01C_EXCLUDED_other_pathology_samples.csv"
  ),
  row.names = FALSE
)

copied <- file.copy(
  FROZEN_EXPR_FILE,
  file.path(
    OUTDIR,
    "02_GSE104140_TMM_logCPM_gene_level.csv"
  ),
  overwrite = TRUE
)

if (!copied) {
  stop(
    "Could not copy frozen GSE104140 expression matrix to output directory."
  )
}

cat(
  "\nFrozen GSE104140 manuscript-reproduction inputs ready.\n"
)
