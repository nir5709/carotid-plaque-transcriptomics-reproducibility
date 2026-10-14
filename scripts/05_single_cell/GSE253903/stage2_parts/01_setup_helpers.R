###############################################################################
# GSE253903 — STAGE 2
# INDEPENDENT CELL-TYPE LOCALIZATION OF THE GSE125771 INFLAMMATORY–FIBROUS AXIS
#
# PRIMARY QUESTIONS
# -----------------
# 1. Does the GSE125771 BROWN inflammatory WGCNA module localize preferentially
#    to macrophage/monocyte populations in independent GSE253903 scRNA-seq?
# 2. Does the GSE125771 BLUE fibrous/structural WGCNA module localize
#    preferentially to smooth-muscle cells (SMCs) and fibroblasts?
# 3. Which broad cell types carry the exact locked 65-gene inflammasome score?
#
# IMPORTANT
# ---------
# * Donor-aware validation: individual cells are not treated as biological replicates.
# * Primary scoring uses donor × broad-cell-type pseudobulk counts.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260923)

###############################################################################
# 0. SETTINGS
###############################################################################

RDS_FILE <- file.path(
  "GSE253903_STAGE1",
  "02_stage1_outputs",
  "GSE253903_STAGE1_clustered_unannotated.rds"
)
MODULE_FILE <- "Whole_transcriptome_module_assignments.csv"
ANNOTATION_FILE <- "GSE253903_LOCKED_CLUSTER_ANNOTATION.csv"
OUTDIR <- "GSE253903_STAGE2_CELLTYPE_VALIDATION"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MIN_CELLS_PER_PSEUDOBULK <- 20L
MIN_DONORS_PER_CELLTYPE <- 4L
BROWN_PRIMARY_CELLTYPES <- c("Macrophage", "Monocyte")
BLUE_PRIMARY_CELLTYPES <- c("SMC", "Fibroblast")
EXCLUDED_LABELS <- c("Ambiguous_exclude", "Mixed_exclude")

###############################################################################
# 1. PACKAGES
###############################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cran_pkgs <- c("Seurat","SeuratObject","Matrix","data.table","ggplot2","dplyr","tidyr","pheatmap")
bioc_pkgs <- c("edgeR")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = FALSE)
}
all_pkgs <- c(cran_pkgs, bioc_pkgs)
missing_pkgs <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(edgeR)
})

###############################################################################
# 2. LOCKED 65-GENE INFLAMMASOME SIGNATURE
###############################################################################

inflammasome_65 <- c(
  "ABHD17A","ABHD8","ADCY7","ATAT1","BRCC3","BTK","CARD8","CASR",
  "CCDC92","CD36","CPTP","CSNK1A1","DDX3X","DHX33","EIF2AK2","ELP6",
  "FBXL2","GBP2","GBP5","GKN2","HSPA8","IFI16","IRAK1","IRAK2",
  "IRAK4","IRGM","KCNJ8","KCNK13","KCNK6","LAMP2","LATS1","LATS2",
  "MAPK8","MARCHF5","MARK4","MAVS","MEFV","MYD88","NEK7","NLRC3",
  "NLRP1","NLRP3","NLRP6","OGT","P2RX7","PLCG2","PPP2CA","PRKD1",
  "PTPN22","PYCARD","PYDC1","PYDC2","SIRT2","STMP1","TLR4","TLR6",
  "TLR9","TREM2","TRIM11","TRIM31","TRIM65","USP50","ZDHHC1",
  "ZDHHC12","ZDHHC5"
)

###############################################################################
# 3. LOCKED CLUSTER ANNOTATION
###############################################################################

locked_annotation <- data.frame(
  Cluster = as.character(0:23),
  Broad_label = c(
    "Macrophage","T_cell","T_cell","Fibroblast","Macrophage","SMC",
    "Endothelial","Macrophage","SMC","Monocyte","Dendritic","T_cell",
    "Endothelial","Ambiguous_exclude","Plasma","Mast","SMC","Neutrophil",
    "Pericyte_mural","NK","pDC","B_cell","T_cell","Mixed_exclude"
  ),
  Detailed_label = c(
    "C1Q_FOLR2_macrophage","CD4_IL7R_T","CD8_cytotoxic_T","Fibroblast",
    "APOC1_macrophage","Contractile_SMC","Endothelial","SPP1_TREM1_macrophage",
    "Modulated_SMC","FCN1_S100A8_monocyte","CD1C_DC","T_cell",
    "Arterial_activated_endothelial","Cycling_immune","Plasma_cell","Mast_cell",
    "SOST_osteogenic_like_SMC","Neutrophil","Pericyte_mural","NK_cell","pDC",
    "B_cell","IFN_responsive_T","SMC_macrophage_mixed"
  ),
  stringsAsFactors = FALSE
)

if (file.exists(ANNOTATION_FILE)) {
  tmp_anno <- read.csv(ANNOTATION_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  required_anno <- c("Cluster", "Broad_label")
  if (all(required_anno %in% names(tmp_anno))) {
    tmp_anno$Cluster <- as.character(tmp_anno$Cluster)
    locked_annotation <- tmp_anno
  }
}

###############################################################################
# 4. HELPERS
###############################################################################

find_one_file <- function(preferred, pattern) {
  if (file.exists(preferred)) return(preferred)
  hits <- list.files(".", pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0L) return(NA_character_)
  if (length(hits) > 1L) message("Multiple candidate files found; using:\n", hits[1])
  hits[1]
}

get_counts <- function(obj) {
  x <- tryCatch(
    Seurat::GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(x)) x <- Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
  x
}

read_modules <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  clean_names <- gsub("[^a-z0-9]", "", tolower(names(x)))
  gene_idx <- which(clean_names %in% c("gene","genes","symbol","genesymbol","geneid"))
  module_idx <- which(clean_names %in% c("module","modulecolor","modulecolour","color","moduleassignment"))
  if (length(gene_idx) == 0L || length(module_idx) == 0L) {
    stop("Could not identify gene/module columns in module file.\nColumns: ", paste(names(x), collapse = ", "))
  }
  gene_col <- names(x)[gene_idx[1]]
  module_col <- names(x)[module_idx[1]]
  x$.__gene__ <- toupper(trimws(as.character(x[[gene_col]])))
  x$.__module__ <- tolower(trimws(as.character(x[[module_col]])))
  list(
    brown = unique(x$.__gene__[x$.__module__ == "brown"]),
    blue = unique(x$.__gene__[x$.__module__ == "blue"])
  )
}

row_z_score <- function(mat) {
  mat <- as.matrix(mat)
  z <- t(apply(mat, 1, function(v) {
    s <- sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(v)))
    (v - mean(v, na.rm = TRUE)) / s
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

score_gene_set <- function(logcpm, genes) {
  genes <- intersect(toupper(genes), rownames(logcpm))
  if (length(genes) < 2L) return(rep(NA_real_, ncol(logcpm)))
  sub <- logcpm[genes, , drop = FALSE]
  variable <- apply(sub, 1, var, na.rm = TRUE)
  genes_use <- names(variable)[!is.na(variable) & variable > 0]
  if (length(genes_use) < 2L) return(rep(NA_real_, ncol(logcpm)))
  z <- row_z_score(sub[genes_use, , drop = FALSE])
  colMeans(z, na.rm = TRUE)
}

paired_target_delta <- function(scores, score_col, target_celltypes, comparison_name) {
  patients <- unique(scores$Sample)
  rows <- list(); rr <- 1L
  for (pat in patients) {
    tmp <- scores[scores$Sample == pat, , drop = FALSE]
    target <- tmp[tmp$Broad_label %in% target_celltypes, score_col, drop = TRUE]
    other <- tmp[!tmp$Broad_label %in% target_celltypes, score_col, drop = TRUE]
    target <- target[is.finite(target)]
    other <- other[is.finite(other)]
    if (length(target) == 0L || length(other) == 0L) next
    rows[[rr]] <- data.frame(
      Comparison = comparison_name,
      Score = score_col,
      Sample = pat,
      Target_mean = mean(target),
      Other_mean = mean(other),
      Delta_target_minus_other = mean(target) - mean(other),
      stringsAsFactors = FALSE
    )
    rr <- rr + 1L
  }
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

save_boxplot <- function(df, score_col, title_text, file_stem) {
  med_order <- df |>
    dplyr::group_by(Broad_label) |>
    dplyr::summarise(med = median(.data[[score_col]], na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(med))
  df$Broad_label <- factor(df$Broad_label, levels = med_order$Broad_label)
  p <- ggplot(df, aes(x = Broad_label, y = .data[[score_col]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.65, size = 1.5) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Pseudobulk gene-set score", title = title_text)
  ggsave(paste0(file_stem, ".png"), p, width = 10, height = 6, dpi = 300)
  ggsave(paste0(file_stem, ".pdf"), p, width = 10, height = 6)
}
