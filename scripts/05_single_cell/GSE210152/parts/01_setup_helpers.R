###############################################################################
# GSE210152 — INDEPENDENT MACROPHAGE-STATE VALIDATION
# IDEA 3C: PLIN2hi/TREM1hi inflammatory lipid-associated macrophage analysis
#
# Primary question: does the GSE125771 inflammatory programme localize to an
# inflammatory macrophage/MNP state in independent GSE210152 plaque immune cells?
# GSE210152 is CD45+ immune-cell data and cannot validate the blue structural
# module. Donor-aware pseudobulk and marker enrichment are prioritized.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260923)

GSE <- "GSE210152"
BASE_DIR <- "GSE210152_IDEA3C"
DOWNLOAD_DIR <- file.path(BASE_DIR, "00_download")
OUTDIR <- file.path(BASE_DIR, "01_results")
dir.create(DOWNLOAD_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

H5AD_GZ <- file.path(DOWNLOAD_DIR, "GSE210152_raw.h5ad.gz")
H5AD <- file.path(DOWNLOAD_DIR, "GSE210152_raw.h5ad")
GEO_URL <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/download/?",
  "acc=GSE210152&file=GSE210152_raw.h5ad.gz&format=file"
)
AUTO_DOWNLOAD <- TRUE
MODULE_FILE <- "Whole_transcriptome_module_assignments.csv"
REQUIRE_BROWN_MODULE <- FALSE
APPLY_EXTRA_QC <- FALSE
MIN_FEATURES_EXTRA_QC <- 200L
MAX_PERCENT_MT_EXTRA_QC <- 25
N_VARIABLE_FEATURES <- 3000L
N_PCS <- 40L
DIMS_USE <- 1:30
BROAD_RESOLUTION <- 0.6
MNP_N_VARIABLE_FEATURES <- 2500L
MNP_N_PCS <- 30L
MNP_DIMS_USE <- 1:20
MNP_RESOLUTION <- 0.6
MIN_CELLS_PER_PSEUDOBULK <- 20L

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cran_pkgs <- c(
  "Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2",
  "patchwork", "harmony", "R.utils", "dplyr"
)
bioc_pkgs <- c(
  "zellkonverter", "SingleCellExperiment", "SummarizedExperiment", "edgeR"
)
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
all_pkgs <- c(cran_pkgs, bioc_pkgs)
missing_pkgs <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(data.table)
  library(ggplot2); library(patchwork); library(harmony); library(R.utils)
  library(dplyr); library(zellkonverter); library(SingleCellExperiment)
  library(SummarizedExperiment); library(edgeR)
})

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

state_sets <- list(
  PLIN2_TREM1_iLAM = c("PLIN2","TREM1","TNF","CEBPB","IL1B","OLR1","CCL2","CCL7","CCL20","CXCL2","CXCL3","CXCL8"),
  TREM2_LAM = c("TREM2","CD9","CD36","FABP4","FABP5","NR1H3","APOE","APOC1","APOC2","LIPA","LGMN","PSAP"),
  C1Q = c("C1QA","C1QB","C1QC","F13A1","MRC1","VSIG4","FOLR2"),
  HMOX1 = c("HMOX1","FTL","SLC40A1","NUPR1","SELENOP","PRDX1","GPNMB","CTSB","CTSD"),
  S100A8_IL1B = c("S100A8","S100A9","S100A12","IL1B","NLRP3","TNF","CCL3","CCL4","CCL20","TLR2"),
  S100A8_HSP = c("S100A8","S100A9","DNAJB1","HSPA1A","HSPA1B"),
  IL10_TNFAIP3 = c("IL10","TNFAIP3","GPR183"),
  IFN_response = c("IFI6","ISG15","IFIT1","IFIT3","MX1")
)

key_genes <- unique(c(
  "PLIN2","TREM1","TREM2","IL1B","TNF","CEBPB","CCL2","CCL7","CCL20",
  "CXCL2","CXCL3","CXCL8","NLRP3","PYCARD","CASP1","IL18","GBP5","SPP1",
  "CD36","FABP4","FABP5","APOE","C1QA","C1QB","C1QC","HMOX1","S100A8",
  "S100A9","S100A12","IL10","TNFAIP3"
))

broad_sets <- list(
  Macrophage = c("LST1","FCER1G","TYROBP","CTSD","CTSB","AIF1","CD68","CSF1R","C1QA","C1QB","C1QC","APOE"),
  Monocyte = c("LYZ","FCN1","S100A8","S100A9","S100A12","CTSS","VCAN"),
  Dendritic = c("FCER1A","CD1C","CLEC10A","CLEC9A","CD74"),
  Neutrophil = c("FCGR3B","CSF3R","CXCR2","S100A8","S100A9"),
  T_cell = c("CD3D","CD3E","TRBC1","TRBC2","IL7R","CD8A"),
  NK_cell = c("NKG7","GNLY","KLRD1","PRF1"),
  B_cell = c("MS4A1","CD79A","CD37","CD74"),
  Plasma_cell = c("MZB1","JCHAIN","SDC1"),
  Mast_cell = c("KIT","TPSAB1","TPSB2","CPA3")
)

save_plot_both <- function(plot_object, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), plot_object, width = width, height = height)
  ggsave(paste0(stem, ".png"), plot_object, width = width, height = height, dpi = 300)
}

get_counts <- function(obj) {
  out <- tryCatch(
    Seurat::GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(out)) out <- Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
  out
}

get_data <- function(obj) {
  out <- tryCatch(
    Seurat::GetAssayData(obj, assay = "RNA", layer = "data"),
    error = function(e) NULL
  )
  if (is.null(out)) out <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
  out
}

map_patient_string <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  gsm_map <- c(
    "GSM6422822" = "P1", "GSM6422823" = "P2", "GSM6422824" = "P3",
    "GSM6422825" = "P5", "GSM6422826" = "P6", "GSM6422827" = "P7"
  )
  for (gsm in names(gsm_map)) out[grepl(gsm, x, fixed = TRUE)] <- gsm_map[[gsm]]
  patient_map <- c(
    "Patient1" = "P1", "Patient2" = "P2", "Patient3" = "P3",
    "Patient5" = "P5", "Patient6" = "P6", "Patient7" = "P7"
  )
  for (pat in names(patient_map)) out[grepl(pat, x, ignore.case = TRUE)] <- patient_map[[pat]]
  for (pat in c("P1","P2","P3","P5","P6","P7")) {
    hit <- grepl(paste0("(^|[^A-Za-z0-9])", pat, "([^A-Za-z0-9]|$)"), x, ignore.case = TRUE)
    out[is.na(out) & hit] <- pat
  }
  out
}

# Detect patient identity from metadata first, then cell names.
detect_patient_vector <- function(meta, cell_names) {
  candidate_results <- list()
  for (nm in colnames(meta)) {
    x <- meta[[nm]]
    if (is.atomic(x) || is.factor(x)) {
      mapped <- map_patient_string(x)
      candidate_results[[nm]] <- list(
        mapped = mapped,
        n_mapped = sum(!is.na(mapped)),
        n_unique = length(unique(na.omit(mapped)))
      )
    }
  }
  if (length(candidate_results) > 0L) {
    scores <- vapply(
      candidate_results,
      function(z) z$n_mapped + 100000 * z$n_unique,
      numeric(1)
    )
    best_name <- names(which.max(scores))
    best <- candidate_results[[best_name]]
    if (best$n_unique >= 4L && best$n_mapped >= 0.8 * nrow(meta)) {
      message("Patient identity detected from metadata column: ", best_name)
      return(list(patient = best$mapped, source = best_name))
    }
  }
  mapped_cells <- map_patient_string(cell_names)
  if (length(unique(na.omit(mapped_cells))) >= 4L &&
      sum(!is.na(mapped_cells)) >= 0.8 * length(cell_names)) {
    message("Patient identity detected from cell names.")
    return(list(patient = mapped_cells, source = "cell_names"))
  }
  list(patient = rep(NA_character_, length(cell_names)), source = NA_character_)
}

scale_rows_safe <- function(mat) {
  if (ncol(mat) < 2L) return(mat * 0)
  z <- t(apply(mat, 1, function(x) {
    s <- sd(x)
    if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x)) / s
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

cluster_gene_set_scores <- function(avg_expr, gene_sets) {
  zexpr <- scale_rows_safe(as.matrix(avg_expr))
  do.call(rbind, lapply(names(gene_sets), function(set_name) {
    genes <- intersect(gene_sets[[set_name]], rownames(zexpr))
    if (length(genes) == 0L) return(NULL)
    data.frame(
      Cluster = colnames(zexpr),
      GeneSet = set_name,
      Score = colMeans(zexpr[genes, , drop = FALSE], na.rm = TRUE),
      N_genes_present = length(genes),
      stringsAsFactors = FALSE
    )
  }))
}

find_module_file <- function(preferred) {
  if (file.exists(preferred)) return(preferred)
  candidates <- list.files(
    ".",
    pattern = "(module.*assign|wgcna.*module).*\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  exact <- candidates[basename(candidates) == "Whole_transcriptome_module_assignments.csv"]
  if (length(exact) >= 1L) return(exact[1])
  if (length(candidates) == 1L) return(candidates[1])
  NA_character_
}

read_brown_module <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  # Critical repair: lowercase BEFORE stripping punctuation.
  nms <- gsub("[^a-z0-9]", "", tolower(names(x)))
  gene_candidates <- which(nms %in% c("gene","genes","symbol","genesymbol","geneid"))
  module_candidates <- which(nms %in% c("module","modulecolor","color","modulecolour","moduleassignment"))
  if (length(gene_candidates) == 0L || length(module_candidates) == 0L) {
    stop("Could not identify gene/module columns in: ", path)
  }
  gene_col <- names(x)[gene_candidates[1]]
  module_col <- names(x)[module_candidates[1]]
  genes <- unique(toupper(trimws(as.character(
    x[[gene_col]][tolower(trimws(as.character(x[[module_col]]))) == "brown"]
  ))))
  genes[!is.na(genes) & nzchar(genes)]
}

score_pseudobulk_gene_set <- function(logcpm, genes) {
  genes <- intersect(genes, rownames(logcpm))
  if (length(genes) < 2L) return(rep(NA_real_, ncol(logcpm)))
  sub <- as.matrix(logcpm[genes, , drop = FALSE])
  z <- t(apply(sub, 1, function(x) {
    s <- sd(x)
    if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x)) / s
  }))
  colMeans(z, na.rm = TRUE)
}

fisher_gene_set <- function(markers, geneset, universe) {
  markers <- intersect(unique(markers), universe)
  geneset <- intersect(unique(geneset), universe)
  a <- length(intersect(markers, geneset))
  b <- length(setdiff(markers, geneset))
  c <- length(setdiff(geneset, markers))
  d <- length(universe) - a - b - c
  if (length(markers) == 0L || length(geneset) == 0L || d < 0) {
    return(c(overlap = a, odds_ratio = NA_real_, p_value = NA_real_))
  }
  ft <- fisher.test(matrix(c(a,b,c,d), nrow = 2, byrow = TRUE), alternative = "greater")
  c(overlap = a, odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}
