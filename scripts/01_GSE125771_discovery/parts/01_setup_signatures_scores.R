###############################################################################
# GSE125771 DISCOVERY ANALYSIS
# Core manuscript-aligned workflow: input QC, locked scores and continuous
# relationships. Raw GEO data are not stored in the repository.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260722)

EXPR_FILE <- "data/frozen/GSE125771_gene_level_matrix_frozen.csv"
OUTDIR <- "GSE125771_DISCOVERY"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

TOP_N <- 10L
FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 2
NMF_FINAL_RANK <- 2L
NMF_NRUN_FINAL <- 100L
NMF_TOP_GENES_PER_MODULE <- 35L

for (d in c(
  "01_signatures", "02_scores", "03_DEG", "04_GO_BP", "05_MCPcounter",
  "06_NMF", "07_reproducibility"
)) dir.create(file.path(OUTDIR, d), recursive = TRUE, showWarnings = FALSE)

cran_pkgs <- c("data.table", "ggplot2", "NMF", "ppcor")
bioc_pkgs <- c("GSVA", "limma", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(GSVA)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL", "-")] <- NA_character_
  x
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  tst <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  data.frame(rho = unname(tst$estimate), P_value = tst$p.value, N = sum(ok))
}

run_ssgsea <- function(expression_matrix, gene_sets) {
  pars <- GSVA::ssgseaParam(
    exprData = as.matrix(expression_matrix),
    geneSets = gene_sets,
    minSize = 2,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE
  )
  GSVA::gsva(pars, verbose = FALSE)
}

make_score_groups <- function(score_vector, samples, n = TOP_N) {
  names(score_vector) <- samples
  ord <- order(score_vector, decreasing = TRUE)
  group <- rep("Middle", length(score_vector)); names(group) <- samples
  group[names(score_vector)[ord[seq_len(n)]]] <- "High"
  group[names(score_vector)[ord[(length(ord)-n+1L):length(ord)]]] <- "Low"
  data.frame(
    Sample = samples,
    Score = as.numeric(score_vector),
    Rank_high_to_low = rank(-score_vector, ties.method = "first"),
    Group = factor(group[samples], levels = c("Low", "Middle", "High")),
    stringsAsFactors = FALSE
  )
}

if (!file.exists(EXPR_FILE)) stop("Missing expression matrix: ", EXPR_FILE)
raw <- data.table::fread(EXPR_FILE, data.table = FALSE, check.names = FALSE)
genes <- clean_gene(raw[[1]])
raw[[1]] <- NULL
expr_full <- as.matrix(raw)
storage.mode(expr_full) <- "numeric"
rownames(expr_full) <- genes

keep <- !is.na(rownames(expr_full)) & nzchar(rownames(expr_full)) &
  rowSums(is.na(expr_full)) < ncol(expr_full)
expr_full <- expr_full[keep, , drop = FALSE]
expr_full <- limma::avereps(expr_full, ID = rownames(expr_full))
if (anyNA(expr_full)) {
  med <- apply(expr_full, 1, median, na.rm = TRUE)
  for (i in seq_len(nrow(expr_full))) {
    miss <- is.na(expr_full[i, ])
    if (any(miss)) expr_full[i, miss] <- med[i]
  }
}
sdv <- apply(expr_full, 1, sd)
expr_full <- expr_full[is.finite(sdv) & sdv > 0, , drop = FALSE]

if (ncol(expr_full) != 40L) warning("Expected 40 plaques; found ", ncol(expr_full))
if (nrow(expr_full) != 31568L) warning("Expected 31,568 genes after final preprocessing; found ", nrow(expr_full))

# Exact locked signatures used in the final manuscript.
INFLAMMASOME_65 <- c(
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
FIBROUS_7 <- c("ACTA2","CNN1","MYH11","MYOCD","TAGLN","LMOD1","ELN")
RUNX2_SPP1 <- c("RUNX2","SPP1")

signature_source <- list(
  Inflammasome = INFLAMMASOME_65,
  Fibrous = FIBROUS_7,
  Osteogenic = RUNX2_SPP1
)
signature_present <- lapply(signature_source, intersect, y = rownames(expr_full))
signature_missing <- lapply(signature_source, setdiff, y = rownames(expr_full))

coverage <- data.frame(
  Signature = names(signature_source),
  Total_genes = lengths(signature_source),
  Present_genes = lengths(signature_present),
  Missing_genes = lengths(signature_missing),
  Missing_gene_names = vapply(signature_missing, paste, collapse = ";", FUN.VALUE = character(1)),
  stringsAsFactors = FALSE
)
write.csv(coverage, file.path(OUTDIR, "01_signatures", "signature_coverage.csv"), row.names = FALSE)

scores_mat <- run_ssgsea(expr_full, signature_present)
score_df <- data.frame(
  Sample = colnames(expr_full),
  Inflammasome = as.numeric(scores_mat["Inflammasome", ]),
  Fibrous = as.numeric(scores_mat["Fibrous", ]),
  Osteogenic = as.numeric(scores_mat["Osteogenic", ]),
  stringsAsFactors = FALSE
)
write.csv(score_df, file.path(OUTDIR, "02_scores", "samplewise_ssGSEA_scores.csv"), row.names = FALSE)

signatures <- c("Inflammasome", "Fibrous", "Osteogenic")
group_tables <- setNames(lapply(signatures, function(s) {
  x <- make_score_groups(score_df[[s]], score_df$Sample)
  x$Signature <- s
  x
}), signatures)
for (s in signatures) {
  write.csv(group_tables[[s]], file.path(OUTDIR, "02_scores", paste0(s, "_Top10_Bottom10_groups.csv")), row.names = FALSE)
}

pairs <- list(c("Inflammasome","Fibrous"), c("Inflammasome","Osteogenic"), c("Fibrous","Osteogenic"))
score_correlations <- do.call(rbind, lapply(pairs, function(z) {
  ct <- safe_cor(score_df[[z[1]]], score_df[[z[2]]])
  data.frame(Signature_1=z[1], Signature_2=z[2], Spearman_rho=ct$rho,
             P_value=ct$P_value, N=ct$N, stringsAsFactors=FALSE)
}))
score_correlations$FDR_BH_3 <- p.adjust(score_correlations$P_value, method = "BH")
write.csv(score_correlations, file.path(OUTDIR, "02_scores", "signature_score_correlations.csv"), row.names = FALSE)
