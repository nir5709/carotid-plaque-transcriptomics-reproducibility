###############################################################################
# GSE21545 — EXTERNAL FIXED-SIGNATURE VALIDATION
# Final manuscript-aligned workflow derived from the BioMart-assisted analysis
# used to generate the GSE21545 portion of manuscript Figure 4.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20261001)

OUTDIR <- "GSE21545_FIXED_SIGNATURE_VALIDATION"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
cran_pkgs <- c("ggplot2")
bioc_pkgs <- c("GEOquery","Biobase","GSVA","limma")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(GSVA); library(limma)
})

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL", "-", "---")] <- NA_character_
  x
}

run_ssgsea <- function(expression_matrix, gene_sets) {
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    p <- GSVA::ssgseaParam(
      exprData = as.matrix(expression_matrix), geneSets = gene_sets,
      minSize = 1, maxSize = Inf, alpha = 0.25, normalize = TRUE
    )
    return(GSVA::gsva(p, verbose = FALSE))
  }
  GSVA::gsva(
    as.matrix(expression_matrix), gene_sets, method = "ssgsea",
    kcdf = "Gaussian", min.sz = 1, max.sz = Inf,
    abs.ranking = FALSE, ssgsea.norm = TRUE, verbose = FALSE
  )
}

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
OSTEOGENIC_2 <- c("RUNX2","SPP1")
SIGNATURES_21545 <- list(
  Inflammasome = INFLAMMASOME_65,
  Fibrous = FIBROUS_7,
  Osteogenic = OSTEOGENIC_2
)

cat("Loading GSE21545...\n")
gse <- GEOquery::getGEO("GSE21545", GSEMatrix = TRUE, AnnotGPL = FALSE)
eset <- if (inherits(gse, "ExpressionSet")) gse else gse[[1]]
expr_probe_all <- Biobase::exprs(eset)
pd_all <- Biobase::pData(eset)
storage.mode(expr_probe_all) <- "numeric"

# Reconstruct the same 126 carotid plaque samples from GEO metadata.
plaque_flag <- rep(FALSE, nrow(pd_all))
if ("title" %in% colnames(pd_all)) {
  plaque_flag <- plaque_flag | grepl("plaque", as.character(pd_all$title), ignore.case = TRUE)
}
source_cols <- grep("source", colnames(pd_all), ignore.case = TRUE, value = TRUE)
for (cc in source_cols) {
  plaque_flag <- plaque_flag | grepl("plaque", as.character(pd_all[[cc]]), ignore.case = TRUE)
}
sample_type_cols <- grep("sample.*type|characteristics", colnames(pd_all), ignore.case = TRUE, value = TRUE)
for (cc in sample_type_cols) {
  plaque_flag <- plaque_flag | grepl("carotid plaque", as.character(pd_all[[cc]]), ignore.case = TRUE)
}
if (sum(plaque_flag) != 126L) {
  stop("Expected 126 GSE21545 plaque samples; detected ", sum(plaque_flag), ".")
}

pd <- pd_all[plaque_flag, , drop = FALSE]
expr_probe <- expr_probe_all[, rownames(pd), drop = FALSE]
probe_ids <- rownames(expr_probe)
write.csv(pd, file.path(OUTDIR, "01_GSE21545_126plaque_metadata.csv"), row.names = TRUE)

# Probe-to-gene annotation is supplied by the frozen manuscript mapping
# in the next workflow module.
