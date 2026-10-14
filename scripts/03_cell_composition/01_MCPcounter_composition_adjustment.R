###############################################################################
# IDEA 3A
# DOES THE INFLAMMASOME SIGNAL PERSIST AFTER ACCOUNTING FOR MONOCYTIC ABUNDANCE?
#
# Dataset: GSE125771
#
# INPUTS
#   1) data/frozen/GSE125771_gene_level_matrix_frozen.csv
#   2) IDEA2_GSE125771_WHOLE_TRANSCRIPTOME_WGCNA/
#      08_Reproducibility/IDEA2_WGCNA_core_objects.RData
#
# OUTPUTS
#   IDEA3A_COMPOSITION_ADJUSTMENT/
#
# IMPORTANT
# This is a composition-adjustment analysis.
# It does NOT prove macrophage-intrinsic activation.
# Cell-specific localization is evaluated separately in the single-cell datasets.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

EXPR_FILE <- "data/frozen/GSE125771_gene_level_matrix_frozen.csv"

WGCNA_RDATA <- file.path(
  "IDEA2_GSE125771_WHOLE_TRANSCRIPTOME_WGCNA",
  "08_Reproducibility",
  "IDEA2_WGCNA_core_objects.RData"
)

OUTDIR <- "IDEA3A_COMPOSITION_ADJUSTMENT"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1. PACKAGES
###############################################################################

cran_pkgs <- c(
  "data.table",
  "ggplot2",
  "ppcor"
)

bioc_pkgs <- c(
  "limma",
  "GSVA"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

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

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ppcor)
  library(limma)
  library(GSVA)
})

###############################################################################
# 2. HELPERS
###############################################################################

clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "-")] <- NA_character_
  toupper(x)
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) {
    return(data.frame(rho = NA_real_, P_value = NA_real_, N = sum(ok)))
  }
  z <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  data.frame(rho = unname(z$estimate), P_value = z$p.value, N = sum(ok))
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

###############################################################################
# 3. READ EXPRESSION MATRIX
###############################################################################

if (!file.exists(EXPR_FILE)) stop("Missing expression matrix: ", EXPR_FILE)

tab <- data.table::fread(EXPR_FILE, data.table = FALSE, check.names = FALSE)
gene_col <- names(tab)[1]
genes <- clean_gene(tab[[gene_col]])
tab[[gene_col]] <- NULL
expr <- as.matrix(tab)
storage.mode(expr) <- "numeric"
rownames(expr) <- genes

keep <- !is.na(rownames(expr)) & rownames(expr) != "" & rowSums(is.na(expr)) < ncol(expr)
expr <- expr[keep, , drop = FALSE]
expr <- limma::avereps(expr, ID = rownames(expr))

if (anyNA(expr)) {
  med <- apply(expr, 1, median, na.rm = TRUE)
  for (i in seq_len(nrow(expr))) {
    miss <- is.na(expr[i, ])
    if (any(miss)) expr[i, miss] <- med[i]
  }
}

sdv <- apply(expr, 1, sd)
expr <- expr[is.finite(sdv) & sdv > 0, , drop = FALSE]
cat("Expression:", nrow(expr), "genes x", ncol(expr), "samples\n")

###############################################################################
# 4. RECREATE THE LOCKED THREE SIGNATURES
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

fibrous_7 <- c("ACTA2","CNN1","MYH11","MYOCD","TAGLN","LMOD1","ELN")
osteogenic_2 <- c("RUNX2","SPP1")

gene_sets <- list(
  Inflammasome = intersect(inflammasome_65, rownames(expr)),
  Fibrous = intersect(fibrous_7, rownames(expr)),
  Osteogenic = intersect(osteogenic_2, rownames(expr))
)

scores <- run_ssgsea(expr, gene_sets)
scores_df <- data.frame(Sample = colnames(expr), t(scores), check.names = FALSE)

###############################################################################
# 5. MCP-COUNTER
###############################################################################

if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  stop(
    "MCPcounter is required for the canonical composition analysis. ",
    "Run scripts/00_install_dependencies.R and retry."
  )
}

mcp_matrix <- MCPcounter::MCPcounter.estimate(
  expr,
  featuresType = "HUGO_symbols"
)

cat("\nMCP-counter cell populations:\n")
print(rownames(mcp_matrix))

###############################################################################
# 6. IDENTIFY MONOCYTIC LINEAGE
###############################################################################

monocyte_name <- grep("monocytic|monocyte|macrophage", rownames(mcp_matrix), ignore.case = TRUE, value = TRUE)[1]
if (is.na(monocyte_name)) stop("Could not identify monocytic/macrophage MCP-counter row.")
cat("\nUsing abundance proxy:", monocyte_name, "\n")

common_samples <- Reduce(intersect, list(scores_df$Sample, colnames(mcp_matrix), colnames(expr)))

if (file.exists(WGCNA_RDATA)) {
  load(WGCNA_RDATA)
  common_samples <- intersect(common_samples, rownames(MEs))
} else {
  warning(paste0("WGCNA RData not found: ", WGCNA_RDATA, "\nBrown/Blue analyses will be skipped."))
}

###############################################################################
# 7. BUILD SAMPLE-LEVEL ANALYSIS TABLE
###############################################################################

analysis_df <- data.frame(Sample = common_samples, stringsAsFactors = FALSE)
analysis_df$Inflammasome <- scores_df$Inflammasome[match(common_samples, scores_df$Sample)]
analysis_df$Fibrous <- scores_df$Fibrous[match(common_samples, scores_df$Sample)]
analysis_df$Osteogenic <- scores_df$Osteogenic[match(common_samples, scores_df$Sample)]
analysis_df$Monocytic_lineage <- as.numeric(mcp_matrix[monocyte_name, common_samples])
analysis_df$Mean_MCPcounter <- colMeans(mcp_matrix[, common_samples, drop = FALSE], na.rm = TRUE)

if ("PTPRC" %in% rownames(expr)) {
  analysis_df$PTPRC <- as.numeric(expr["PTPRC", common_samples])
} else {
  analysis_df$PTPRC <- NA_real_
}

if (exists("MEs")) {
  if ("MEbrown" %in% colnames(MEs)) analysis_df$Brown_ME <- MEs[common_samples, "MEbrown"]
  if ("MEblue" %in% colnames(MEs)) analysis_df$Blue_ME <- MEs[common_samples, "MEblue"]
}

###############################################################################
# 8. PRIMARY CORRELATIONS
###############################################################################

rows <- list()
add_result <- function(label, x, y) {
  z <- safe_cor(x, y)
  data.frame(Test = label, rho = z$rho, P_value = z$P_value, N = z$N)
}

rows[[1]] <- add_result("Inflammasome_vs_Monocytic_lineage", analysis_df$Inflammasome, analysis_df$Monocytic_lineage)
if ("Brown_ME" %in% names(analysis_df)) rows[[2]] <- add_result("Brown_ME_vs_Monocytic_lineage", analysis_df$Brown_ME, analysis_df$Monocytic_lineage)
if ("Blue_ME" %in% names(analysis_df)) rows[[3]] <- add_result("Blue_ME_vs_Monocytic_lineage", analysis_df$Blue_ME, analysis_df$Monocytic_lineage)

cor_results <- do.call(rbind, rows)
cor_results$FDR <- p.adjust(cor_results$P_value, method = "BH")

###############################################################################
# 9. PARTIAL SPEARMAN
###############################################################################

partial_rows <- list()
if (all(is.finite(analysis_df$PTPRC))) {
  p1 <- ppcor::pcor.test(
    analysis_df$Inflammasome,
    analysis_df$Monocytic_lineage,
    analysis_df$PTPRC,
    method = "spearman"
  )
  partial_rows[[1]] <- data.frame(
    Test = "Inflammasome_vs_Monocytic_adjusted_for_PTPRC",
    Partial_rho = unname(p1$estimate), P_value = p1$p.value, N = p1$n
  )
  p2 <- ppcor::pcor.test(
    analysis_df$Inflammasome,
    analysis_df$Monocytic_lineage,
    cbind(analysis_df$PTPRC, analysis_df$Mean_MCPcounter),
    method = "spearman"
  )
  partial_rows[[2]] <- data.frame(
    Test = "Inflammasome_vs_Monocytic_adjusted_for_PTPRC_and_mean_MCPcounter",
    Partial_rho = unname(p2$estimate), P_value = p2$p.value, N = p2$n
  )
}

partial_results <- if (length(partial_rows) > 0) do.call(rbind, partial_rows) else data.frame()
if (nrow(partial_results) > 0) partial_results$FDR <- p.adjust(partial_results$P_value, method = "BH")

###############################################################################
# 10. RESIDUALISE INFLAMMASOME SCORE FOR COMPOSITION
###############################################################################

if (all(is.finite(analysis_df$PTPRC))) {
  abundance_model <- lm(Inflammasome ~ Monocytic_lineage + PTPRC + Mean_MCPcounter, data = analysis_df)
} else {
  abundance_model <- lm(Inflammasome ~ Monocytic_lineage + Mean_MCPcounter, data = analysis_df)
}

analysis_df$Inflammasome_residual <- residuals(abundance_model)
analysis_df$Inflammasome_fitted_by_composition <- fitted(abundance_model)

###############################################################################
# 11. DOES THE BROWN PROGRAM PERSIST AFTER COMPOSITION ADJUSTMENT?
###############################################################################

resid_results <- list()
if ("Brown_ME" %in% names(analysis_df)) {
  resid_results[[1]] <- add_result("Brown_ME_vs_composition_adjusted_Inflammasome", analysis_df$Brown_ME, analysis_df$Inflammasome_residual)
}
if ("Blue_ME" %in% names(analysis_df)) {
  resid_results[[2]] <- add_result("Blue_ME_vs_composition_adjusted_Inflammasome", analysis_df$Blue_ME, analysis_df$Inflammasome_residual)
}
resid_results <- if (length(resid_results) > 0) do.call(rbind, resid_results) else data.frame()
if (nrow(resid_results) > 0) resid_results$FDR <- p.adjust(resid_results$P_value, method = "BH")

###############################################################################
# 12. CORE INFLAMMATORY GENES AFTER COMPOSITION ADJUSTMENT
###############################################################################

core_genes <- c("NLRP3","PYCARD","CASP1","IL1B","IL18","GBP5","TREM1","MMP9","MMP12","SPP1")
core_genes <- intersect(core_genes, rownames(expr))
gene_rows <- list()

for (g in core_genes) {
  y <- as.numeric(expr[g, common_samples])
  if (all(is.finite(analysis_df$PTPRC))) {
    model_g <- lm(y ~ analysis_df$Monocytic_lineage + analysis_df$PTPRC + analysis_df$Mean_MCPcounter)
  } else {
    model_g <- lm(y ~ analysis_df$Monocytic_lineage + analysis_df$Mean_MCPcounter)
  }
  gene_resid <- residuals(model_g)
  ct <- safe_cor(gene_resid, analysis_df$Inflammasome_residual)
  gene_rows[[g]] <- data.frame(Gene = g, Residual_vs_residual_rho = ct$rho, P_value = ct$P_value, N = ct$N)
}

core_gene_results <- do.call(rbind, gene_rows)
rownames(core_gene_results) <- NULL
core_gene_results$FDR <- p.adjust(core_gene_results$P_value, method = "BH")

###############################################################################
# 13. VARIANCE EXPLAINED BY COMPOSITION MODEL
###############################################################################

model_summary <- summary(abundance_model)
variance_summary <- data.frame(
  N = nrow(analysis_df),
  R_squared = model_summary$r.squared,
  Adjusted_R_squared = model_summary$adj.r.squared,
  Model_P = pf(model_summary$fstatistic[1], model_summary$fstatistic[2], model_summary$fstatistic[3], lower.tail = FALSE)
)

###############################################################################
# 14. EXPORT
###############################################################################

data.table::fwrite(analysis_df, file.path(OUTDIR, "Sample_level_composition_adjusted_data.csv"))
data.table::fwrite(cor_results, file.path(OUTDIR, "Primary_MCPcounter_correlations.csv"))
if (nrow(partial_results) > 0) data.table::fwrite(partial_results, file.path(OUTDIR, "Partial_correlations.csv"))
if (nrow(resid_results) > 0) data.table::fwrite(resid_results, file.path(OUTDIR, "Residualized_module_results.csv"))
data.table::fwrite(core_gene_results, file.path(OUTDIR, "Core_gene_residual_results.csv"))
data.table::fwrite(variance_summary, file.path(OUTDIR, "Inflammasome_variance_explained_by_composition.csv"))

###############################################################################
# 15. FIGURE
###############################################################################

p <- ggplot(analysis_df, aes(x = Monocytic_lineage, y = Inflammasome)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic(base_size = 12) +
  labs(
    title = "Inflammasome activity vs monocytic-lineage abundance proxy",
    x = "MCP-counter monocytic-lineage score",
    y = "Inflammasome ssGSEA score"
  )

ggsave(file.path(OUTDIR, "Inflammasome_vs_Monocytic_lineage.pdf"), p, width = 6, height = 5)
ggsave(file.path(OUTDIR, "Inflammasome_vs_Monocytic_lineage.png"), p, width = 6, height = 5, dpi = 300)

###############################################################################
# 16. FINAL SUMMARY
###############################################################################

sink(file.path(OUTDIR, "FINAL_IDEA3A_SUMMARY.txt"))
cat("IDEA 3A: COMPOSITION-ADJUSTED BULK ANALYSIS\n")
cat("===========================================\n\n")
cat("Samples:", nrow(analysis_df), "\n")
cat("Monocytic abundance proxy:", monocyte_name, "\n\n")
cat("PRIMARY CORRELATIONS\n"); print(cor_results)
cat("\n\nPARTIAL CORRELATIONS\n"); print(partial_results)
cat("\n\nCOMPOSITION MODEL VARIANCE EXPLAINED\n"); print(variance_summary)
cat("\n\nRESIDUALIZED MODULE RESULTS\n"); print(resid_results)
cat("\n\nCORE GENE RESIDUAL RESULTS\n"); print(core_gene_results)
sink()

cat("\n=============================================\n")
cat("IDEA 3A COMPLETED\n")
cat("=============================================\n\n")
