###############################################################################
# GSE163154 — PRESPECIFIED IPH / PLAQUE-INSTABILITY VALIDATION
# Repository copy prepared for the carotid-plaque-transcriptomics project.
# The primary 43-segment analysis is segment-level and not patient-adjusted
# unless an explicit sample-to-patient map is supplied.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260924)

GSE <- "GSE163154"
MODULE_FILE <- "Whole_transcriptome_module_assignments.csv"
HELPER_FILE <- if (file.exists("VALIDATION_COMMON_HELPERS_v4_FIXED.R")) {
  "VALIDATION_COMMON_HELPERS_v4_FIXED.R"
} else {
  file.path("scripts", "04_external_validation", "VALIDATION_COMMON_HELPERS_v4_FIXED.R")
}
PATIENT_MAP_FILE <- "GSE163154_patient_map.csv"
FROZEN_MAP <- file.path(
  "data", "frozen", "GSE163154_GPL6104_probe_to_gene_frozen.csv"
)
OUTDIR <- "GSE163154_PRESPECIFIED_VALIDATION"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(HELPER_FILE)) stop("Missing helper file: ", HELPER_FILE)
source(HELPER_FILE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

cran_pkgs <- c("ggplot2","dplyr","tidyr")
bioc_pkgs <- c("GEOquery","limma")

for (p in cran_pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
for (p in bioc_pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, ask = FALSE, update = FALSE)

missing <- c(cran_pkgs, bioc_pkgs)[
  !vapply(c(cran_pkgs, bioc_pkgs), requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(GEOquery); library(limma)
  library(ggplot2); library(dplyr); library(tidyr)
})

module_path <- find_module_file(MODULE_FILE)
if (is.na(module_path) || !file.exists(module_path))
  stop("Required module assignment file not found: ", MODULE_FILE)
mods <- read_brown_blue(module_path)
brown_genes <- mods$brown
blue_genes <- mods$blue

message("Loading ", GSE, " from GEO...")
g <- GEOquery::getGEO(GSE, GSEMatrix = TRUE, getGPL = FALSE)
eset <- if (inherits(g, "ExpressionSet")) g else g[[1]]
expr_probe <- Biobase::exprs(eset)
pd <- Biobase::pData(eset)
write.csv(pd, file.path(OUTDIR, "01_GEO_sample_metadata_raw.csv"), row.names = TRUE)

titles <- trimws(as.character(pd$title))
condition <- rep(NA_character_, length(titles))
condition[grepl("^non-IPH[[:space:]]*[0-9]+$", titles, ignore.case = TRUE)] <- "non_IPH"
condition[grepl("^IPH[[:space:]]*[0-9]+$", titles, ignore.case = TRUE)] <- "IPH"

sample_meta <- data.frame(
  Sample = rownames(pd), Title = titles, Condition = condition,
  stringsAsFactors = FALSE
)

cat("\nDetected GSE163154 samples:\n")
print(sample_meta[, c("Sample","Title","Condition")])
cat("\nCondition counts:\n")
print(table(sample_meta$Condition, useNA = "ifany"))

if (anyNA(sample_meta$Condition)) {
  bad <- sample_meta[is.na(sample_meta$Condition), , drop = FALSE]
  write.csv(bad, file.path(OUTDIR, "ERROR_unclassified_GSE163154_samples.csv"), row.names = FALSE)
  stop("One or more samples could not be classified. See ERROR_unclassified_GSE163154_samples.csv.")
}

write.csv(sample_meta, file.path(OUTDIR, "02_GSE163154_group_table.csv"), row.names = FALSE)
group_counts <- table(sample_meta$Condition)
if (!all(c("IPH","non_IPH") %in% names(group_counts)))
  stop("IPH/non-IPH groups were not reconstructed correctly.")
if (as.integer(group_counts["IPH"]) != 27L || as.integer(group_counts["non_IPH"]) != 16L)
  stop(paste0("Unexpected group counts. Expected IPH=27 and non_IPH=16; detected: ",
              paste(names(group_counts), as.integer(group_counts), sep = "=", collapse = ", ")))

probe_ids <- rownames(expr_probe)

if (!file.exists(FROZEN_MAP))
  stop("Frozen GSE163154 probe-to-gene mapping not found: ", FROZEN_MAP)

message("Using frozen manuscript GSE163154 probe-to-gene mapping: ", FROZEN_MAP)

frozen_map <- read.csv(
  FROZEN_MAP,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

if (!all(c("Probe", "Gene_symbol") %in% names(frozen_map)))
  stop("Frozen mapping must contain columns Probe and Gene_symbol.")

if (nrow(frozen_map) != 22184L)
  stop("Unexpected frozen mapping row count: ", nrow(frozen_map))

if (anyDuplicated(frozen_map$Probe))
  stop("Frozen mapping contains duplicate probe IDs.")

idx <- match(probe_ids, frozen_map$Probe)

if (anyNA(idx))
  stop(sum(is.na(idx)), " current GEO probes are absent from the frozen mapping.")

symbols <- frozen_map$Gene_symbol[idx]

mapped_n <- sum(!is.na(symbols) & nzchar(trimws(symbols)))
unmapped_n <- length(symbols) - mapped_n

if (mapped_n != 21199L || unmapped_n != 985L)
  stop(
    "Unexpected frozen mapping counts. Expected mapped=21199 and unmapped=985; detected mapped=",
    mapped_n, " and unmapped=", unmapped_n
  )

message("Frozen mapping rows aligned: ", length(symbols))
message("Mapped probes: ", mapped_n)
message("Unmapped probes: ", unmapped_n)

expr_gene <- collapse_to_gene(expr_probe, symbols)

write.csv(
  data.frame(Probe = probe_ids, Gene_symbol = symbols),
  file.path(OUTDIR, "03_probe_to_gene_mapping.csv"),
  row.names = FALSE
)

sb <- score_gene_set(expr_gene, brown_genes)
sl <- score_gene_set(expr_gene, blue_genes)
si <- score_gene_set(expr_gene, inflammasome_65)
coverage <- data.frame(
  GeneSet = c("Brown_WGCNA","Blue_WGCNA","Inflammasome_65"),
  Total_genes = c(length(brown_genes), length(blue_genes), length(inflammasome_65)),
  Present_genes = c(length(sb$present), length(sl$present), length(si$present))
)
coverage$Coverage_fraction <- coverage$Present_genes / coverage$Total_genes
write.csv(coverage, file.path(OUTDIR, "04_gene_set_coverage.csv"), row.names = FALSE)

scores <- sample_meta
scores$Brown_WGCNA <- sb$score[match(scores$Sample, colnames(expr_gene))]
scores$Blue_WGCNA <- sl$score[match(scores$Sample, colnames(expr_gene))]
scores$Inflammasome_65 <- si$score[match(scores$Sample, colnames(expr_gene))]
write.csv(scores, file.path(OUTDIR, "05_sample_scores.csv"), row.names = FALSE)

score_names <- c("Brown_WGCNA","Blue_WGCNA","Inflammasome_65")
results <- do.call(rbind, lapply(score_names, function(sc) {
  iph <- scores[scores$Condition == "IPH", sc, drop = TRUE]
  non <- scores[scores$Condition == "non_IPH", sc, drop = TRUE]
  wt <- suppressWarnings(wilcox.test(iph, non, paired = FALSE, exact = FALSE))
  tt <- t.test(iph, non, var.equal = FALSE)
  predicted <- if (sc == "Blue_WGCNA") "IPH lower than non-IPH" else "IPH higher than non-IPH"
  observed <- if (median(iph, na.rm = TRUE) > median(non, na.rm = TRUE))
    "IPH higher than non-IPH" else "IPH lower than non-IPH"
  data.frame(
    Score = sc,
    N_IPH = sum(is.finite(iph)), N_non_IPH = sum(is.finite(non)),
    Median_IPH = median(iph, na.rm = TRUE),
    Median_non_IPH = median(non, na.rm = TRUE),
    Median_difference = median(iph, na.rm = TRUE) - median(non, na.rm = TRUE),
    Wilcoxon_W = unname(wt$statistic), Wilcoxon_P = wt$p.value,
    Welch_t = unname(tt$statistic), Welch_t_P = tt$p.value,
    Predicted_direction = predicted, Observed_direction = observed,
    Direction_matches_prediction = predicted == observed
  )
}))
results$Wilcoxon_FDR <- p.adjust(results$Wilcoxon_P, method = "BH")
results$Welch_t_FDR <- p.adjust(results$Welch_t_P, method = "BH")
write.csv(results, file.path(OUTDIR, "06_prespecified_IPH_tests.csv"), row.names = FALSE)

blocked_results <- NULL
if (file.exists(PATIENT_MAP_FILE)) {
  message("Patient map found; running patient-blocked limma sensitivity.")
  pm <- read.csv(PATIENT_MAP_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("Sample","Patient") %in% names(pm)))
    stop(PATIENT_MAP_FILE, " must contain columns Sample and Patient.")
  scores2 <- merge(scores, pm[, c("Sample","Patient")], by = "Sample", all.x = TRUE, sort = FALSE)
  if (anyNA(scores2$Patient)) stop("Patient map is missing one or more GEO samples.")
  score_mat <- t(as.matrix(scores2[, score_names, drop = FALSE]))
  colnames(score_mat) <- scores2$Sample
  group <- factor(scores2$Condition, levels = c("non_IPH","IPH"))
  design <- model.matrix(~ group)
  block <- factor(scores2$Patient)
  dup <- limma::duplicateCorrelation(score_mat, design, block = block)
  fit <- limma::lmFit(score_mat, design, block = block, correlation = dup$consensus.correlation)
  fit <- limma::eBayes(fit, robust = TRUE)
  blocked_results <- limma::topTable(fit, coef = "groupIPH", number = Inf, sort.by = "none")
  blocked_results$Score <- rownames(blocked_results); rownames(blocked_results) <- NULL
  blocked_results$BH_FDR <- p.adjust(blocked_results$P.Value, method = "BH")
  write.csv(blocked_results, file.path(OUTDIR, "07_patient_blocked_limma_sensitivity.csv"), row.names = FALSE)
}

plot_df <- scores |>
  tidyr::pivot_longer(cols = dplyr::all_of(score_names), names_to = "Score", values_to = "Value")
p <- ggplot(plot_df, aes(x = Condition, y = Value)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, alpha = 0.7, size = 1.6) +
  facet_wrap(~Score, scales = "free_y", nrow = 1) + theme_bw(base_size = 11) +
  labs(title = "GSE163154: frozen programme scores in IPH and non-IPH plaque segments",
       x = NULL, y = "Gene-set score")
ggsave(file.path(OUTDIR, "Figure_GSE163154_IPH_scores.png"), p, width = 11, height = 4.5, dpi = 300)

sink(file.path(OUTDIR, "FINAL_GSE163154_VALIDATION_SUMMARY.txt"))
cat("GSE163154 — PRESPECIFIED IPH VALIDATION\n")
cat("=======================================\n\n")
cat("Samples:", nrow(scores), "\nIPH:", sum(scores$Condition == "IPH"),
    "\nnon-IPH:", sum(scores$Condition == "non_IPH"), "\n\n")
cat("FROZEN HYPOTHESES\nBrown_WGCNA higher in IPH\nInflammasome_65 higher in IPH\nBlue_WGCNA lower in IPH\n\n")
cat("GENE-SET COVERAGE\n"); print(coverage, row.names = FALSE)
cat("\nPRESPECIFIED IPH TESTS\n"); print(results, row.names = FALSE)
if (!is.null(blocked_results)) {
  cat("\nOPTIONAL PATIENT-BLOCKED LIMMA SENSITIVITY\n")
  print(blocked_results, row.names = FALSE)
}
cat("\nINTERPRETATION\n")
cat("Primary evidence: Wilcoxon IPH-vs-non-IPH tests with BH correction across 3 frozen tests.\n")
cat("Public GEO metadata does not expose a clean patient-ID mapping for all retained transcriptomes; the primary 43-sample analysis is therefore not described as fully paired.\n")
cat("Association with IPH is histopathological association, not prospective risk prediction.\n")
cat("No module redefinition or post-hoc gene selection should be performed in this cohort.\n")
sink()

capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))
cat("\nGSE163154 COMPLETE\n")
