###############################################################################
# GSE43292 — PRESPECIFIED PAIRED DISEASE-TRANSITION VALIDATION
# Repository copy prepared for the carotid-plaque-transcriptomics project.
# Original analysis logic retained; helper lookup made repository-relative.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260924)

GSE <- "GSE43292"
MODULE_FILE <- "Whole_transcriptome_module_assignments.csv"
HELPER_FILE <- if (file.exists("VALIDATION_COMMON_HELPERS_v4_FIXED.R")) {
  "VALIDATION_COMMON_HELPERS_v4_FIXED.R"
} else {
  file.path("scripts", "04_external_validation", "VALIDATION_COMMON_HELPERS_v4_FIXED.R")
}
OUTDIR <- "GSE43292_PRESPECIFIED_VALIDATION"
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
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
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

titles <- as.character(pd$title)
patient <- sub(".*Patient[[:space:]]+([0-9]+).*", "\\1", titles, ignore.case = TRUE)

condition <- ifelse(
  grepl("Atheroma plaque", titles, ignore.case = TRUE),
  "Plaque",
  ifelse(grepl("Macroscopically intact", titles, ignore.case = TRUE), "Intact", NA_character_)
)

sample_meta <- data.frame(
  Sample = rownames(pd), Title = titles, Patient = patient, Condition = condition,
  stringsAsFactors = FALSE
)

if (anyNA(sample_meta$Condition)) stop("Could not classify all samples. Inspect metadata output.")

pair_check <- sample_meta |>
  dplyr::count(Patient, Condition) |>
  tidyr::pivot_wider(names_from = Condition, values_from = n, values_fill = 0)

if (nrow(pair_check) != 32L || any(pair_check$Plaque != 1) || any(pair_check$Intact != 1))
  stop("Expected 32 complete pairs; reconstructed pairing is not complete.")

write.csv(sample_meta, file.path(OUTDIR, "02_GSE43292_pairing_table.csv"), row.names = FALSE)

###############################################################################
# FROZEN MANUSCRIPT PROBE-TO-GENE MAPPING
#
# The original manuscript run used the then-current Bioconductor annotation
# stack. Gene-symbol resolution can change between annotation-database
# releases, even when the chip annotation package version is unchanged.
#
# The exact final probe-to-gene assignments from the historical manuscript
# analysis are therefore frozen in the repository.
###############################################################################

probe_ids <- rownames(expr_probe)

FROZEN_MAP <- file.path(
  "data",
  "frozen",
  "GSE43292_GPL6244_probe_to_gene_frozen.csv"
)

if (!file.exists(FROZEN_MAP)) {
  stop(
    "Frozen GSE43292 probe-to-gene mapping is missing: ",
    FROZEN_MAP
  )
}

message(
  "Using frozen manuscript GSE43292 probe-to-gene mapping: ",
  FROZEN_MAP
)

frozen_map <- read.csv(
  FROZEN_MAP,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

if (!all(c("Probe", "Gene_symbol") %in% names(frozen_map))) {
  stop("Frozen GSE43292 mapping must contain Probe and Gene_symbol columns.")
}

frozen_map$Probe <- as.character(frozen_map$Probe)

if (nrow(frozen_map) != 33297L) {
  stop(
    "Frozen mapping row-count mismatch: expected 33297, observed ",
    nrow(frozen_map)
  )
}

if (anyDuplicated(frozen_map$Probe)) {
  stop("Frozen GSE43292 mapping contains duplicated probe IDs.")
}

idx <- match(
  probe_ids,
  frozen_map$Probe
)

if (anyNA(idx)) {
  stop(
    "Current GEO matrix contains ",
    sum(is.na(idx)),
    " probe(s) absent from the frozen GSE43292 mapping."
  )
}

symbols <- frozen_map$Gene_symbol[idx]

mapped_n <- sum(!is.na(symbols) & nzchar(symbols))
unmapped_n <- length(symbols) - mapped_n

message("Frozen mapping rows aligned: ", length(symbols))
message("Mapped probes: ", mapped_n)
message("Unmapped probes: ", unmapped_n)

if (mapped_n != 22024L) {
  stop(
    "Mapped-probe count mismatch: expected 22024, observed ",
    mapped_n
  )
}

if (unmapped_n != 11273L) {
  stop(
    "Unmapped-probe count mismatch: expected 11273, observed ",
    unmapped_n
  )
}

expr_gene <- collapse_to_gene(expr_probe, symbols)

write.csv(
  data.frame(
    Probe = probe_ids,
    Gene_symbol = symbols,
    check.names = FALSE
  ),
  file.path(OUTDIR, "03_probe_to_gene_mapping.csv"),
  row.names = FALSE,
  na = ""
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
delta_rows <- list(); k <- 1L
for (pat in unique(scores$Patient)) {
  z <- scores[scores$Patient == pat, , drop = FALSE]
  plaque <- z[z$Condition == "Plaque", , drop = FALSE]
  intact <- z[z$Condition == "Intact", , drop = FALSE]
  for (sc in score_names) {
    delta_rows[[k]] <- data.frame(
      Patient = pat, Score = sc, Intact = intact[[sc]][1], Plaque = plaque[[sc]][1],
      Plaque_minus_Intact = plaque[[sc]][1] - intact[[sc]][1]
    )
    k <- k + 1L
  }
}
deltas <- do.call(rbind, delta_rows)
write.csv(deltas, file.path(OUTDIR, "06_paired_patient_deltas.csv"), row.names = FALSE)

results <- do.call(rbind, lapply(score_names, function(sc) {
  x <- deltas[deltas$Score == sc, , drop = FALSE]
  wt <- suppressWarnings(wilcox.test(x$Plaque, x$Intact, paired = TRUE, exact = FALSE))
  tt <- t.test(x$Plaque, x$Intact, paired = TRUE)
  predicted <- if (sc == "Blue_WGCNA") "Plaque lower than Intact" else "Plaque higher than Intact"
  observed <- if (median(x$Plaque_minus_Intact, na.rm = TRUE) > 0)
    "Plaque higher than Intact" else "Plaque lower than Intact"
  data.frame(
    Score = sc, N_pairs = nrow(x),
    Median_Intact = median(x$Intact, na.rm = TRUE),
    Median_Plaque = median(x$Plaque, na.rm = TRUE),
    Median_delta = median(x$Plaque_minus_Intact, na.rm = TRUE),
    Positive_delta_pairs = sum(x$Plaque_minus_Intact > 0, na.rm = TRUE),
    Negative_delta_pairs = sum(x$Plaque_minus_Intact < 0, na.rm = TRUE),
    Wilcoxon_V = unname(wt$statistic), Wilcoxon_P = wt$p.value,
    Paired_t = unname(tt$statistic), Paired_t_P = tt$p.value,
    Predicted_direction = predicted, Observed_direction = observed,
    Direction_matches_prediction = predicted == observed
  )
}))

results$Wilcoxon_FDR <- p.adjust(results$Wilcoxon_P, method = "BH")
results$Paired_t_FDR <- p.adjust(results$Paired_t_P, method = "BH")
write.csv(results, file.path(OUTDIR, "07_prespecified_paired_tests.csv"), row.names = FALSE)

plot_df <- scores |>
  tidyr::pivot_longer(cols = dplyr::all_of(score_names), names_to = "Score", values_to = "Value")

p <- ggplot(plot_df, aes(x = Condition, y = Value, group = Patient)) +
  geom_line(alpha = 0.35) + geom_point(size = 1.7) +
  facet_wrap(~Score, scales = "free_y", nrow = 1) +
  theme_bw(base_size = 11) +
  labs(title = "GSE43292: paired validation of frozen plaque programmes", x = NULL, y = "Gene-set score")

ggsave(file.path(OUTDIR, "Figure_GSE43292_paired_scores.png"), p, width = 11, height = 4.5, dpi = 300)

sink(file.path(OUTDIR, "FINAL_GSE43292_VALIDATION_SUMMARY.txt"))
cat("GSE43292 — PRESPECIFIED PAIRED DISEASE-TRANSITION VALIDATION\n")
cat("================================================================\n\n")
cat("Patients:", length(unique(scores$Patient)), "\n")
cat("Samples:", nrow(scores), "\n")
cat("Design: paired advanced plaque vs distant macroscopically intact/early-stage tissue.\n\n")
cat("FROZEN HYPOTHESES\nBrown_WGCNA higher in Plaque\nInflammasome_65 higher in Plaque\nBlue_WGCNA lower in Plaque\n\n")
cat("GENE-SET COVERAGE\n"); print(coverage, row.names = FALSE)
cat("\nPRESPECIFIED PAIRED TESTS\n"); print(results, row.names = FALSE)
cat("\nINTERPRETATION\n")
cat("Primary evidence: paired Wilcoxon tests with BH correction across 3 frozen tests.\n")
cat("The intact comparator is not a truly healthy artery; it is distant stage-I/II tissue from patients with atheroma.\n")
cat("No module redefinition or post-hoc gene selection should be performed in this cohort.\n")
sink()

capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))
cat("\nGSE43292 COMPLETE\n")
