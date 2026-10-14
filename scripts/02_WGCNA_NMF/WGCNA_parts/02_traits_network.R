###############################################################################
# 7. LOAD THE LOCKED SIGNATURES USED IN THE MANUSCRIPT
#
# The original development script queried MSigDB at runtime. For the public
# reproducibility workflow the exact membership is frozen in signatures/ so a
# future msigdbr/MSigDB update cannot silently change the scored gene set.
###############################################################################

read_locked_signature <- function(path) {
  if (!file.exists(path)) stop("Missing locked signature file: ", path)
  x <- clean_gene(readLines(path, warn = FALSE))
  unique(x[!is.na(x) & nzchar(x)])
}

inflammasome_source <- read_locked_signature(
  file.path("signatures", "Inflammasome_source67.txt")
)
FIBROUS_GENES <- read_locked_signature(
  file.path("signatures", "Fibrous7.txt")
)
OSTEOGENIC_2GENE <- read_locked_signature(
  file.path("signatures", "RUNX2_SPP1.txt")
)

signature_source <- list(
  Inflammasome = inflammasome_source,
  Fibrous_remodelling = FIBROUS_GENES,
  Osteogenic_RUNX2_SPP1 = OSTEOGENIC_2GENE
)

signature_present <- lapply(
  signature_source,
  function(genes) intersect(genes, rownames(expr_full))
)

signature_missing <- lapply(
  signature_source,
  function(genes) setdiff(genes, rownames(expr_full))
)

signature_summary <- data.frame(
  Signature = names(signature_source),
  Source_genes = lengths(signature_source),
  Present_genes = lengths(signature_present),
  Missing_genes = lengths(signature_missing),
  Coverage_percent = round(
    100 * lengths(signature_present) / lengths(signature_source),
    1
  ),
  Missing_gene_names = vapply(
    signature_missing,
    paste,
    collapse = ";",
    FUN.VALUE = character(1)
  )
)

write_csv_safe(
  signature_summary,
  file.path(OUTDIR, "01_Trait_scores", "Signature_coverage.csv")
)

###############################################################################
# 8. CALCULATE THE THREE ssGSEA SCORES
###############################################################################

ssgsea_matrix <- run_ssgsea(expr_full, signature_present)

score_df <- data.frame(
  Sample = colnames(expr_full),
  t(ssgsea_matrix),
  check.names = FALSE
)

write_csv_safe(
  score_df,
  file.path(OUTDIR, "01_Trait_scores", "Samplewise_ssGSEA_scores.csv")
)

datTraits <- data.frame(
  Inflammasome = score_df$Inflammasome,
  Fibrous = score_df$Fibrous_remodelling,
  Osteogenic = score_df$Osteogenic_RUNX2_SPP1,
  row.names = score_df$Sample,
  check.names = FALSE
)

###############################################################################
# 9. TRANSPOSE TO WGCNA FORMAT
###############################################################################

datExpr <- as.data.frame(t(expr_full))

###############################################################################
# 10. QUALITY CHECK USING goodSamplesGenes()
###############################################################################

gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  bad_genes <- colnames(datExpr)[!gsg$goodGenes]
  bad_samples <- rownames(datExpr)[!gsg$goodSamples]

  write_csv_safe(
    data.frame(Gene = bad_genes),
    file.path(OUTDIR, "00_Input_QC", "WGCNA_removed_bad_genes.csv")
  )

  write_csv_safe(
    data.frame(Sample = bad_samples),
    file.path(OUTDIR, "00_Input_QC", "WGCNA_removed_bad_samples.csv")
  )

  datExpr <- datExpr[
    gsg$goodSamples,
    gsg$goodGenes,
    drop = FALSE
  ]
}

datTraits <- datTraits[rownames(datExpr), , drop = FALSE]

if (!identical(rownames(datExpr), rownames(datTraits))) {
  stop("Sample order in expression and trait data does not match.")
}

cat("\nWGCNA input:\n")
print(dim(datExpr))

###############################################################################
# 11. SAMPLE CLUSTERING
###############################################################################

sampleTree <- hclust(dist(datExpr), method = "average")

pdf(
  file.path(OUTDIR, "00_Input_QC", "Sample_clustering_dendrogram.pdf"),
  width = 10,
  height = 6
)
plot(
  sampleTree,
  main = "GSE125771 sample clustering before WGCNA",
  xlab = "",
  sub = "",
  cex = 0.75
)
dev.off()

###############################################################################
# 12. TEST SOFT-THRESHOLD POWERS
###############################################################################

cat("\nSelecting soft-thresholding power...\n")

sft <- WGCNA::pickSoftThreshold(
  datExpr,
  powerVector = POWERS,
  networkType = NETWORK_TYPE,
  corFnc = WGCNA::bicor,
  corOptions = list(maxPOutliers = MAX_P_OUTLIERS),
  verbose = 5
)

fit_table <- as.data.frame(sft$fitIndices)
fit_table$Signed_R2 <- -sign(fit_table[["slope"]]) * fit_table[["SFT.R.sq"]]

eligible_power <- fit_table[
  is.finite(fit_table$Signed_R2) &
    fit_table$Signed_R2 >= SFT_R2_TARGET,
  ,
  drop = FALSE
]

if (nrow(eligible_power) > 0) {
  softPower <- min(eligible_power[["Power"]])
  power_reason <- paste0(
    "Lowest tested signed-network power with signed scale-free R2 >= ",
    SFT_R2_TARGET
  )
} else {
  softPower <- FALLBACK_SIGNED_POWER_30_40
  power_reason <- paste0(
    "No tested power reached signed scale-free R2 >= ",
    SFT_R2_TARGET,
    "; using WGCNA FAQ fallback for a signed network with 30-40 samples."
  )
}

cat("\nSelected soft power:", softPower, "\n")
cat("Reason:", power_reason, "\n")

fit_table$Selected <- fit_table[["Power"]] == softPower

write_csv_safe(
  fit_table,
  file.path(OUTDIR, "02_Soft_threshold", "Soft_threshold_results.csv")
)

writeLines(
  c(
    paste("Selected_power:", softPower),
    paste("Reason:", power_reason)
  ),
  file.path(OUTDIR, "02_Soft_threshold", "Selected_power.txt")
)

###############################################################################
# 13. SOFT-THRESHOLD FIGURE
###############################################################################

pdf(
  file.path(OUTDIR, "02_Soft_threshold", "Soft_threshold_diagnostics.pdf"),
  width = 10,
  height = 5
)

par(mfrow = c(1, 2))

plot(
  fit_table[["Power"]],
  fit_table$Signed_R2,
  xlab = "Soft-thresholding power",
  ylab = "Signed scale-free topology fit (R2)",
  main = "Scale-free topology fit",
  type = "n"
)
text(fit_table[["Power"]], fit_table$Signed_R2, labels = fit_table[["Power"]])
abline(h = SFT_R2_TARGET, lty = 2)
abline(v = softPower, lty = 3)

plot(
  fit_table[["Power"]],
  fit_table[["mean(k)"]],
  xlab = "Soft-thresholding power",
  ylab = "Mean connectivity",
  main = "Mean connectivity",
  type = "n"
)
text(
  fit_table[["Power"]],
  fit_table[["mean(k)"]],
  labels = fit_table[["Power"]]
)
abline(v = softPower, lty = 3)

dev.off()

###############################################################################
# 14. BUILD SIGNED WGCNA NETWORK
###############################################################################

cat("\nBuilding whole-transcriptome WGCNA network...\n")

net <- WGCNA::blockwiseModules(
  datExpr,
  power = softPower,
  corType = COR_TYPE,
  maxPOutliers = MAX_P_OUTLIERS,
  networkType = NETWORK_TYPE,
  TOMType = TOM_TYPE,
  maxBlockSize = MAX_BLOCK_SIZE,
  minModuleSize = MIN_MODULE_SIZE,
  deepSplit = DEEP_SPLIT,
  mergeCutHeight = MERGE_CUT_HEIGHT,
  pamRespectsDendro = FALSE,
  reassignThreshold = 0,
  numericLabels = FALSE,
  saveTOMs = FALSE,
  randomSeed = SEED,
  useCorOptionsThroughout = TRUE,
  verbose = 3
)

moduleColors <- as.character(net$colors)
names(moduleColors) <- colnames(datExpr)
