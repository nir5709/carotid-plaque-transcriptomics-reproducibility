###############################################################################
# 5. SIGNATURE COVERAGE AND ssGSEA
###############################################################################

coverage_21545 <- do.call(rbind, lapply(names(SIGNATURES_21545), function(nm) {
  gs <- unique(clean_gene(SIGNATURES_21545[[nm]]))
  present <- intersect(gs, rownames(expr_gene_21545))
  missing <- setdiff(gs, rownames(expr_gene_21545))
  data.frame(
    Signature = nm,
    Total_genes = length(gs),
    Present_genes = length(present),
    Missing_genes = length(missing),
    Percent_present = 100 * length(present) / length(gs),
    Missing_gene_names = paste(missing, collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
write.csv(coverage_21545, file.path(OUTDIR, "04_signature_coverage.csv"), row.names = FALSE)

# Only genes actually represented after probe-level mapping contribute to ssGSEA.
sets_present <- lapply(SIGNATURES_21545, intersect, y = rownames(expr_gene_21545))
scores_21545 <- run_ssgsea(expr_gene_21545, sets_present)
score_df_21545 <- data.frame(
  GSM = colnames(scores_21545),
  Inflammasome = as.numeric(scores_21545["Inflammasome", ]),
  Fibrous = as.numeric(scores_21545["Fibrous", ]),
  Osteogenic = as.numeric(scores_21545["Osteogenic", ]),
  stringsAsFactors = FALSE
)
write.csv(score_df_21545, file.path(OUTDIR, "05_samplewise_ssGSEA.csv"), row.names = FALSE)

###############################################################################
# 6. THREE PRESPECIFIED PAIRWISE SPEARMAN CORRELATIONS
###############################################################################

pairs <- list(
  c("Inflammasome", "Fibrous"),
  c("Inflammasome", "Osteogenic"),
  c("Fibrous", "Osteogenic")
)

stats <- do.call(rbind, lapply(pairs, function(pr) {
  ct <- suppressWarnings(
    cor.test(
      score_df_21545[[pr[1]]],
      score_df_21545[[pr[2]]],
      method = "spearman",
      exact = FALSE
    )
  )
  data.frame(
    X = pr[1], Y = pr[2],
    N = sum(complete.cases(score_df_21545[, pr])),
    Spearman_rho = unname(ct$estimate),
    P_value = ct$p.value,
    stringsAsFactors = FALSE
  )
}))
stats$FDR_BH_3correlations <- p.adjust(stats$P_value, method = "BH")
write.csv(stats, file.path(OUTDIR, "06_signature_correlations.csv"), row.names = FALSE)

###############################################################################
# 7. FINAL SUMMARY
###############################################################################

sink(file.path(OUTDIR, "FINAL_GSE21545_SUMMARY.txt"))
cat("GSE21545 EXTERNAL FIXED-SIGNATURE VALIDATION\n")
cat("============================================\n\n")
cat("Plaques:", nrow(score_df_21545), "\n\n")
cat("SIGNATURE COVERAGE AFTER ACTUAL GPL570 PROBE MAPPING\n")
print(coverage_21545, row.names = FALSE)
cat("\nPAIRWISE SPEARMAN CORRELATIONS\n")
print(stats, row.names = FALSE)
cat("\nINTERPRETATION\n")
cat("1. The three signature relationships were prespecified and BH-adjusted as one three-test family.\n")
cat("2. Probe-to-gene assignments use the frozen mapping derived from the final BioMart-assisted GPL570 annotation; only genes represented by mapped probes contribute to ssGSEA.\n")
cat("3. GSE21545 is used as an external molecular comparison; participant independence from GSE125771 is not assumed because both derive from the broader BiKE resource.\n")
sink()

capture.output(sessionInfo(), file=file.path(OUTDIR, "sessionInfo.txt"))

cat("GSE21545 fixed-signature validation complete.\n")
