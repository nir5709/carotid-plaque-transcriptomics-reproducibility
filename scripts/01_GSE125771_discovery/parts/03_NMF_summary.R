###############################################################################
# 6. NMF ON THE STRICT INFLAMMASOME TOP10-vs-BOTTOM10 DEGs
#
# The manuscript-reported hard assignment is the k=2 solution. Components are
# renamed after fitting so NMF_Module_1 is the component most positively
# associated with the continuous Inflammasome score. This relabelling does not
# alter the NMF solution.
###############################################################################

if (!requireNamespace("NMF", quietly = TRUE)) {
  install.packages("NMF", repos = "https://cloud.r-project.org")
}

# NMF's multi-run backend may call path.package() internally. Attach NMF
# explicitly so clean non-interactive Rscript runs work reliably.
suppressPackageStartupMessages(library(NMF))

inflam_degs <- unique(extreme_results$Inflammasome$significant$Gene)
nmf_matrix <- expr_full[intersect(inflam_degs, rownames(expr_full)), , drop = FALSE]
nmf_matrix <- nmf_matrix[apply(nmf_matrix, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

if (nrow(nmf_matrix) < 30L) stop("Too few strict Inflammasome DEGs for NMF.")
if (min(nmf_matrix) < 0) nmf_matrix <- nmf_matrix - min(nmf_matrix) + 1e-8

write.csv(
  data.frame(Gene=rownames(nmf_matrix)),
  file.path(OUTDIR, "06_NMF", "NMF_input_inflammasome_DEGs.csv"),
  row.names = FALSE
)

# Optional rank survey retained for transparency; the final hard assignment is k=2.
# Force sequential execution for portability across clean Rscript environments.
rank_rows <- list()
for (k in 2:6) {
  fit_k <- NMF::nmf(
    nmf_matrix,
    rank = k,
    method = "brunet",
    nrun = 20,
    seed = 20260722 + k,
    .options = "vp1"
  )
  sm <- summary(fit_k)
  Hk <- NMF::coef(fit_k)
  groups_k <- apply(Hk, 2, which.max)
  rank_rows[[as.character(k)]] <- data.frame(
    Rank = k,
    Cophenetic = if ("cophenetic" %in% names(sm)) as.numeric(sm[["cophenetic"]]) else NA_real_,
    Dispersion = if ("dispersion" %in% names(sm)) as.numeric(sm[["dispersion"]]) else NA_real_,
    RSS = if ("rss" %in% names(sm)) as.numeric(sm[["rss"]]) else NA_real_,
    Minimum_sample_group = min(tabulate(groups_k, nbins = k)),
    stringsAsFactors = FALSE
  )
}
rank_metrics <- do.call(rbind, rank_rows)
rank_metrics$Manuscript_primary_rank <- rank_metrics$Rank == NMF_FINAL_RANK
write.csv(rank_metrics, file.path(OUTDIR, "06_NMF", "NMF_rank_survey.csv"), row.names = FALSE)

final_fit <- NMF::nmf(
  nmf_matrix,
  rank = NMF_FINAL_RANK,
  method = "brunet",
  nrun = NMF_NRUN_FINAL,
  seed = 20260722 + 20000 + NMF_FINAL_RANK,
  .options = "vp1"
)

W <- NMF::basis(final_fit)
H <- NMF::coef(final_fit)

provisional_scores <- t(sapply(seq_len(ncol(W)), function(i) {
  weights <- W[, i]
  as.numeric(crossprod(weights / sum(weights), nmf_matrix))
}))
colnames(provisional_scores) <- colnames(nmf_matrix)
inflam_score <- score_df$Inflammasome[match(colnames(nmf_matrix), score_df$Sample)]
component_rho <- apply(provisional_scores, 1, function(x) cor(x, inflam_score, method = "spearman"))
ord <- order(component_rho, decreasing = TRUE, na.last = TRUE)
W <- W[, ord, drop = FALSE]
H <- H[ord, , drop = FALSE]
colnames(W) <- paste0("NMF_Module_", seq_len(ncol(W)))
rownames(H) <- colnames(W)

hard_module <- apply(W, 1, which.max)
gene_modules <- data.frame(
  Gene = rownames(W),
  Module = colnames(W)[hard_module],
  stringsAsFactors = FALSE
)

# Preserve the simple two-column file expected by the WGCNA overlap workflow.
write.csv(gene_modules, file.path(OUTDIR, "06_NMF", "NMF_gene_modules.csv"), row.names = FALSE)
write.csv(gene_modules, "NMF_gene_modules.csv", row.names = FALSE)

module_counts <- as.data.frame(table(gene_modules$Module), stringsAsFactors = FALSE)
names(module_counts) <- c("Module", "Genes")
write.csv(module_counts, file.path(OUTDIR, "06_NMF", "NMF_module_gene_counts.csv"), row.names = FALSE)

sample_group_index <- apply(H, 2, which.max)
sample_groups <- data.frame(
  Sample = colnames(H),
  NMF_group = paste0("NMF_Group_", sample_group_index),
  stringsAsFactors = FALSE
)
write.csv(sample_groups, file.path(OUTDIR, "06_NMF", "NMF_sample_groups.csv"), row.names = FALSE)

module_scores <- t(sapply(seq_len(ncol(W)), function(i) {
  weights <- W[, i]
  as.numeric(crossprod(weights / sum(weights), nmf_matrix))
}))
rownames(module_scores) <- colnames(W)
colnames(module_scores) <- colnames(nmf_matrix)
module_score_df <- data.frame(Sample=colnames(module_scores), t(module_scores), check.names=FALSE)
module_score_df <- merge(module_score_df, score_df, by="Sample", all.x=TRUE, sort=FALSE)

nmf_cor <- do.call(rbind, lapply(colnames(W), function(m) {
  do.call(rbind, lapply(signatures, function(sig) {
    ct <- safe_cor(module_score_df[[m]], module_score_df[[sig]])
    data.frame(Module=m, Signature=sig, Spearman_rho=ct$rho,
               P_value=ct$P_value, N=ct$N, stringsAsFactors=FALSE)
  }))
}))
nmf_cor$FDR <- p.adjust(nmf_cor$P_value, method="BH")
write.csv(nmf_cor, file.path(OUTDIR, "06_NMF", "NMF_module_signature_correlations.csv"), row.names = FALSE)
write.csv(data.frame(Gene=rownames(W), W, check.names=FALSE),
          file.path(OUTDIR, "06_NMF", "NMF_W_basis.csv"), row.names = FALSE)
write.csv(data.frame(Module=rownames(H), H, check.names=FALSE),
          file.path(OUTDIR, "06_NMF", "NMF_H_coefficients.csv"), row.names = FALSE)

###############################################################################
# 7. FINAL AUDIT SUMMARY AND SESSION INFORMATION
###############################################################################

summary_file <- file.path(OUTDIR, "FINAL_GSE125771_DISCOVERY_SUMMARY.txt")
sink(summary_file)
cat("GSE125771 DISCOVERY ANALYSIS\n")
cat("============================\n\n")
cat("Samples:", ncol(expr_full), "\n")
cat("Genes:", nrow(expr_full), "\n\n")
cat("LOCKED SIGNATURE COVERAGE\n"); print(coverage, row.names=FALSE)
cat("\nSIGNATURE CORRELATIONS\n"); print(score_correlations, row.names=FALSE)
cat("\nTOP10-vs-BOTTOM10 DEG SUMMARY\n"); print(deg_summary, row.names=FALSE)
cat("\nNMF INPUT GENES:", nrow(nmf_matrix), "\n")
cat("NMF FINAL RANK:", NMF_FINAL_RANK, "\n")
cat("\nNMF HARD-ASSIGNMENT COUNTS\n"); print(module_counts, row.names=FALSE)
if (exists("mcp_cor")) {
  cat("\nMCP-COUNTER CORRELATIONS\n")
  print(mcp_cor, row.names=FALSE)
}
if (exists("partial_table")) {
  cat("\nMCP-COUNTER COMPOSITION-ADJUSTED SENSITIVITY\n")
  print(partial_table, row.names=FALSE)
}
sink()

capture.output(sessionInfo(), file=file.path(OUTDIR, "07_reproducibility", "sessionInfo.txt"))

save(
  expr_full, score_df, group_tables, score_correlations,
  extreme_results, deg_summary, overlap, gene_modules, module_scores,
  file=file.path(OUTDIR, "07_reproducibility", "GSE125771_discovery_core_objects.RData")
)

cat("GSE125771 discovery workflow complete.\n")
