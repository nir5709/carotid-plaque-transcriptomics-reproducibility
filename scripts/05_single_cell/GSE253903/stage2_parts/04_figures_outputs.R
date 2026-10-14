###############################################################################
# 16. FIGURES
###############################################################################

save_boxplot(
  scores,
  "Brown_WGCNA",
  "GSE253903: Brown inflammatory WGCNA module by cell type",
  file.path(OUTDIR, "Figure_Brown_by_celltype")
)

save_boxplot(
  scores,
  "Blue_WGCNA",
  "GSE253903: Blue fibrous WGCNA module by cell type",
  file.path(OUTDIR, "Figure_Blue_by_celltype")
)

save_boxplot(
  scores,
  "Inflammasome_65",
  "GSE253903: locked 65-gene inflammasome score by cell type",
  file.path(OUTDIR, "Figure_Inflammasome65_by_celltype")
)

heat <- score_summary[, c("Broad_label","Brown_median","Blue_median","Inflammasome65_median")]
rownames(heat) <- heat$Broad_label
heat$Broad_label <- NULL

png(
  file.path(OUTDIR,"Figure_pseudobulk_module_heatmap.png"),
  width = 1800,
  height = 1400,
  res = 180
)
pheatmap::pheatmap(
  as.matrix(heat),
  scale = "column",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "GSE253903 donor-pseudobulk module localization"
)
dev.off()

pdf(
  file.path(OUTDIR,"Figure_pseudobulk_module_heatmap.pdf"),
  width = 8,
  height = 7
)
pheatmap::pheatmap(
  as.matrix(heat),
  scale = "column",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "GSE253903 donor-pseudobulk module localization"
)
dev.off()

###############################################################################
# 17. COVERAGE TABLE
###############################################################################

coverage <- data.frame(
  GeneSet = c("Brown_WGCNA","Blue_WGCNA","Inflammasome_65"),
  Total_genes = c(length(brown_genes),length(blue_genes),length(inflammasome_65)),
  Present_in_filtered_pseudobulk = c(
    length(intersect(brown_genes,rownames(logcpm))),
    length(intersect(blue_genes,rownames(logcpm))),
    length(intersect(inflammasome_65,rownames(logcpm)))
  ),
  stringsAsFactors = FALSE
)
coverage$Coverage_fraction <- coverage$Present_in_filtered_pseudobulk / coverage$Total_genes

write.csv(
  coverage,
  file.path(OUTDIR,"08_gene_set_coverage.csv"),
  row.names = FALSE
)

###############################################################################
# 18. SAVE ANNOTATED OBJECT
###############################################################################

saveRDS(
  obj,
  file.path(OUTDIR,"GSE253903_STAGE2_annotated_object.rds"),
  compress = FALSE
)

###############################################################################
# 19. FINAL SUMMARY
###############################################################################

summary_file <- file.path(OUTDIR,"FINAL_GSE253903_STAGE2_SUMMARY.txt")
sink(summary_file)
cat("GSE253903 STAGE 2 — CELL-TYPE LOCALIZATION VALIDATION\n")
cat("====================================================\n\n")
cat("Stage-1 object:",rds_path,"\n")
cat("Module file:",module_path,"\n")
cat("Total Stage-1 cells:",ncol(obj),"\n")
cat("Primary-analysis cells after ambiguous/mixed exclusions:",ncol(obj_primary),"\n")
cat("Donor/sample column:",sample_col,"\n")
cat("Pseudobulks retained (>=" , MIN_CELLS_PER_PSEUDOBULK, " cells): ", nrow(scores), "\n\n", sep = "")
cat("GENE-SET COVERAGE\n")
print(coverage,row.names=FALSE)
cat("\n\nLOCKED BROAD CELL-TYPE COUNTS\n")
print(annotation_counts,row.names=FALSE)
cat("\n\nCELL-TYPE SCORE SUMMARY\n")
print(score_summary,row.names=FALSE)
cat("\n\nPRIMARY DONOR-BALANCED TESTS\n")
print(paired_tests,row.names=FALSE)
cat("\n\nLOCALIZATION RANKING\n")
print(ranking,row.names=FALSE)
cat("\n\nINTERPRETATION RULES\n")
cat("1. Brown-module validation is strongest if Macrophage/Monocyte rank near the top and the donor-balanced Brown target-vs-other delta is consistently positive.\n")
cat("2. Blue-module validation is strongest if SMC/Fibroblast rank near the top and the donor-balanced Blue target-vs-other delta is consistently positive.\n")
cat("3. Pericyte/mural cells are kept separate from SMCs in the primary analysis and may be discussed as secondary structural localization.\n")
cat("4. Cluster 13 (cycling immune) and cluster 23 (SMC–macrophage mixed/doublet-like) are excluded from primary validation.\n")
cat("5. Cell-level observations are not treated as independent biological replicates; donor × cell-type pseudobulk is the primary evidence.\n")
cat("6. This analysis validates cell-type localization, not clinical prognosis or causality.\n")
sink()

capture.output(
  sessionInfo(),
  file = file.path(OUTDIR,"sessionInfo.txt")
)

cat("\n============================================================\n")
cat("GSE253903 STAGE 2 COMPLETE\n")
cat("============================================================\n\n")
cat("Please upload:\n")
cat("1. FINAL_GSE253903_STAGE2_SUMMARY.txt\n")
cat("2. 03_pseudobulk_gene_set_scores.csv\n")
cat("3. 04_celltype_score_summary.csv\n")
cat("4. 05_primary_target_vs_other_donor_deltas.csv\n")
cat("5. 06_primary_target_vs_other_paired_tests.csv\n")
cat("6. 07_module_localization_ranking.csv\n")
cat("7. 08_gene_set_coverage.csv\n")
cat("8. Figure_pseudobulk_module_heatmap.png\n")
cat("9. Figure_Brown_by_celltype.png\n")
cat("10. Figure_Blue_by_celltype.png\n")
