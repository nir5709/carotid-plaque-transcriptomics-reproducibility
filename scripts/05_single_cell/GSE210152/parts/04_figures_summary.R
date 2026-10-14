###############################################################################
# 24. FIGURES
###############################################################################

p_mnp <- Seurat::DimPlot(
  mnp,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("GSE210152 MNP/macrophage-state clustering")

save_plot_both(
  p_mnp,
  file.path(OUTDIR, "Figure_MNP_UMAP_clusters"),
  width = 9,
  height = 7
)

state_marker_vec <- intersect(unique(unlist(state_sets)), rownames(mnp))
p_state_dot <- Seurat::DotPlot(
  mnp,
  features = state_marker_vec,
  group.by = "seurat_clusters",
  dot.scale = 5
) + RotatedAxis() + ggtitle("Published macrophage-state markers in GSE210152")

save_plot_both(
  p_state_dot,
  file.path(OUTDIR, "Figure_MNP_state_marker_dotplot"),
  width = 18,
  height = 8
)

candidate_features <- intersect(
  c(
    "PLIN2","TREM1","TREM2","IL1B","NLRP3","PYCARD",
    "CASP1","GBP5","C1QA","HMOX1","S100A8"
  ),
  rownames(mnp)
)

if (length(candidate_features) > 0L) {
  p_feature <- Seurat::FeaturePlot(
    mnp,
    features = candidate_features,
    reduction = "umap",
    ncol = 3,
    order = TRUE
  )
  save_plot_both(
    p_feature,
    file.path(OUTDIR, "Figure_MNP_key_gene_FeaturePlots"),
    width = 12,
    height = 12
  )
}

###############################################################################
# 25. SUMMARY TABLES FOR CANDIDATE CLUSTER
###############################################################################

normalize_cluster_label <- function(x) {
  sub("^g(?=[0-9]+$)", "", as.character(x), perl = TRUE)
}

candidate_state_scores <- state_scores[
  normalize_cluster_label(state_scores$Cluster) == candidate_cluster,
  , drop = FALSE
]

candidate_enrichment <- enrichment[
  normalize_cluster_label(enrichment$Cluster) == candidate_cluster,
  , drop = FALSE
]

candidate_key_avg <- key_avg[
  normalize_cluster_label(key_avg$Cluster) == candidate_cluster,
  , drop = FALSE
]

write.csv(
  candidate_state_scores,
  file.path(OUTDIR, "18_candidate_iLAM_state_scores.csv"),
  row.names = FALSE
)
write.csv(
  candidate_enrichment,
  file.path(OUTDIR, "19_candidate_iLAM_gene_set_enrichment.csv"),
  row.names = FALSE
)
write.csv(
  candidate_key_avg,
  file.path(OUTDIR, "20_candidate_iLAM_key_gene_expression.csv"),
  row.names = FALSE
)

###############################################################################
# 26. SAVE ANALYSIS OBJECT
###############################################################################

saveRDS(
  mnp,
  file.path(OUTDIR, "GSE210152_MNP_clustered_scored.rds"),
  compress = FALSE
)

###############################################################################
# 27. FINAL SUMMARY
###############################################################################

summary_file <- file.path(OUTDIR, "FINAL_GSE210152_MACROPHAGE_STATE_SUMMARY.txt")

sink(summary_file)
cat("GSE210152 — MACROPHAGE-STATE VALIDATION\n")
cat("========================================\n\n")
cat("Input:", basename(H5AD), "\n")
cat("H5AD assay used:", chosen_assay, "\n")
cat("Patient metadata source:", patient_info$source, "\n")
cat("Total immune cells retained:", ncol(obj), "\n")
cat("MNP cells analysed:", ncol(mnp), "\n")
cat("Number of MNP clusters:", length(unique(mnp$seurat_clusters)), "\n")
cat("Candidate PLIN2/TREM1 inflammatory-LAM-like cluster:", candidate_cluster, "\n\n")
cat("PATIENT COUNTS IN DATASET\n")
print(table(obj$patient))
cat("\n\nSYMPTOM STATUS\n")
print(table(obj$patient, obj$symptom_status))
cat("\n\nCANDIDATE iLAM-LIKE CLUSTER — DONOR REPRESENTATION\n")
print(candidate_donor, row.names = FALSE)
cat("\n\nCANDIDATE iLAM-LIKE CLUSTER — PUBLISHED STATE SCORES\n")
print(candidate_state_scores, row.names = FALSE)
cat("\n\nCANDIDATE iLAM-LIKE CLUSTER — KEY GENE EXPRESSION\n")
print(candidate_key_avg, row.names = FALSE)
cat("\n\nCANDIDATE iLAM-LIKE CLUSTER — GENE-SET ENRICHMENT\n")
print(candidate_enrichment, row.names = FALSE)
cat("\n\nDONOR-AWARE PSEUDOBULK PAIRED TESTS\n")
print(paired_results, row.names = FALSE)
cat("\n\nLOCKED 65-GENE COVERAGE\n")
cat("Present:", length(intersect(inflammasome_65, rownames(mnp))), "of 65\n")
cat("Missing:", paste(setdiff(inflammasome_65, rownames(mnp)), collapse = ", "), "\n")
cat("\n\nBROWN MODULE\n")
if (length(brown_genes) > 0L) {
  cat("Brown module source:", module_path, "\n")
  cat("Brown genes loaded:", length(brown_genes), "\n")
  cat(
    "Brown genes present in GSE210152 MNP expression matrix:",
    length(intersect(brown_genes, rownames(mnp))),
    "\n"
  )
} else {
  cat("Brown module analysis skipped because the module-assignment file was not found.\n")
}
cat("\n\nINTERPRETATION RULES\n")
cat("1. Treat the automatically selected cluster as a candidate PLIN2hi/TREM1hi iLAM-like cluster until its marker profile and donor representation are inspected.\n")
cat("2. Strong validation requires the candidate cluster to show PLIN2/TREM1 inflammatory-LAM markers across multiple donors, not only a single patient.\n")
cat("3. The primary external-validation question is whether GSE125771 Brown_WGCNA and/or Inflammasome_65 are enriched in the candidate inflammatory macrophage state.\n")
cat("4. Prefer donor-aware pseudobulk effect direction and cluster-marker enrichment over cell-level significance.\n")
cat("5. Do not interpret symptomatic-vs-asymptomatic differences as definitive because the dataset contains only 4 symptomatic and 2 asymptomatic donors.\n")
cat("6. GSE210152 contains CD45+ immune cells and cannot validate the blue structural/fibrous module.\n")
sink()

capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))

cat("\n============================================================\n")
cat("GSE210152 IDEA 3C COMPLETE\n")
cat("============================================================\n\n")
cat("Main outputs:\n")
cat("1. FINAL_GSE210152_MACROPHAGE_STATE_SUMMARY.txt\n")
cat("2. 07_MNP_cluster_state_scores.csv\n")
cat("3. 08_MNP_cluster_key_gene_average_expression.csv\n")
cat("4. 10_candidate_iLAM_donor_representation.csv\n")
cat("5. 12_gene_set_enrichment_by_MNP_cluster.csv\n")
cat("6. 15_pseudobulk_gene_set_scores.csv\n")
cat("7. 16_candidate_vs_other_within_donor_deltas.csv\n")
cat("8. 17_candidate_vs_other_paired_tests.csv\n")
cat("9. Figure_MNP_state_marker_dotplot.png\n")
cat("10. Figure_MNP_UMAP_clusters.png\n")
