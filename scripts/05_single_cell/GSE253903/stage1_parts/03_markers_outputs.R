###############################################################################
# 11. CLUSTER MARKER GENES
###############################################################################

message("Finding positive markers for each cluster...")
Seurat::Idents(merged) <- "seurat_clusters"

markers_all <- Seurat::FindAllMarkers(
  merged,
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

write.csv(
  markers_all,
  file.path(OUTDIR, "06_all_positive_markers_per_cluster.csv"),
  row.names = FALSE
)

fc_col <- if ("avg_log2FC" %in% colnames(markers_all)) {
  "avg_log2FC"
} else if ("avg_logFC" %in% colnames(markers_all)) {
  "avg_logFC"
} else {
  stop("Could not identify the marker fold-change column.")
}

markers_all <- markers_all[
  order(markers_all$cluster, -markers_all[[fc_col]]),
  ,
  drop = FALSE
]

top20 <- do.call(
  rbind,
  lapply(
    split(markers_all, markers_all$cluster),
    function(x) {
      head(
        x[
          order(x$p_val_adj, -x[[fc_col]]),
          ,
          drop = FALSE
        ],
        20
      )
    }
  )
)

write.csv(
  top20,
  file.path(OUTDIR, "06_top20_markers_per_cluster.csv"),
  row.names = FALSE
)

###############################################################################
# 12. CLUSTER DISTRIBUTION BY SUBJECT
###############################################################################

cluster_sample_table <- as.data.frame(
  table(
    Cluster = merged$seurat_clusters,
    Sample = merged$sample_id
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cluster_sample_table,
  file.path(OUTDIR, "05_cluster_cell_counts_by_sample.csv"),
  row.names = FALSE
)

cluster_total <- as.data.frame(
  table(Cluster = merged$seurat_clusters),
  stringsAsFactors = FALSE
)

write.csv(
  cluster_total,
  file.path(OUTDIR, "04_cluster_cell_counts_total.csv"),
  row.names = FALSE
)

###############################################################################
# 13. SAVE STAGE-1 OBJECT
###############################################################################

saveRDS(
  merged,
  file.path(OUTDIR, "GSE253903_STAGE1_clustered_unannotated.rds"),
  compress = FALSE
)

###############################################################################
# 14. FINAL STAGE-1 SUMMARY
###############################################################################

summary_file <- file.path(OUTDIR, "FINAL_STAGE1_SUMMARY.txt")
sink(summary_file)
cat("GSE253903 STAGE 1: QC AND CLUSTERING\n")
cat("===================================\n\n")
cat("GEO accession:", GSE, "\n")
cat("Raw scRNA samples detected:", length(matrix_files), "\n")
cat("Total cells retained after per-subject QC:", ncol(merged), "\n")
cat("Number of clusters:", length(unique(merged$seurat_clusters)), "\n\n")
cat("QC THRESHOLDS\n")
cat("Minimum detected genes:", MIN_FEATURES, "\n")
cat("Maximum detected genes:", MAX_FEATURES, "\n")
cat("Maximum UMI counts:", MAX_COUNTS, "\n")
cat("Maximum mitochondrial percentage:", MAX_PERCENT_MT, "\n\n")
cat("SAMPLE QC\n")
print(qc_summary, row.names = FALSE)
cat("\n\nCLUSTER TOTALS\n")
print(cluster_total, row.names = FALSE)
cat("\n\nNEXT STEP\n")
cat("Do NOT assign the final cell-type labels from cluster numbers alone.\n")
cat("Review the canonical-marker dotplot and top20 marker table.\n")
cat(
  "Lock broad labels (Macrophage/Myeloid, T, NK, B, Plasma, DC, ",
  "Neutrophil, Mast, Endothelial, Fibroblast, SMC) from marker evidence.\n"
)
cat(
  "Only after those labels are locked should Stage 2 test the discovery ",
  "Brown and Blue WGCNA modules and the locked inflammasome signature.\n"
)
sink()

capture.output(
  sessionInfo(),
  file = file.path(OUTDIR, "sessionInfo.txt")
)

cat("\n")
cat("============================================================\n")
cat("GSE253903 STAGE 1 COMPLETE\n")
cat("============================================================\n\n")
cat("Main review outputs:\n")
cat("1. FINAL_STAGE1_SUMMARY.txt\n")
cat("2. 03_sample_QC_summary.csv\n")
cat("3. 04_cluster_cell_counts_total.csv\n")
cat("4. 05_cluster_cell_counts_by_sample.csv\n")
cat("5. 06_top20_markers_per_cluster.csv\n")
cat("6. Figure_cluster_canonical_marker_dotplot.png\n")
cat("7. Figure_UMAP_clusters.png\n")
