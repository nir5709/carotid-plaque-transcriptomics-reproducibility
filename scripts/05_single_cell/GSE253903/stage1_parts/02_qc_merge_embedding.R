###############################################################################
# 6. READ EACH SUBJECT AND APPLY GEO QC
###############################################################################

object_list <- list()
qc_list <- list()

for (i in seq_along(matrix_files)) {
  matrix_file <- matrix_files[i]
  gsm <- gsm_from_file[i]
  message("")
  message("[", i, "/", length(matrix_files), "] Reading ", gsm)
  counts <- read_one_10x_triplet(matrix_file)
  colnames(counts) <- paste0(gsm, "_", colnames(counts))
  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = gsm,
    min.cells = 3,
    min.features = MIN_FEATURES
  )
  obj$sample_id <- gsm
  meta_row <- sample_metadata[sample_metadata$GSM == gsm, , drop = FALSE]
  obj$age <- if (nrow(meta_row) == 1L) meta_row$Age else NA_character_
  obj$sex <- if (nrow(meta_row) == 1L) meta_row$Sex else NA_character_
  obj$symptom <- if (nrow(meta_row) == 1L) meta_row$Symptom else NA_character_
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  before_n <- ncol(obj)
  qc_before <- obj[[]]
  obj <- subset(
    obj,
    subset =
      nFeature_RNA >= MIN_FEATURES &
      nFeature_RNA <= MAX_FEATURES &
      nCount_RNA <= MAX_COUNTS &
      percent.mt <= MAX_PERCENT_MT
  )
  after_n <- ncol(obj)
  qc_list[[gsm]] <- data.frame(
    Sample = gsm,
    Cells_before_final_QC = before_n,
    Cells_after_final_QC = after_n,
    Fraction_retained = after_n / before_n,
    Median_nFeature_after = median(obj$nFeature_RNA),
    Median_nCount_after = median(obj$nCount_RNA),
    Median_percent_mt_after = median(obj$percent.mt),
    stringsAsFactors = FALSE
  )
  object_list[[gsm]] <- obj
  rm(counts, obj, qc_before)
  gc()
}

qc_summary <- do.call(rbind, qc_list)
write.csv(qc_summary, file.path(OUTDIR, "03_sample_QC_summary.csv"), row.names = FALSE)

###############################################################################
# 7. MERGE THE 12 SUBJECTS
###############################################################################

message("")
message("Merging subjects...")
merged <- object_list[[1]]
if (length(object_list) > 1L) {
  for (i in 2:length(object_list)) {
    merged <- merge(merged, y = object_list[[i]])
  }
}
rm(object_list)
gc()

# Seurat v5 may create separate RNA layers after merge.
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  merged <- SeuratObject::JoinLayers(merged, assay = "RNA")
}

message("Merged cells after QC: ", ncol(merged))

###############################################################################
# 8. NORMALIZATION AND DONOR-AWARE EMBEDDING
###############################################################################

merged <- Seurat::NormalizeData(
  merged,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

merged <- Seurat::FindVariableFeatures(
  merged,
  selection.method = "vst",
  nfeatures = N_VARIABLE_FEATURES,
  verbose = FALSE
)

merged <- Seurat::ScaleData(
  merged,
  features = Seurat::VariableFeatures(merged),
  verbose = FALSE
)

merged <- Seurat::RunPCA(
  merged,
  features = Seurat::VariableFeatures(merged),
  npcs = N_PCS,
  verbose = FALSE
)

message("Running Harmony using subject/sample as the batch variable...")
merged <- harmony::RunHarmony(
  merged,
  group.by.vars = "sample_id",
  dims.use = DIMS_USE,
  verbose = FALSE
)

merged <- Seurat::FindNeighbors(
  merged,
  reduction = "harmony",
  dims = DIMS_USE,
  verbose = FALSE
)

merged <- Seurat::FindClusters(
  merged,
  resolution = CLUSTER_RESOLUTION,
  verbose = FALSE
)

merged <- Seurat::RunUMAP(
  merged,
  reduction = "harmony",
  dims = DIMS_USE,
  verbose = FALSE
)

###############################################################################
# 9. UMAP OUTPUTS
###############################################################################

p_clusters <- Seurat::DimPlot(
  merged,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("GSE253903: broad clustering after donor-aware embedding")

save_plot_both(
  p_clusters,
  file.path(OUTDIR, "Figure_UMAP_clusters"),
  width = 9,
  height = 7
)

p_sample <- Seurat::DimPlot(
  merged,
  reduction = "umap",
  group.by = "sample_id"
) + ggtitle("GSE253903: cells by subject")

save_plot_both(
  p_sample,
  file.path(OUTDIR, "Figure_UMAP_subjects"),
  width = 10,
  height = 7
)

if (any(!is.na(merged$symptom))) {
  p_symptom <- Seurat::DimPlot(
    merged,
    reduction = "umap",
    group.by = "symptom"
  ) + ggtitle("GSE253903: symptomatic vs asymptomatic subjects")
  save_plot_both(
    p_symptom,
    file.path(OUTDIR, "Figure_UMAP_symptom_status"),
    width = 9,
    height = 7
  )
}

###############################################################################
# 10. CANONICAL BROAD-LINEAGE MARKER DOTPLOT
###############################################################################

canonical_markers <- c(
  "LYZ","LST1","FCER1G","TYROBP","CTSD","CD68","C1QA","C1QB","C1QC",
  "CD3D","CD3E","TRBC1","IL7R","CD8A",
  "NKG7","GNLY","KLRD1",
  "MS4A1","CD79A","CD37","CD74",
  "MZB1","JCHAIN","SDC1",
  "FCER1A","CD1C","CLEC10A",
  "FCGR3B","CSF3R","S100A8","S100A9",
  "KIT","TPSAB1","TPSB2",
  "PECAM1","VWF","KDR","EMCN",
  "COL1A1","COL1A2","DCN","LUM","LOXL1",
  "ACTA2","TAGLN","MYH11","CNN1","LMOD1"
)

canonical_markers <- canonical_markers[canonical_markers %in% rownames(merged)]

p_dot <- Seurat::DotPlot(
  merged,
  features = canonical_markers,
  group.by = "seurat_clusters",
  dot.scale = 5
) + RotatedAxis() + ggtitle("Canonical markers for broad carotid-plaque cell types")

save_plot_both(
  p_dot,
  file.path(OUTDIR, "Figure_cluster_canonical_marker_dotplot"),
  width = 18,
  height = 8
)
