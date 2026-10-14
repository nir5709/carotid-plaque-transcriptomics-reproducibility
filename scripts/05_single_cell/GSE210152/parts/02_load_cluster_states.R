###############################################################################
# 5. DOWNLOAD / PREPARE GSE210152 H5AD
###############################################################################

if (!file.exists(H5AD) && !file.exists(H5AD_GZ)) {
  if (!AUTO_DOWNLOAD) stop("Input file not found: ", H5AD_GZ)
  message("Downloading GSE210152 H5AD from GEO...")
  download.file(GEO_URL, destfile = H5AD_GZ, mode = "wb", quiet = FALSE)
}

if (!file.exists(H5AD) && file.exists(H5AD_GZ)) {
  message("Decompressing H5AD...")
  R.utils::gunzip(
    H5AD_GZ,
    destname = H5AD,
    remove = FALSE,
    overwrite = TRUE
  )
}

if (!file.exists(H5AD)) stop("Could not locate or create: ", H5AD)

###############################################################################
# 6. READ H5AD
###############################################################################

message("Reading H5AD...")
sce <- zellkonverter::readH5AD(H5AD, use_hdf5 = FALSE)

capture.output({
  cat("Assays:\n")
  print(SummarizedExperiment::assayNames(sce))
  cat("\nCell metadata columns:\n")
  print(colnames(SummarizedExperiment::colData(sce)))
  cat("\nGene metadata columns:\n")
  print(colnames(SummarizedExperiment::rowData(sce)))
}, file = file.path(OUTDIR, "01_H5AD_structure.txt"))

assay_names <- SummarizedExperiment::assayNames(sce)
preferred_assays <- c("counts", "raw_counts", "X")
chosen_assay <- preferred_assays[preferred_assays %in% assay_names][1]
if (length(chosen_assay) == 0L || is.na(chosen_assay)) chosen_assay <- assay_names[1]
message("Using H5AD assay: ", chosen_assay)

counts <- SummarizedExperiment::assay(sce, chosen_assay)
if (!inherits(counts, "Matrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)

###############################################################################
# 7. GENE NAMES
###############################################################################

gene_names <- rownames(sce)

if (is.null(gene_names) || mean(grepl("^ENSG", gene_names)) > 0.5) {
  rd <- as.data.frame(SummarizedExperiment::rowData(sce))
  candidate_symbol_cols <- c(
    "gene_symbol", "gene_symbols", "symbol", "gene_name", "feature_name", "name"
  )
  hit <- candidate_symbol_cols[candidate_symbol_cols %in% colnames(rd)][1]
  if (length(hit) > 0L && !is.na(hit)) gene_names <- as.character(rd[[hit]])
}

if (is.null(gene_names) || length(gene_names) != nrow(counts)) {
  stop("Could not establish gene names from the H5AD.")
}

gene_names <- toupper(trimws(gene_names))
bad <- is.na(gene_names) | !nzchar(gene_names)
gene_names[bad] <- paste0("UNKNOWN_", which(bad))
rownames(counts) <- make.unique(gene_names)

###############################################################################
# 8. CELL METADATA + PATIENT DETECTION
###############################################################################

meta_raw <- as.data.frame(SummarizedExperiment::colData(sce))
keep_meta <- vapply(meta_raw, function(x) is.atomic(x) || is.factor(x), logical(1))
meta <- meta_raw[, keep_meta, drop = FALSE]
rownames(meta) <- colnames(sce)
if (is.null(colnames(counts))) colnames(counts) <- rownames(meta)

patient_info <- detect_patient_vector(meta, colnames(counts))
meta$patient <- patient_info$patient

write.csv(
  data.frame(metadata_column = colnames(meta_raw)),
  file.path(OUTDIR, "02_H5AD_metadata_columns.csv"),
  row.names = FALSE
)

if (length(unique(na.omit(meta$patient))) < 4L) {
  write.csv(
    head(meta, 100),
    file.path(OUTDIR, "STOP_patient_metadata_preview.csv"),
    row.names = TRUE
  )
  stop(
    "Could not reliably detect patient identities from H5AD metadata. ",
    "Inspect 01_H5AD_structure.txt and STOP_patient_metadata_preview.csv."
  )
}

meta$symptom_status <- ifelse(
  meta$patient %in% c("P1","P2","P6","P7"),
  "Symptomatic",
  ifelse(meta$patient %in% c("P3","P5"), "Asymptomatic", NA_character_)
)

###############################################################################
# 9. CREATE SEURAT OBJECT
###############################################################################

obj <- Seurat::CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  project = GSE,
  min.cells = 3,
  min.features = 0
)
rm(counts, sce, meta_raw, meta)
gc()

obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
qc_before <- ncol(obj)

if (APPLY_EXTRA_QC) {
  obj <- subset(
    obj,
    subset = nFeature_RNA >= MIN_FEATURES_EXTRA_QC &
      percent.mt <= MAX_PERCENT_MT_EXTRA_QC
  )
}
qc_after <- ncol(obj)

write.csv(
  data.frame(
    Cells_before_optional_QC = qc_before,
    Cells_after_optional_QC = qc_after,
    Extra_QC_applied = APPLY_EXTRA_QC,
    stringsAsFactors = FALSE
  ),
  file.path(OUTDIR, "03_QC_summary.csv"),
  row.names = FALSE
)

###############################################################################
# 10. NORMALIZE + DONOR-AWARE BROAD IMMUNE CLUSTERING
###############################################################################

obj <- Seurat::NormalizeData(obj, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = N_VARIABLE_FEATURES,
  verbose = FALSE
)
obj <- Seurat::ScaleData(
  obj,
  features = Seurat::VariableFeatures(obj),
  verbose = FALSE
)
obj <- Seurat::RunPCA(
  obj,
  features = Seurat::VariableFeatures(obj),
  npcs = N_PCS,
  verbose = FALSE
)

harmony_ok <- TRUE
obj <- tryCatch(
  harmony::RunHarmony(
    obj,
    group.by.vars = "patient",
    reduction.use = "pca",
    dims.use = DIMS_USE,
    verbose = FALSE
  ),
  error = function(e) {
    warning("Harmony failed; using PCA for clustering. Error: ", conditionMessage(e))
    harmony_ok <<- FALSE
    obj
  }
)

reduction_use <- if (harmony_ok && "harmony" %in% names(obj@reductions)) "harmony" else "pca"
obj <- Seurat::FindNeighbors(obj, reduction = reduction_use, dims = DIMS_USE, verbose = FALSE)
obj <- Seurat::FindClusters(obj, resolution = BROAD_RESOLUTION, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, reduction = reduction_use, dims = DIMS_USE, verbose = FALSE)

###############################################################################
# 11. BROAD CLUSTER MARKER-SCORE ANNOTATION
###############################################################################

avg_broad <- Seurat::AverageExpression(
  obj,
  assays = "RNA",
  group.by = "seurat_clusters",
  slot = "data",
  verbose = FALSE
)$RNA

broad_scores <- cluster_gene_set_scores(avg_broad, broad_sets)
write.csv(
  broad_scores,
  file.path(OUTDIR, "04_broad_cluster_marker_scores.csv"),
  row.names = FALSE
)

broad_wide <- reshape(
  broad_scores[, c("Cluster","GeneSet","Score")],
  idvar = "Cluster",
  timevar = "GeneSet",
  direction = "wide"
)
score_cols <- grep("^Score\\.", names(broad_wide), value = TRUE)
broad_wide$Broad_label <- apply(
  broad_wide[, score_cols, drop = FALSE],
  1,
  function(x) sub("^Score\\.", "", score_cols[which.max(x)])
)
broad_cluster_key <- sub(
  "^g(?=[0-9]+$)",
  "",
  as.character(broad_wide$Cluster),
  perl = TRUE
)
cluster_to_label <- setNames(broad_wide$Broad_label, broad_cluster_key)

broad_labels <- unname(
  cluster_to_label[as.character(obj$seurat_clusters)]
)

if (anyNA(broad_labels)) {
  stop(
    "Broad-cluster annotation failed for cluster(s): ",
    paste(
      unique(as.character(obj$seurat_clusters)[is.na(broad_labels)]),
      collapse = ", "
    )
  )
}

obj$Broad_label <- broad_labels

write.csv(
  broad_wide,
  file.path(OUTDIR, "05_broad_cluster_assignments.csv"),
  row.names = FALSE
)

###############################################################################
# 12. BROAD FIGURES
###############################################################################

p_broad <- Seurat::DimPlot(
  obj,
  reduction = "umap",
  group.by = "Broad_label",
  label = TRUE,
  repel = TRUE
) + ggtitle("GSE210152 broad immune-cell annotation")

save_plot_both(
  p_broad,
  file.path(OUTDIR, "Figure_broad_immune_UMAP"),
  width = 9,
  height = 7
)

broad_marker_vec <- intersect(unique(unlist(broad_sets)), rownames(obj))
p_broad_dot <- Seurat::DotPlot(
  obj,
  features = broad_marker_vec,
  group.by = "Broad_label",
  dot.scale = 5
) + RotatedAxis() + ggtitle("Broad immune-cell marker expression")

save_plot_both(
  p_broad_dot,
  file.path(OUTDIR, "Figure_broad_immune_marker_dotplot"),
  width = 18,
  height = 7
)

###############################################################################
# 13. SUBSET MONONUCLEAR PHAGOCYTES
###############################################################################

mnp_cells <- colnames(obj)[obj$Broad_label %in% c("Macrophage", "Monocyte")]
if (length(mnp_cells) < 200L) {
  stop("Fewer than 200 Macrophage/Monocyte cells were identified. Inspect broad annotation first.")
}
mnp <- subset(obj, cells = mnp_cells)
message("MNP cells selected: ", ncol(mnp))

###############################################################################
# 14. RECLUSTER MNP COMPARTMENT
###############################################################################

mnp <- Seurat::NormalizeData(mnp, verbose = FALSE)
mnp <- Seurat::FindVariableFeatures(
  mnp,
  selection.method = "vst",
  nfeatures = MNP_N_VARIABLE_FEATURES,
  verbose = FALSE
)
mnp <- Seurat::ScaleData(
  mnp,
  features = Seurat::VariableFeatures(mnp),
  verbose = FALSE
)
mnp <- Seurat::RunPCA(
  mnp,
  features = Seurat::VariableFeatures(mnp),
  npcs = MNP_N_PCS,
  verbose = FALSE
)

mnp_harmony_ok <- TRUE
mnp <- tryCatch(
  harmony::RunHarmony(
    mnp,
    group.by.vars = "patient",
    reduction.use = "pca",
    dims.use = MNP_DIMS_USE,
    verbose = FALSE
  ),
  error = function(e) {
    warning("MNP Harmony failed; using PCA. Error: ", conditionMessage(e))
    mnp_harmony_ok <<- FALSE
    mnp
  }
)

mnp_reduction <- if (
  mnp_harmony_ok && "harmony" %in% names(mnp@reductions)
) "harmony" else "pca"

mnp <- Seurat::FindNeighbors(mnp, reduction = mnp_reduction, dims = MNP_DIMS_USE, verbose = FALSE)
mnp <- Seurat::FindClusters(mnp, resolution = MNP_RESOLUTION, verbose = FALSE)
mnp <- Seurat::RunUMAP(mnp, reduction = mnp_reduction, dims = MNP_DIMS_USE, verbose = FALSE)

###############################################################################
# 15. MNP STATE SCORES + CANDIDATE iLAM-LIKE CLUSTER
###############################################################################

avg_mnp <- Seurat::AverageExpression(
  mnp,
  assays = "RNA",
  group.by = "seurat_clusters",
  slot = "data",
  verbose = FALSE
)$RNA

state_scores <- cluster_gene_set_scores(avg_mnp, state_sets)
write.csv(
  state_scores,
  file.path(OUTDIR, "07_MNP_cluster_state_scores.csv"),
  row.names = FALSE
)

key_present <- intersect(key_genes, rownames(avg_mnp))
key_avg <- as.data.frame(t(as.matrix(avg_mnp[key_present, , drop = FALSE])))
key_avg$Cluster <- rownames(key_avg)
write.csv(
  key_avg,
  file.path(OUTDIR, "08_MNP_cluster_key_gene_average_expression.csv"),
  row.names = FALSE
)

ilam_rows <- state_scores[state_scores$GeneSet == "PLIN2_TREM1_iLAM", , drop = FALSE]
candidate_cluster_display <- as.character(ilam_rows$Cluster[which.max(ilam_rows$Score)])
# Integrated repair of the Seurat-v5 g5-versus-5 mismatch found in the final run.
candidate_cluster <- sub("^g(?=[0-9]+$)", "", candidate_cluster_display, perl = TRUE)

message(
  "Candidate PLIN2/TREM1 inflammatory-LAM-like cluster: ",
  candidate_cluster_display,
  " -> normalized ",
  candidate_cluster
)

mnp$Candidate_iLAM <- ifelse(
  as.character(mnp$seurat_clusters) == candidate_cluster,
  "Candidate_PLIN2_TREM1_iLAM",
  "Other_MNP"
)

###############################################################################
# 16. DONOR REPRESENTATION OF MNP CLUSTERS
###############################################################################

donor_cluster <- as.data.frame(
  table(Patient = mnp$patient, Cluster = mnp$seurat_clusters),
  stringsAsFactors = FALSE
)
names(donor_cluster)[3] <- "Cells"
patient_total <- aggregate(Cells ~ Patient, donor_cluster, sum)
names(patient_total)[2] <- "Total_MNP_cells"
donor_cluster <- merge(donor_cluster, patient_total, by = "Patient", all.x = TRUE)
donor_cluster$Fraction_of_patient_MNP <- donor_cluster$Cells / donor_cluster$Total_MNP_cells

write.csv(
  donor_cluster,
  file.path(OUTDIR, "09_MNP_cluster_donor_representation.csv"),
  row.names = FALSE
)

candidate_donor <- donor_cluster[
  as.character(donor_cluster$Cluster) == candidate_cluster,
  , drop = FALSE
]
candidate_donor$Symptom_status <- ifelse(
  candidate_donor$Patient %in% c("P1","P2","P6","P7"),
  "Symptomatic",
  "Asymptomatic"
)
write.csv(
  candidate_donor,
  file.path(OUTDIR, "10_candidate_iLAM_donor_representation.csv"),
  row.names = FALSE
)

###############################################################################
# 17. MNP MARKERS
###############################################################################

Seurat::Idents(mnp) <- "seurat_clusters"
mnp_markers <- Seurat::FindAllMarkers(
  mnp,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)
write.csv(
  mnp_markers,
  file.path(OUTDIR, "11_MNP_all_positive_cluster_markers.csv"),
  row.names = FALSE
)
fc_col <- if ("avg_log2FC" %in% colnames(mnp_markers)) "avg_log2FC" else "avg_logFC"