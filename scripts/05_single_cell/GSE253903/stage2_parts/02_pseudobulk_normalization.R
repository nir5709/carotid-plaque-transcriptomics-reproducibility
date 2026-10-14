###############################################################################
# 5. LOCATE INPUT FILES
###############################################################################

rds_path <- find_one_file(
  RDS_FILE,
  "GSE253903_STAGE1_clustered_unannotated\\.rds$"
)

module_path <- find_one_file(
  MODULE_FILE,
  "Whole_transcriptome_module_assignments\\.csv$"
)

if (is.na(rds_path) || !file.exists(rds_path)) {
  stop("Could not find GSE253903_STAGE1_clustered_unannotated.rds")
}
if (is.na(module_path) || !file.exists(module_path)) {
  stop("Could not find Whole_transcriptome_module_assignments.csv")
}

message("Stage-1 object: ", rds_path)
message("Module assignment file: ", module_path)

###############################################################################
# 6. LOAD STAGE-1 OBJECT AND LOCK ANNOTATIONS
###############################################################################

obj <- readRDS(rds_path)
if (!"seurat_clusters" %in% colnames(obj[[]])) {
  stop("seurat_clusters missing from Stage-1 Seurat object.")
}

cluster_chr <- as.character(obj$seurat_clusters)
cluster_map <- setNames(locked_annotation$Broad_label, as.character(locked_annotation$Cluster))
detail_map <- setNames(locked_annotation$Detailed_label, as.character(locked_annotation$Cluster))
broad <- unname(cluster_map[cluster_chr])
detail <- unname(detail_map[cluster_chr])

if (anyNA(broad)) {
  stop("Unmapped cluster(s): ", paste(unique(cluster_chr[is.na(broad)]), collapse = ", "))
}

obj$Broad_label_locked <- broad
obj$Detailed_label_locked <- detail

annotation_counts <- as.data.frame(
  table(Cluster = cluster_chr, Broad_label = broad),
  stringsAsFactors = FALSE
)
annotation_counts <- annotation_counts[annotation_counts$Freq > 0, , drop = FALSE]
write.csv(annotation_counts, file.path(OUTDIR, "01_locked_annotation_cell_counts.csv"), row.names = FALSE)

###############################################################################
# 7. LOAD BROWN / BLUE MODULES
###############################################################################

mods <- read_modules(module_path)
brown_genes <- mods$brown
blue_genes <- mods$blue
message("Brown genes loaded: ", length(brown_genes))
message("Blue genes loaded: ", length(blue_genes))

###############################################################################
# 8. REMOVE AMBIGUOUS / MIXED CLUSTERS FROM PRIMARY VALIDATION
###############################################################################

obj_primary <- subset(
  obj,
  cells = colnames(obj)[!obj$Broad_label_locked %in% EXCLUDED_LABELS]
)
message("Cells in primary analysis: ", ncol(obj_primary))

###############################################################################
# 9. IDENTIFY DONOR/SAMPLE COLUMN
###############################################################################

meta <- obj_primary[[]]
candidate_sample_cols <- c("sample_id", "orig.ident", "Sample", "sample")
sample_col <- candidate_sample_cols[candidate_sample_cols %in% colnames(meta)][1]
if (length(sample_col) == 0L || is.na(sample_col)) stop("Could not identify donor/sample metadata column.")
obj_primary$Stage2_sample <- as.character(meta[[sample_col]])
message("Using donor/sample column: ", sample_col)

###############################################################################
# 10. DONOR × BROAD-CELL-TYPE PSEUDOBULK
###############################################################################

counts <- get_counts(obj_primary)
rownames(counts) <- toupper(sub("\\.\\d+$", "", rownames(counts)))

if (anyDuplicated(rownames(counts))) {
  message("Collapsing duplicated gene symbols using sparse matrix multiplication...")
  gene_group <- factor(rownames(counts), levels = unique(rownames(counts)))
  collapse_mat <- Matrix::sparseMatrix(
    i = as.integer(gene_group),
    j = seq_along(gene_group),
    x = 1,
    dims = c(nlevels(gene_group), length(gene_group)),
    dimnames = list(levels(gene_group), NULL)
  )
  counts <- collapse_mat %*% counts
  rownames(counts) <- levels(gene_group)
  rm(collapse_mat, gene_group)
  gc()
  message("Duplicate collapse complete. Matrix remains sparse: ", inherits(counts, "sparseMatrix"))
}

pb_group <- paste(obj_primary$Stage2_sample, obj_primary$Broad_label_locked, sep = "__")
design_pb <- Matrix::sparse.model.matrix(~ 0 + pb_group)
colnames(design_pb) <- sub("^pb_group", "", colnames(design_pb))
pb_counts <- counts %*% design_pb

cell_counts <- as.data.frame(table(Pseudobulk = pb_group), stringsAsFactors = FALSE)
names(cell_counts)[2] <- "N_cells"
pb_meta <- data.frame(Pseudobulk = colnames(pb_counts), stringsAsFactors = FALSE)
split_index <- regexpr("__", pb_meta$Pseudobulk, fixed = TRUE)
pb_meta$Sample <- substr(pb_meta$Pseudobulk, 1, split_index - 1)
pb_meta$Broad_label <- substr(pb_meta$Pseudobulk, split_index + 2, nchar(pb_meta$Pseudobulk))
pb_meta <- merge(pb_meta, cell_counts, by = "Pseudobulk", all.x = TRUE, sort = FALSE)
pb_meta <- pb_meta[match(colnames(pb_counts), pb_meta$Pseudobulk), , drop = FALSE]

keep_pb <- pb_meta$N_cells >= MIN_CELLS_PER_PSEUDOBULK
pb_counts <- pb_counts[, keep_pb, drop = FALSE]
pb_meta <- pb_meta[keep_pb, , drop = FALSE]
pb_meta <- pb_meta[match(colnames(pb_counts), pb_meta$Pseudobulk), , drop = FALSE]
write.csv(pb_meta, file.path(OUTDIR, "02_pseudobulk_metadata.csv"), row.names = FALSE)

###############################################################################
# 11. edgeR NORMALIZATION
###############################################################################

dge <- edgeR::DGEList(counts = pb_counts)
keep_gene <- edgeR::filterByExpr(dge, group = pb_meta$Broad_label)
dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge)
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
