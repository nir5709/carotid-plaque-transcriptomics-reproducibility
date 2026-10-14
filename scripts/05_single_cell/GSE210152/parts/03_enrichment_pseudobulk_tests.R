###############################################################################
# 18. LOAD GSE125771 BROWN MODULE
###############################################################################

module_path <- find_module_file(MODULE_FILE)
brown_genes <- character(0)

if (!is.na(module_path) && file.exists(module_path)) {
  brown_genes <- read_brown_module(module_path)
  message("Brown-module genes loaded: ", length(brown_genes))
} else {
  warning("Brown WGCNA module file not found. Brown-module analyses will be skipped.")
  if (REQUIRE_BROWN_MODULE) stop("Set MODULE_FILE to your Idea-2 module-assignment CSV.")
}

###############################################################################
# 19. GENE-SET ENRICHMENT AMONG MNP CLUSTER MARKERS
###############################################################################

mnp_counts <- get_counts(mnp)

universe <- rownames(mnp_counts)[
  Matrix::rowSums(mnp_counts > 0) >= max(10L, ceiling(0.01 * ncol(mnp)))
]
universe <- sub("\\.\\d+$", "", universe)

gene_sets_for_enrichment <- c(
  list(Inflammasome_65 = inflammasome_65),
  state_sets
)
if (length(brown_genes) > 0L) {
  gene_sets_for_enrichment <- c(
    list(Brown_WGCNA = brown_genes),
    gene_sets_for_enrichment
  )
}

clusters <- sort(unique(as.character(mnp$seurat_clusters)))
enrichment_rows <- list()
rr <- 1L

for (cl in clusters) {
  mk <- mnp_markers[
    as.character(mnp_markers$cluster) == cl &
      mnp_markers$p_val_adj < 0.05 &
      mnp_markers[[fc_col]] > 0.25,
    , drop = FALSE
  ]

  marker_genes <- unique(
    sub("\\.\\d+$", "", toupper(mk$gene))
  )

  for (set_name in names(gene_sets_for_enrichment)) {
    tst <- fisher_gene_set(
      markers = marker_genes,
      geneset = gene_sets_for_enrichment[[set_name]],
      universe = universe
    )

    enrichment_rows[[rr]] <- data.frame(
      Cluster = cl,
      GeneSet = set_name,
      N_cluster_markers = length(intersect(marker_genes, universe)),
      N_geneset_in_universe = length(intersect(gene_sets_for_enrichment[[set_name]], universe)),
      Overlap = unname(tst["overlap"]),
      Odds_ratio = unname(tst["odds_ratio"]),
      P_value = unname(tst["p_value"]),
      stringsAsFactors = FALSE
    )
    rr <- rr + 1L
  }
}

enrichment <- do.call(rbind, enrichment_rows)
enrichment$FDR_within_gene_set <- NA_real_

for (gs in unique(enrichment$GeneSet)) {
  ii <- which(enrichment$GeneSet == gs)
  enrichment$FDR_within_gene_set[ii] <- p.adjust(
    enrichment$P_value[ii],
    method = "BH"
  )
}

write.csv(
  enrichment,
  file.path(OUTDIR, "12_gene_set_enrichment_by_MNP_cluster.csv"),
  row.names = FALSE
)

###############################################################################
# 20. PSEUDOBULK BY PATIENT × MNP CLUSTER
###############################################################################

pb_group <- paste(mnp$patient, mnp$seurat_clusters, sep = "__CL")
pb_design <- Matrix::sparse.model.matrix(~ 0 + pb_group)
colnames(pb_design) <- sub("^pb_group", "", colnames(pb_design))
pb_counts <- mnp_counts %*% pb_design

pb_meta <- data.frame(
  Pseudobulk = colnames(pb_counts),
  stringsAsFactors = FALSE
)
pb_meta$Patient <- sub("__CL.*$", "", pb_meta$Pseudobulk)
pb_meta$Cluster <- sub("^.*__CL", "", pb_meta$Pseudobulk)

cell_counts_pb <- as.data.frame(
  table(Pseudobulk = pb_group),
  stringsAsFactors = FALSE
)
names(cell_counts_pb)[2] <- "N_cells"
pb_meta <- merge(pb_meta, cell_counts_pb, by = "Pseudobulk", all.x = TRUE, sort = FALSE)

keep_pb <- pb_meta$N_cells >= MIN_CELLS_PER_PSEUDOBULK
pb_counts <- pb_counts[, keep_pb, drop = FALSE]
pb_meta <- pb_meta[keep_pb, , drop = FALSE]
pb_meta <- pb_meta[match(colnames(pb_counts), pb_meta$Pseudobulk), , drop = FALSE]

###############################################################################
# 21. NORMALIZED PSEUDOBULK GENE-SET SCORES
###############################################################################

dge <- edgeR::DGEList(counts = pb_counts)
keep_genes_pb <- edgeR::filterByExpr(dge, group = pb_meta$Cluster)
dge <- dge[keep_genes_pb, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge)
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

pb_scores <- pb_meta
pb_scores$Inflammasome65_score <- score_pseudobulk_gene_set(logcpm, inflammasome_65)
pb_scores$iLAM_score <- score_pseudobulk_gene_set(logcpm, state_sets$PLIN2_TREM1_iLAM)
pb_scores$TREM2_LAM_score <- score_pseudobulk_gene_set(logcpm, state_sets$TREM2_LAM)

if (length(brown_genes) > 0L) {
  pb_scores$Brown_WGCNA_score <- score_pseudobulk_gene_set(logcpm, brown_genes)
}

pb_scores$Candidate_iLAM <- ifelse(
  as.character(pb_scores$Cluster) == candidate_cluster,
  "Candidate_PLIN2_TREM1_iLAM",
  "Other_MNP"
)

pb_scores$Symptom_status <- ifelse(
  pb_scores$Patient %in% c("P1","P2","P6","P7"),
  "Symptomatic",
  "Asymptomatic"
)

write.csv(
  pb_scores,
  file.path(OUTDIR, "15_pseudobulk_gene_set_scores.csv"),
  row.names = FALSE
)

###############################################################################
# 22. WITHIN-DONOR CANDIDATE-vs-OTHER STATE DELTAS
###############################################################################

score_columns <- c(
  "Inflammasome65_score",
  "iLAM_score",
  "TREM2_LAM_score"
)
if ("Brown_WGCNA_score" %in% names(pb_scores)) {
  score_columns <- c(score_columns, "Brown_WGCNA_score")
}

delta_rows <- list()
dd <- 1L

for (pat in unique(pb_scores$Patient)) {
  tmp <- pb_scores[pb_scores$Patient == pat, , drop = FALSE]
  cand <- tmp[
    tmp$Candidate_iLAM == "Candidate_PLIN2_TREM1_iLAM",
    , drop = FALSE
  ]
  other <- tmp[tmp$Candidate_iLAM == "Other_MNP", , drop = FALSE]
  if (nrow(cand) < 1L || nrow(other) < 1L) next

  for (sc in score_columns) {
    candidate_value <- mean(cand[[sc]], na.rm = TRUE)
    other_mean <- mean(other[[sc]], na.rm = TRUE)
    delta_rows[[dd]] <- data.frame(
      Patient = pat,
      Score = sc,
      Candidate_value = candidate_value,
      Other_MNP_mean = other_mean,
      Delta_candidate_minus_other = candidate_value - other_mean,
      stringsAsFactors = FALSE
    )
    dd <- dd + 1L
  }
}

if (length(delta_rows) > 0L) {
  deltas <- do.call(rbind, delta_rows)
} else {
  deltas <- data.frame(
    Patient = character(0),
    Score = character(0),
    Candidate_value = numeric(0),
    Other_MNP_mean = numeric(0),
    Delta_candidate_minus_other = numeric(0)
  )
}

write.csv(
  deltas,
  file.path(OUTDIR, "16_candidate_vs_other_within_donor_deltas.csv"),
  row.names = FALSE
)

###############################################################################
# 23. DONOR-LEVEL PAIRED TESTS
#
# The final manuscript used the repaired exact two-sided signed-rank results.
# Both exact and approximate P values are retained here; BH adjustment is across
# the tested pseudobulk scores.
###############################################################################

paired_results <- do.call(
  rbind,
  lapply(unique(deltas$Score), function(sc) {
    x <- deltas[deltas$Score == sc, , drop = FALSE]
    positive_n <- sum(x$Delta_candidate_minus_other > 0, na.rm = TRUE)
    negative_n <- sum(x$Delta_candidate_minus_other < 0, na.rm = TRUE)

    if (nrow(x) < 3L) {
      return(data.frame(
        Score = sc,
        N_donors = nrow(x),
        Positive_delta_donors = positive_n,
        Negative_delta_donors = negative_n,
        Median_delta = if (nrow(x) > 0L) median(x$Delta_candidate_minus_other, na.rm = TRUE) else NA_real_,
        Wilcoxon_V = NA_real_,
        P_value_exact = NA_real_,
        P_value_approx = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    tst_exact <- suppressWarnings(
      wilcox.test(
        x$Candidate_value,
        x$Other_MNP_mean,
        paired = TRUE,
        exact = TRUE,
        alternative = "two.sided"
      )
    )
    tst_approx <- suppressWarnings(
      wilcox.test(
        x$Candidate_value,
        x$Other_MNP_mean,
        paired = TRUE,
        exact = FALSE,
        alternative = "two.sided"
      )
    )

    data.frame(
      Score = sc,
      N_donors = nrow(x),
      Positive_delta_donors = positive_n,
      Negative_delta_donors = negative_n,
      Median_delta = median(x$Delta_candidate_minus_other, na.rm = TRUE),
      Wilcoxon_V = unname(tst_exact$statistic),
      P_value_exact = tst_exact$p.value,
      P_value_approx = tst_approx$p.value,
      stringsAsFactors = FALSE
    )
  })
)

if (nrow(paired_results) > 0L) {
  paired_results$FDR_exact <- p.adjust(paired_results$P_value_exact, method = "BH")
  paired_results$FDR_approx <- p.adjust(paired_results$P_value_approx, method = "BH")
}

write.csv(
  paired_results,
  file.path(OUTDIR, "17_candidate_vs_other_paired_tests.csv"),
  row.names = FALSE
)
