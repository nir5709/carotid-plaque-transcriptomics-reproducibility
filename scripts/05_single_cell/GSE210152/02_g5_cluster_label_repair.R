###############################################################################
# GSE210152 — FIX FOR g5 vs 5 CLUSTER-NAME MISMATCH
#
# The upstream state-score table uses Seurat-v5 labels such as g0...g5,
# whereas downstream tables use numeric labels 0...5. This repair normalizes
# the candidate cluster label and regenerates the candidate-specific outputs.
# It does NOT rerun clustering.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

find_one <- function(filename) {
  if (file.exists(filename)) {
    return(normalizePath(filename, winslash = "/", mustWork = TRUE))
  }
  hits <- list.files(
    ".",
    pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )
  if (length(hits) == 0L) return(NA_character_)
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

strip_g <- function(x) {
  sub("^g(?=[0-9]+$)", "", as.character(x), perl = TRUE)
}

STATE_FILE <- find_one("07_MNP_cluster_state_scores.csv")
KEY_FILE <- find_one("08_MNP_cluster_key_gene_average_expression.csv")
DONOR_FILE <- find_one("09_MNP_cluster_donor_representation.csv")
ENRICH_FILE <- find_one("12_gene_set_enrichment_by_MNP_cluster.csv")
PB_FILE <- find_one("15_pseudobulk_gene_set_scores.csv")
RDS_FILE <- find_one("GSE210152_MNP_clustered_scored.rds")

required <- c(STATE_FILE, KEY_FILE, ENRICH_FILE, PB_FILE)
if (any(is.na(required))) {
  stop(
    "Could not find all required GSE210152 output files.\n",
    "Required:\n",
    "07_MNP_cluster_state_scores.csv\n",
    "08_MNP_cluster_key_gene_average_expression.csv\n",
    "12_gene_set_enrichment_by_MNP_cluster.csv\n",
    "15_pseudobulk_gene_set_scores.csv"
  )
}

OUTDIR <- dirname(STATE_FILE)
message("Using result directory: ", OUTDIR)

state_scores <- read.csv(STATE_FILE, check.names = FALSE, stringsAsFactors = FALSE)
key_avg <- read.csv(KEY_FILE, check.names = FALSE, stringsAsFactors = FALSE)
enrichment <- read.csv(ENRICH_FILE, check.names = FALSE, stringsAsFactors = FALSE)
pb_scores <- read.csv(PB_FILE, check.names = FALSE, stringsAsFactors = FALSE)

ilam_rows <- state_scores[
  state_scores$GeneSet == "PLIN2_TREM1_iLAM",
  , drop = FALSE
]
if (nrow(ilam_rows) == 0L) {
  stop("PLIN2_TREM1_iLAM rows were not found in the state-score table.")
}

candidate_cluster_display <- as.character(
  ilam_rows$Cluster[which.max(ilam_rows$Score)]
)
candidate_cluster_numeric <- strip_g(candidate_cluster_display)

message("Candidate cluster in state-score table: ", candidate_cluster_display)
message("Normalized cluster ID for downstream matching: ", candidate_cluster_numeric)

candidate_state_scores <- state_scores[
  strip_g(state_scores$Cluster) == candidate_cluster_numeric,
  , drop = FALSE
]
candidate_key_avg <- key_avg[
  strip_g(key_avg$Cluster) == candidate_cluster_numeric,
  , drop = FALSE
]
candidate_enrichment <- enrichment[
  strip_g(enrichment$Cluster) == candidate_cluster_numeric,
  , drop = FALSE
]

write.csv(
  candidate_state_scores,
  file.path(OUTDIR, "18_candidate_iLAM_state_scores_FIXED.csv"),
  row.names = FALSE
)
write.csv(
  candidate_enrichment,
  file.path(OUTDIR, "19_candidate_iLAM_gene_set_enrichment_FIXED.csv"),
  row.names = FALSE
)
write.csv(
  candidate_key_avg,
  file.path(OUTDIR, "20_candidate_iLAM_key_gene_expression_FIXED.csv"),
  row.names = FALSE
)

candidate_donor <- NULL

if (!is.na(DONOR_FILE) && file.exists(DONOR_FILE)) {
  donor_cluster <- read.csv(DONOR_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  candidate_donor <- donor_cluster[
    strip_g(donor_cluster$Cluster) == candidate_cluster_numeric,
    , drop = FALSE
  ]
} else if (!is.na(RDS_FILE) && file.exists(RDS_FILE)) {
  message("09_MNP_cluster_donor_representation.csv not found; rebuilding from RDS.")
  mnp <- readRDS(RDS_FILE)
  if (!"patient" %in% colnames(mnp[[]])) {
    stop("patient metadata is missing from the saved MNP object.")
  }
  cluster_values <- as.character(mnp$seurat_clusters)
  donor_cluster <- as.data.frame(
    table(Patient = mnp$patient, Cluster = cluster_values),
    stringsAsFactors = FALSE
  )
  names(donor_cluster)[3] <- "Cells"
  patient_total <- aggregate(Cells ~ Patient, donor_cluster, sum)
  names(patient_total)[2] <- "Total_MNP_cells"
  donor_cluster <- merge(donor_cluster, patient_total, by = "Patient", all.x = TRUE)
  donor_cluster$Fraction_of_patient_MNP <- donor_cluster$Cells / donor_cluster$Total_MNP_cells
  candidate_donor <- donor_cluster[
    strip_g(donor_cluster$Cluster) == candidate_cluster_numeric,
    , drop = FALSE
  ]
} else {
  warning(
    "Could not rebuild exact donor representation because neither ",
    "09_MNP_cluster_donor_representation.csv nor the saved RDS was found."
  )
  candidate_donor <- data.frame()
}

if (nrow(candidate_donor) > 0L) {
  candidate_donor$Symptom_status <- ifelse(
    candidate_donor$Patient %in% c("P1", "P2", "P6", "P7"),
    "Symptomatic",
    ifelse(
      candidate_donor$Patient %in% c("P3", "P5"),
      "Asymptomatic",
      NA_character_
    )
  )
}

write.csv(
  candidate_donor,
  file.path(OUTDIR, "10_candidate_iLAM_donor_representation_FIXED.csv"),
  row.names = FALSE
)

pb_scores$Cluster_norm <- strip_g(pb_scores$Cluster)
pb_scores$Candidate_iLAM <- ifelse(
  pb_scores$Cluster_norm == candidate_cluster_numeric,
  "Candidate_PLIN2_TREM1_iLAM",
  "Other_MNP"
)

write.csv(
  pb_scores,
  file.path(OUTDIR, "15_pseudobulk_gene_set_scores_FIXED.csv"),
  row.names = FALSE
)

score_columns <- intersect(
  c(
    "Inflammasome65_score",
    "iLAM_score",
    "TREM2_LAM_score",
    "Brown_WGCNA_score"
  ),
  names(pb_scores)
)

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
  file.path(OUTDIR, "16_candidate_vs_other_within_donor_deltas_FIXED.csv"),
  row.names = FALSE
)

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
        Median_delta = if (nrow(x) > 0L) {
          median(x$Delta_candidate_minus_other, na.rm = TRUE)
        } else {
          NA_real_
        },
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
  paired_results$FDR_exact <- p.adjust(
    paired_results$P_value_exact,
    method = "BH"
  )
  paired_results$FDR_approx <- p.adjust(
    paired_results$P_value_approx,
    method = "BH"
  )
}

write.csv(
  paired_results,
  file.path(OUTDIR, "17_candidate_vs_other_paired_tests_FIXED.csv"),
  row.names = FALSE
)

summary_file <- file.path(
  OUTDIR,
  "FINAL_GSE210152_MACROPHAGE_STATE_SUMMARY_FIXED.txt"
)

sink(summary_file)
cat("GSE210152 — MACROPHAGE-STATE VALIDATION (FIXED)\n")
cat("================================================\n\n")
cat("Original Seurat-v5 candidate label:", candidate_cluster_display, "\n")
cat("Normalized downstream cluster ID:", candidate_cluster_numeric, "\n\n")
cat("CANDIDATE iLAM-LIKE STATE SCORES\n")
print(candidate_state_scores, row.names = FALSE)
cat("\n\nCANDIDATE iLAM-LIKE KEY GENE EXPRESSION\n")
print(candidate_key_avg, row.names = FALSE)
cat("\n\nCANDIDATE iLAM-LIKE GENE-SET ENRICHMENT\n")
print(candidate_enrichment, row.names = FALSE)
cat("\n\nCANDIDATE iLAM-LIKE DONOR REPRESENTATION\n")
print(candidate_donor, row.names = FALSE)
cat("\n\nWITHIN-DONOR CANDIDATE-vs-OTHER DELTAS\n")
print(deltas, row.names = FALSE)
cat("\n\nDONOR-LEVEL PAIRED TESTS\n")
print(paired_results, row.names = FALSE)
cat("\n\nINTERPRETATION\n")
cat("1. The previous empty candidate-specific tables were caused by the label mismatch g5 versus 5, not by absence of the biological state.\n")
cat("2. Treat the cluster as a PLIN2/TREM1 inflammatory-LAM-like state only if marker expression, enrichment, and donor representation all support it.\n")
cat("3. With only six total donors, prioritize direction consistency and effect size; exact Wilcoxon tests have limited attainable significance.\n")
cat("4. The Brown WGCNA module need not be exclusive to this macrophage state. Evidence of enrichment in the candidate state supports overlap, not one-to-one identity.\n")
cat("5. Inflammasome-65 enrichment and Brown-module enrichment are separate questions and should not be conflated.\n")
sink()

cat("\n============================================================\n")
cat("GSE210152 g5/5 FIX COMPLETE\n")
cat("============================================================\n\n")
cat("Candidate cluster:", candidate_cluster_display, "->", candidate_cluster_numeric, "\n\n")
cat("Main fixed outputs:\n")
cat("1. FINAL_GSE210152_MACROPHAGE_STATE_SUMMARY_FIXED.txt\n")
cat("2. 10_candidate_iLAM_donor_representation_FIXED.csv\n")
cat("3. 16_candidate_vs_other_within_donor_deltas_FIXED.csv\n")
cat("4. 17_candidate_vs_other_paired_tests_FIXED.csv\n")
cat("5. 19_candidate_iLAM_gene_set_enrichment_FIXED.csv\n")
