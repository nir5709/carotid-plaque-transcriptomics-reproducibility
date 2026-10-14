###############################################################################
# MANUSCRIPT REPRODUCIBILITY REGRESSION CHECKS
# Run from repository root after completing all manuscript workflows.
###############################################################################

options(stringsAsFactors = FALSE)
REPORT_DIR <- file.path("results", "validation")
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

checks <- list(); k <- 0L
record <- function(dataset, item, metric, observed, expected, tol = 0, note = "") {
  k <<- k + 1L
  o <- suppressWarnings(as.numeric(observed)); e <- suppressWarnings(as.numeric(expected))
  status <- if (is.finite(o) && is.finite(e) && abs(o - e) <= tol) "PASS" else "FAIL"
  checks[[k]] <<- data.frame(
    Dataset=dataset, Item=item, Metric=metric,
    Observed=as.character(observed), Expected=as.character(expected),
    Tolerance=tol, Status=status, Note=note, stringsAsFactors=FALSE
  )
}
skip <- function(dataset, item, note) {
  k <<- k + 1L
  checks[[k]] <<- data.frame(
    Dataset=dataset, Item=item, Metric="file", Observed=NA_character_,
    Expected=NA_character_, Tolerance=NA_real_, Status="SKIP", Note=note,
    stringsAsFactors=FALSE
  )
}

find_output <- function(name, required=TRUE) {
  if (file.exists(name)) return(normalizePath(name, winslash="/"))
  esc <- gsub("\\.", "\\\\.", name)
  hits <- list.files(".", pattern=paste0("^", esc, "$"), recursive=TRUE, full.names=TRUE)
  hits <- hits[!grepl("(^|/)(scripts|docs)/", gsub("\\\\", "/", hits))]
  if (!length(hits)) {
    if (required) stop("Required output not found: ", name)
    return(NA_character_)
  }
  if (length(hits) > 1L) {
    hits <- hits[order(file.info(hits)$mtime, decreasing=TRUE, na.last=TRUE)]
    warning("Multiple copies of ", name, " found; using newest: ", hits[1])
  }
  normalizePath(hits[1], winslash="/")
}
read_required <- function(name) read.csv(find_output(name), check.names=FALSE, stringsAsFactors=FALSE)
row1 <- function(x, col, value) {
  z <- x[as.character(x[[col]]) == value,,drop=FALSE]
  if (nrow(z) != 1L) stop("Expected one row: ", col, " = ", value, "; found ", nrow(z))
  z
}
pair1 <- function(x, a, b) {
  if (all(c("X","Y") %in% names(x))) {
    z <- x[(x$X==a & x$Y==b)|(x$X==b & x$Y==a),,drop=FALSE]
  } else if (all(c("Signature_1","Signature_2") %in% names(x))) {
    z <- x[(x$Signature_1==a & x$Signature_2==b)|(x$Signature_1==b & x$Signature_2==a),,drop=FALSE]
  } else stop("Unsupported correlation-table columns.")
  if (nrow(z) != 1L) stop("Expected one correlation row for ", a, " / ", b)
  z
}

cat("\nCarotid plaque transcriptomics — manuscript regression validation\n")
cat("================================================================\n")

###############################################################################
# GSE125771 discovery
###############################################################################

f <- find_output("signature_score_correlations.csv", required=FALSE)
if (!is.na(f)) {
  x <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
  rho_col <- intersect(c("Spearman_rho","rho"), names(x))[1]
  fdr_col <- intersect(c("FDR_BH_3","FDR","FDR_BH_3correlations"), names(x))[1]
  exp <- list(
    c("Inflammasome","Fibrous",-0.7589118199,2.0816055683e-08),
    c("Inflammasome","Osteogenic",0.7320825516,7.9522742107e-08),
    c("Fibrous","Osteogenic",-0.8223264540,2.2958275518e-10)
  )
  for (v in exp) {
    r <- pair1(x,v[1],v[2])
    record("GSE125771",paste(v[1],v[2],sep="_vs_"),"Spearman_rho",r[[rho_col]],as.numeric(v[3]),1e-6)
    record("GSE125771",paste(v[1],v[2],sep="_vs_"),"FDR_BH_3",r[[fdr_col]],as.numeric(v[4]),1e-10)
  }
} else skip("GSE125771","signature_score_correlations","Run discovery workflow first.")

f <- find_output("DEG_summary.csv", required=FALSE)
if (!is.na(f)) {
  x <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
  exp <- list(Inflammasome=c(326,225,101), Fibrous=c(395,199,196), Osteogenic=c(325,178,147))
  for (nm in names(exp)) {
    r <- row1(x,"Signature",nm)
    record("GSE125771",paste0(nm,"_DEG"),"Significant_total",r$Significant_total,exp[[nm]][1],0)
    record("GSE125771",paste0(nm,"_DEG"),"Upregulated",r$Upregulated,exp[[nm]][2],0)
    record("GSE125771",paste0(nm,"_DEG"),"Downregulated",r$Downregulated,exp[[nm]][3],0)
  }
} else skip("GSE125771","DEG_summary","Run discovery workflow first.")

f <- find_output("Whole_transcriptome_module_assignments.csv", required=FALSE)
if (!is.na(f)) {
  x <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
  nms <- gsub("[^a-z0-9]","",tolower(names(x)))
  idx <- which(nms %in% c("module","modulecolor","modulecolour","color","colour","moduleassignment"))[1]
  if (is.na(idx)) stop("Could not identify WGCNA module column.")
  mods <- tolower(trimws(as.character(x[[idx]])))
  record("GSE125771","WGCNA_Brown","Module_genes",sum(mods=="brown",na.rm=TRUE),4262,0)
  record("GSE125771","WGCNA_Blue","Module_genes",sum(mods=="blue",na.rm=TRUE),5197,0)
} else skip("GSE125771","WGCNA_modules","Run WGCNA workflow first.")

f <- find_output("NMF_module_gene_counts.csv", required=FALSE)
if (!is.na(f)) {
  x <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
  record("GSE125771","NMF_Module_1","Hard_assignment_genes",row1(x,"Module","NMF_Module_1")$Genes,225,0)
  record("GSE125771","NMF_Module_2","Hard_assignment_genes",row1(x,"Module","NMF_Module_2")$Genes,101,0)
} else skip("GSE125771","NMF_module_gene_counts","Run discovery NMF workflow first.")

f <- find_output("MCPcounter_monocytic_partial_correlation.csv", required=FALSE)
if (!is.na(f)) {
  x <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
  record("GSE125771","MCPcounter_partial","Partial_rho",x$Partial_Spearman_rho[1],0.453374,1e-5)
  record("GSE125771","MCPcounter_partial","P_value",x$P_value[1],0.00425496,1e-7)
} else skip("GSE125771","MCPcounter_partial","Partial-correlation output not found.")

###############################################################################
# CIBERSORT-LM22 — canonical historical Windows reproduction
###############################################################################

x <- read_required(
  file.path(
    "data",
    "frozen",
    "GSE125771_CIBERSORT_macrophage_myeloid_composite_correlations_frozen.csv"
  )
)

r <- x[
  x$Variable == "Macrophage_total" &
  x$Signature == "Inflammasome",
  ,
  drop = FALSE
]

if (nrow(r) != 1L) {
  stop("Historical CIBERSORT macrophage/inflammasome row not unique.")
}

record(
  "GSE125771",
  "CIBERSORT_Macrophage_Inflammasome",
  "Spearman_rho",
  r$Spearman_rho,
  0.499437148217636,
  1e-12
)

record(
  "GSE125771",
  "CIBERSORT_Macrophage_Inflammasome",
  "FDR",
  r$FDR,
  0.004656342871635,
  1e-12
)

r <- x[
  x$Variable == "Macrophage_total" &
  x$Signature == "Fibrous",
  ,
  drop = FALSE
]

if (nrow(r) != 1L) {
  stop("Historical CIBERSORT macrophage/fibrous row not unique.")
}

record(
  "GSE125771",
  "CIBERSORT_Macrophage_Fibrous",
  "Spearman_rho",
  r$Spearman_rho,
  -0.396622889305816,
  1e-12
)

record(
  "GSE125771",
  "CIBERSORT_Macrophage_Fibrous",
  "FDR",
  r$FDR,
  0.0203100807856598,
  1e-12
)

x <- read_required(
  file.path(
    "data",
    "frozen",
    "GSE125771_CIBERSORT_signature_correlations_Plt0.05_samples_frozen.csv"
  )
)

record(
  "GSE125771",
  "CIBERSORT_good_fit",
  "N_tests",
  nrow(x),
  66,
  0
)

record(
  "GSE125771",
  "CIBERSORT_good_fit",
  "Minimum_FDR_global",
  min(x$FDR_global, na.rm = TRUE),
  0.23234454681722,
  1e-12
)

###############################################################################
# GSE21545
###############################################################################

x <- read_required("06_signature_correlations.csv")
exp <- list(
  c("Inflammasome","Fibrous",-0.6539933,3.046063e-16),
  c("Inflammasome","Osteogenic",0.2149921,0.01562451),
  c("Fibrous","Osteogenic",-0.3792696,1.783132e-05)
)
for (v in exp) {
  r <- pair1(x,v[1],v[2])
  record("GSE21545",paste(v[1],v[2],sep="_vs_"),"Spearman_rho",r$Spearman_rho,as.numeric(v[3]),1e-6)
  tol <- if (as.numeric(v[4]) < 1e-10) 1e-18 else 1e-7
  record("GSE21545",paste(v[1],v[2],sep="_vs_"),"FDR_BH_3correlations",r$FDR_BH_3correlations,as.numeric(v[4]),tol)
}
x <- read_required("04_signature_coverage.csv")
for (v in list(c("Inflammasome",63),c("Fibrous",7),c("Osteogenic",2))) {
  record("GSE21545",paste0(v[1],"_coverage"),"Present_genes",row1(x,"Signature",v[1])$Present_genes,as.numeric(v[2]),0)
}

###############################################################################
# GSE43292
###############################################################################

x <- read_required("07_prespecified_paired_tests.csv")
exp <- list(Brown_WGCNA=c(0.6106384680,1.3431381438e-05),
            Blue_WGCNA=c(-0.7660466050,7.6590446824e-05),
            Inflammasome_65=c(0.2093919901,1.3431381438e-05))
for (nm in names(exp)) {
  r <- row1(x,"Score",nm)
  record("GSE43292",nm,"N_pairs",r$N_pairs,32,0)
  record("GSE43292",nm,"Median_delta",r$Median_delta,exp[[nm]][1],1e-6)
  record("GSE43292",nm,"Wilcoxon_FDR",r$Wilcoxon_FDR,exp[[nm]][2],1e-9)
}

###############################################################################
# GSE163154 — segment-level arithmetic only; not patient-adjusted confirmation
###############################################################################

x <- read_required("06_prespecified_IPH_tests.csv")
exp <- list(Brown_WGCNA=c(0.775729846,5.427777832e-07),
            Blue_WGCNA=c(-0.851045055,5.427777832e-07),
            Inflammasome_65=c(0.404485595,6.116556968e-07))
for (nm in names(exp)) {
  r <- row1(x,"Score",nm)
  note <- "Exploratory segment-level arithmetic; patient linkage is unavailable."
  record("GSE163154",nm,"Median_difference",r$Median_difference,exp[[nm]][1],1e-6,note)
  record("GSE163154",nm,"Wilcoxon_FDR",r$Wilcoxon_FDR,exp[[nm]][2],1e-10,note)
}

###############################################################################
# GSE104140
###############################################################################

x <- read_required("07_SECONDARY_advanced_vs_DIT_tests.csv")
exp <- list(Brown_WGCNA=c(0.8192411,8.129964e-05),
            Blue_WGCNA=c(-0.6239728,8.129964e-05),
            Inflammasome_65=c(0.3334155,1.182586e-04),
            RUNX2_SPP1=c(1.7012593,8.129964e-05))
for (nm in names(exp)) {
  r <- row1(x,"Score",nm)
  record("GSE104140",paste0(nm,"_Advanced_vs_DIT"),"Median_difference",r$Median_difference_advanced_minus_DIT,exp[[nm]][1],1e-6)
  record("GSE104140",paste0(nm,"_Advanced_vs_DIT"),"FDR_BH_4scores",r$FDR_BH_4scores,exp[[nm]][2],1e-9)
}
x <- read_required("05_PRIMARY_calcified_vs_noncalcified_tests.csv")
r <- row1(x,"Score","RUNX2_SPP1")
record("GSE104140","RUNX2_SPP1_calcified_vs_noncalcified","Median_difference",r$Median_difference_calc_minus_non,0.88589,1e-5)
record("GSE104140","RUNX2_SPP1_calcified_vs_noncalcified","FDR_BH_4scores",r$FDR_BH_4scores,0.69223,1e-5)

###############################################################################
# GSE253903
###############################################################################

x <- read_required("06_primary_target_vs_other_paired_tests.csv")
exp <- list(Brown_WGCNA=c(0.6997114838,0.0025261743),
            Blue_WGCNA=c(0.9187133003,0.0025261743),
            Inflammasome_65=c(0.5802742761,0.0025261743))
for (nm in names(exp)) {
  r <- row1(x,"Score",nm)
  record("GSE253903",nm,"N_donors",r$N_donors,12,0)
  record("GSE253903",nm,"Median_delta",r$Median_delta,exp[[nm]][1],1e-6)
  record("GSE253903",nm,"FDR",r$FDR,exp[[nm]][2],1e-8)
}
x <- read_required("04_celltype_score_summary.csv")
record("GSE253903","pseudobulk_summary","N_broad_cell_types",nrow(x),14,0)
x <- read_required("03_pseudobulk_gene_set_scores.csv")
record("GSE253903","pseudobulk_summary","N_pseudobulks",nrow(x),141,0)

###############################################################################
# GSE210152 — canonical frozen historical reproduction
###############################################################################

gse210152_dir <- "GSE210152_FROZEN_REPRODUCTION"

x <- read_required(
  file.path(
    gse210152_dir,
    "17_candidate_vs_other_paired_tests_FIXED.csv"
  )
)

exp <- list(
  Inflammasome65_score = c(
    0.19549249733257,
    0.1875,
    0.1875
  ),
  iLAM_score = c(
    1.13952315919155,
    0.0625,
    0.0833333333333333
  ),
  TREM2_LAM_score = c(
    0.732711969557989,
    0.0625,
    0.0833333333333333
  ),
  Brown_WGCNA_score = c(
    0.326827205244852,
    0.0625,
    0.0833333333333333
  )
)

for (nm in names(exp)) {

  r <- row1(
    x,
    "Score",
    nm
  )

  record(
    "GSE210152",
    nm,
    "N_donors",
    r$N_donors,
    5,
    0
  )

  record(
    "GSE210152",
    nm,
    "Median_delta",
    r$Median_delta,
    exp[[nm]][1],
    1e-12
  )

  record(
    "GSE210152",
    nm,
    "P_value_exact",
    r$P_value_exact,
    exp[[nm]][2],
    1e-12
  )

  record(
    "GSE210152",
    nm,
    "FDR_exact",
    r$FDR_exact,
    exp[[nm]][3],
    1e-12
  )
}

x <- read_required(
  file.path(
    gse210152_dir,
    "10_candidate_iLAM_donor_representation_FIXED.csv"
  )
)

record(
  "GSE210152",
  "candidate_iLAM",
  "Candidate_cells_total",
  sum(x$Cells),
  604,
  0
)

r <- row1(
  x,
  "Patient",
  "P3"
)

record(
  "GSE210152",
  "candidate_iLAM",
  "P3_candidate_cells",
  r$Cells,
  17,
  0
)

x <- read_required(
  file.path(
    gse210152_dir,
    "19_candidate_iLAM_gene_set_enrichment_FIXED.csv"
  )
)

r <- row1(
  x,
  "GeneSet",
  "Brown_WGCNA"
)

record(
  "GSE210152",
  "candidate_Brown_marker_enrichment",
  "N_cluster_markers",
  r$N_cluster_markers,
  748,
  0
)

record(
  "GSE210152",
  "candidate_Brown_marker_enrichment",
  "Overlap",
  r$Overlap,
  414,
  0
)

record(
  "GSE210152",
  "candidate_Brown_marker_enrichment",
  "Odds_ratio",
  r$Odds_ratio,
  2.87681789118138,
  1e-12
)

record(
  "GSE210152",
  "candidate_Brown_marker_enrichment",
  "FDR_within_gene_set",
  r$FDR_within_gene_set,
  2.580779355979e-42,
  1e-50
)

r <- row1(
  x,
  "GeneSet",
  "Inflammasome_65"
)

record(
  "GSE210152",
  "candidate_Inflammasome_marker_enrichment",
  "Overlap",
  r$Overlap,
  5,
  0
)

record(
  "GSE210152",
  "candidate_Inflammasome_marker_enrichment",
  "Odds_ratio",
  r$Odds_ratio,
  1.44990385375493,
  1e-12
)

record(
  "GSE210152",
  "candidate_Inflammasome_marker_enrichment",
  "FDR_within_gene_set",
  r$FDR_within_gene_set,
  0.519844907699352,
  1e-12
)

###############################################################################
# Report
###############################################################################

report <- do.call(rbind, checks)
write.csv(report, file.path(REPORT_DIR,"reference_output_validation.csv"), row.names=FALSE)
cat("\nValidation summary\n")
print(table(report$Status, useNA="ifany"))
cat("\nReport: ", file.path(REPORT_DIR,"reference_output_validation.csv"), "\n", sep="")

failed <- report[report$Status=="FAIL",,drop=FALSE]
if (nrow(failed)) {
  cat("\nFAILED CHECKS\n"); print(failed,row.names=FALSE)
  stop("One or more manuscript reference checks failed. Investigate before release.")
}
if (any(report$Status=="SKIP")) {
  warning("No checked values failed, but one or more checks were skipped. Resolve SKIP entries before public release.")
} else {
  cat("\nAll configured manuscript reference checks passed.\n")
}
