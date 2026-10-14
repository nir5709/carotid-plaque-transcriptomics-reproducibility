###############################################################################
# GSE210152 — FROZEN HISTORICAL MANUSCRIPT REPRODUCTION
#
# Purpose:
#   Reconstruct the manuscript-level GSE210152 results deterministically from
#   frozen historical intermediate tables.
#
# Historical upstream environment:
#   R 4.4.2
#   Seurat 5.5.0
#   Harmony 1.2.4
#
# The historical clustering itself is NOT rerun here.
# Instead, frozen outputs from that clustering are staged and passed through
# the validated g5 -> 5 cluster-label repair.
###############################################################################

options(stringsAsFactors = FALSE)

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)

if (length(file_arg) != 1L) {
  stop("Could not determine script location.")
}

script_path <- normalizePath(
  sub("^--file=", "", file_arg),
  mustWork = TRUE
)

script_dir <- dirname(script_path)

repo_root <- normalizePath(
  file.path(script_dir, "..", "..", ".."),
  mustWork = TRUE
)

frozen_dir <- file.path(
  repo_root,
  "data",
  "frozen"
)

repair_script <- file.path(
  script_dir,
  "02_g5_cluster_label_repair.R"
)

if (!file.exists(repair_script)) {
  stop("Repair script not found: ", repair_script)
}

###############################################################################
# Frozen historical inputs
###############################################################################

inputs <- c(
  "07_MNP_cluster_state_scores.csv" =
    "GSE210152_MNP_cluster_state_scores_frozen.csv",

  "08_MNP_cluster_key_gene_average_expression.csv" =
    "GSE210152_MNP_cluster_key_gene_average_expression_frozen.csv",

  "09_MNP_cluster_donor_representation.csv" =
    "GSE210152_MNP_cluster_donor_representation_frozen.csv",

  "12_gene_set_enrichment_by_MNP_cluster.csv" =
    "GSE210152_gene_set_enrichment_by_MNP_cluster_frozen.csv",

  "15_pseudobulk_gene_set_scores.csv" =
    "GSE210152_pseudobulk_gene_set_scores_frozen.csv"
)

source_paths <- file.path(
  frozen_dir,
  unname(inputs)
)

if (any(!file.exists(source_paths))) {
  stop(
    "Missing frozen input(s):\n",
    paste(
      source_paths[!file.exists(source_paths)],
      collapse = "\n"
    )
  )
}

###############################################################################
# Output directory
#
# Can be overridden, for example:
#
# GSE210152_FROZEN_OUTDIR=/tmp/GSE210152_test \
#   Rscript scripts/05_single_cell/GSE210152/03_frozen_historical_reproduction.R
###############################################################################

outdir_env <- Sys.getenv(
  "GSE210152_FROZEN_OUTDIR",
  unset = ""
)

if (nzchar(outdir_env)) {
  outdir <- path.expand(outdir_env)
} else {
  outdir <- file.path(
    repo_root,
    "GSE210152_FROZEN_REPRODUCTION"
  )
}

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################################
# Stage frozen inputs using the filenames expected by the historical repair
###############################################################################

for (i in seq_along(inputs)) {

  src <- source_paths[i]

  dest <- file.path(
    outdir,
    names(inputs)[i]
  )

  ok <- file.copy(
    src,
    dest,
    overwrite = TRUE
  )

  if (!ok) {
    stop(
      "Failed to stage frozen input: ",
      src
    )
  }
}

###############################################################################
# Preserve historical session information alongside reproduction outputs
###############################################################################

historical_session <- file.path(
  frozen_dir,
  "GSE210152_historical_sessionInfo.txt"
)

if (!file.exists(historical_session)) {
  stop(
    "Historical sessionInfo not found: ",
    historical_session
  )
}

file.copy(
  historical_session,
  file.path(
    outdir,
    "historical_sessionInfo.txt"
  ),
  overwrite = TRUE
)

###############################################################################
# Run validated g5 -> 5 repair in an isolated environment
###############################################################################

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)

setwd(outdir)

repair_env <- new.env(
  parent = globalenv()
)

sys.source(
  repair_script,
  envir = repair_env
)

###############################################################################
# Confirm expected manuscript outputs were generated
###############################################################################

expected_outputs <- c(
  "10_candidate_iLAM_donor_representation_FIXED.csv",
  "15_pseudobulk_gene_set_scores_FIXED.csv",
  "16_candidate_vs_other_within_donor_deltas_FIXED.csv",
  "17_candidate_vs_other_paired_tests_FIXED.csv",
  "18_candidate_iLAM_state_scores_FIXED.csv",
  "19_candidate_iLAM_gene_set_enrichment_FIXED.csv",
  "20_candidate_iLAM_key_gene_expression_FIXED.csv"
)

missing_outputs <- expected_outputs[
  !file.exists(
    file.path(
      outdir,
      expected_outputs
    )
  )
]

if (length(missing_outputs) > 0L) {
  stop(
    "Frozen reproduction did not generate:\n",
    paste(
      missing_outputs,
      collapse = "\n"
    )
  )
}

cat("\n")
cat("============================================================\n")
cat("GSE210152 FROZEN HISTORICAL REPRODUCTION: COMPLETE\n")
cat("============================================================\n")
cat("Output directory:\n", outdir, "\n")
cat("Generated manuscript tables:", length(expected_outputs), "\n")
