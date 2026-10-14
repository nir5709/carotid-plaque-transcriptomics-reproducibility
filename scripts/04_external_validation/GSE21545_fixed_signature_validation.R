# GSE21545 final fixed-signature validation driver
# Run from repository root. Internet access is required only to retrieve GEO
# data on a fresh run. The manuscript probe-to-gene mapping is frozen under
# data/frozen/ for deterministic reproduction.

parts <- file.path(
  "scripts", "04_external_validation", "GSE21545_parts",
  c(
    "01_setup_download_local_annotation.R",
    "02_biomart_consensus_gene_matrix.R",
    "03_scores_statistics_summary.R"
  )
)

missing <- parts[!file.exists(parts)]
if (length(missing)) stop("Missing GSE21545 module(s): ", paste(missing, collapse = ", "))

analysis_env <- new.env(parent = globalenv())
for (f in parts) {
  message("Running: ", f)
  source(f, local = analysis_env)
}
