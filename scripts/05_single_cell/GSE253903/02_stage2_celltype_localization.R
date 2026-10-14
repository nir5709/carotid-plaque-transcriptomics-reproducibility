# GSE253903 Stage 2 driver
# Run this script from the repository root.
# The original Stage 2 analysis is split at section boundaries for readability
# while preserving the donor-aware analysis sequence.

parts <- file.path(
  "scripts", "05_single_cell", "GSE253903", "stage2_parts",
  c(
    "01_setup_helpers.R",
    "02_pseudobulk_normalization.R",
    "03_scores_tests_ranking.R",
    "04_figures_outputs.R"
  )
)

missing <- parts[!file.exists(parts)]
if (length(missing)) stop("Missing Stage 2 module(s): ", paste(missing, collapse = ", "))

analysis_env <- new.env(parent = globalenv())
for (f in parts) {
  message("Running: ", f)
  source(f, local = analysis_env)
}
