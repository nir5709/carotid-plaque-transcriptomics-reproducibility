# GSE210152 macrophage/MNP-state validation driver
# Run from the repository root.
# The final cleaned workflow integrates the g5-versus-5 label normalization that
# was required during the original analysis; the separate repair script is kept
# for provenance of the final historical run.

parts <- file.path(
  "scripts", "05_single_cell", "GSE210152", "parts",
  c(
    "01_setup_helpers.R",
    "02_load_cluster_states.R",
    "03_enrichment_pseudobulk_tests.R",
    "04_figures_summary.R"
  )
)

missing <- parts[!file.exists(parts)]
if (length(missing)) stop("Missing GSE210152 module(s): ", paste(missing, collapse = ", "))

analysis_env <- new.env(parent = globalenv())
for (f in parts) {
  message("Running: ", f)
  source(f, local = analysis_env)
}
