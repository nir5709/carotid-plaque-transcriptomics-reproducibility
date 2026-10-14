# GSE104140 lesion-stage and calcification validation driver
# Run from the repository root.

parts <- file.path(
  "scripts", "04_external_validation", "GSE104140_parts",
  c(
    "01_setup_download_annotation.R",
    "02_normalization_primary_tests.R",
    "03_secondary_figures_summary.R"
  )
)

missing <- parts[!file.exists(parts)]
if (length(missing)) stop("Missing GSE104140 module(s): ", paste(missing, collapse = ", "))

analysis_env <- new.env(parent = globalenv())
for (f in parts) {
  message("Running: ", f)
  source(f, local = analysis_env)
}
