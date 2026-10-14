# GSE253903 Stage 1 driver
# Run this script from the repository root.
# The original Stage 1 analysis has been split at section boundaries only so the
# repository remains readable while preserving the analysis sequence.

parts <- file.path(
  "scripts", "05_single_cell", "GSE253903", "stage1_parts",
  c(
    "01_setup_download_metadata.R",
    "02_qc_merge_embedding.R",
    "03_markers_outputs.R"
  )
)

missing <- parts[!file.exists(parts)]
if (length(missing)) stop("Missing Stage 1 module(s): ", paste(missing, collapse = ", "))

# Ensure all CRAN dependencies are available before the historical setup module
# runs. This also avoids the original setup script's package-install loop being
# asked to resolve individual missing packages during analysis.
cran_stage1 <- c("Seurat","SeuratObject","Matrix","data.table","ggplot2","patchwork","harmony")
missing_cran <- cran_stage1[
  !vapply(cran_stage1, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# Use one shared analysis environment. This keeps objects available between
# modules while preventing the first module's rm(list = ls()) from deleting the
# driver's own bookkeeping objects.
analysis_env <- new.env(parent = globalenv())

for (f in parts) {
  message("Running: ", f)
  source(f, local = analysis_env)
}
