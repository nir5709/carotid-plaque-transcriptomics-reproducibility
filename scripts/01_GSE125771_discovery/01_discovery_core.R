# GSE125771 discovery-analysis driver
# Run from the repository root with data/frozen/GSE125771_gene_level_matrix_frozen.csv in
# the working directory. Historical CIBERSORT/LM22 results are preserved only
# as author-derived frozen summary tables under data/frozen/; no executable
# CIBERSORT/LM22 workflow is redistributed in this archive.

parts <- file.path(
  "scripts", "01_GSE125771_discovery", "parts",
  c(
    "01_setup_signatures_scores.R",
    "02_DEG_GO_MCPcounter.R",
    "03_NMF_summary.R"
  )
)

missing <- parts[!file.exists(parts)]
if (length(missing)) stop("Missing discovery module(s): ", paste(missing, collapse = ", "))

analysis_env <- new.env(parent = globalenv())
for (f in parts) {
  message("Running: ", f)
  source(f, local = analysis_env)
}
