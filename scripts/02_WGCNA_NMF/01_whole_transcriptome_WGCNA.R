# GSE125771 whole-transcriptome WGCNA driver
# Run from the repository root with data/frozen/GSE125771_gene_level_matrix_frozen.csv
# available in the working directory.

parts <- file.path(
  "scripts", "02_WGCNA_NMF", "WGCNA_parts",
  c(
    "01_setup_expression.R",
    "02_traits_network.R",
    "03_module_traits_hubs.R",
    "04_GO_NMF_reproducibility.R"
  )
)

missing <- parts[!file.exists(parts)]
if (length(missing)) stop("Missing WGCNA module(s): ", paste(missing, collapse = ", "))

analysis_env <- new.env(parent = globalenv())
for (f in parts) {
  message("Running: ", f)
  source(f, local = analysis_env)
}

# Publish the frozen module-assignment table at the repository working root as
# the common downstream input expected by the validation workflows. The same
# file also remains inside the WGCNA output directory for provenance.
module_source <- file.path(
  "IDEA2_GSE125771_WHOLE_TRANSCRIPTOME_WGCNA",
  "03_WGCNA_modules",
  "Whole_transcriptome_module_assignments.csv"
)
if (file.exists(module_source)) {
  file.copy(
    module_source,
    "Whole_transcriptome_module_assignments.csv",
    overwrite = TRUE
  )
  message("Published downstream module file: Whole_transcriptome_module_assignments.csv")
}
