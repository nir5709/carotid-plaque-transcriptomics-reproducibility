# Install R/Bioconductor dependencies used by the manuscript workflows.
# Run once in a clean R environment before executing the analysis drivers.

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "data.table", "ggplot2", "dplyr", "tidyr", "patchwork", "ragg", "png",
  "pheatmap", "ppcor", "NMF", "WGCNA", "msigdbr", "Seurat",
  "SeuratObject", "Matrix", "harmony", "R.utils", "remotes"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "GEOquery", "Biobase", "GSVA", "limma", "edgeR", "AnnotationDbi",
  "org.Hs.eg.db", "hgu133plus2.db", "hugene10sttranscriptcluster.db",
  "illuminaHumanv2.db", "clusterProfiler", "biomaRt", "zellkonverter",
  "SingleCellExperiment", "SummarizedExperiment", "preprocessCore", "impute",
  "GO.db"
)

cran_missing <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(cran_missing)) install.packages(cran_missing)

bioc_missing <- bioc_pkgs[
  !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(bioc_missing)) {
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
}

# MCP-counter provenance for the validated release environment:
# ebecht/MCPcounter commit b6eac73e91c246fcff0bb1a5c68a816cd588fc48.
# Pin the installation so a future GitHub HEAD change cannot silently alter the
# canonical composition analysis.
if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  try(
    remotes::install_github(
      "ebecht/MCPcounter@b6eac73e91c246fcff0bb1a5c68a816cd588fc48",
      subdir = "Source",
      upgrade = "never"
    ),
    silent = TRUE
  )
}

cat("\nDependency check\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs, "MCPcounter")
status <- data.frame(
  Package = all_pkgs,
  Available = vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE),
  stringsAsFactors = FALSE
)
print(status, row.names = FALSE)

if (any(!status$Available)) {
  warning("One or more dependencies are still unavailable; inspect the table above.")
}
