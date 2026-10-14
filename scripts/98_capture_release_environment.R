# Capture the software environment used for final release verification.
# Run from the repository root after the release-candidate workflows and
# regression validator have completed successfully.

options(stringsAsFactors = FALSE)

OUTDIR <- "environment"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

capture_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

get_git_commit <- function() {
  out <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  )
  if (length(out) != 1L || !grepl("^[0-9a-fA-F]{40}$", out)) return(NA_character_)
  out
}

git_commit <- get_git_commit()

bioc_version <- if (requireNamespace("BiocManager", quietly = TRUE)) {
  as.character(BiocManager::version())
} else {
  NA_character_
}

# Keep only non-identifying operating-system fields. In particular, do not
# archive Sys.info() fields such as nodename, login, user or effective_user.
sys <- Sys.info()
safe_system <- sys[intersect(c("sysname", "release", "version", "machine"), names(sys))]

system_info_file <- file.path(OUTDIR, "release_environment.txt")

capture.output(
  {
    cat("Carotid plaque transcriptomics — release verification environment\n")
    cat("================================================================\n\n")
    cat("Captured: ", capture_time, "\n", sep = "")
    cat("Git commit: ", ifelse(is.na(git_commit), "unavailable", git_commit), "\n", sep = "")
    cat("R: ", R.version.string, "\n", sep = "")
    cat("R platform: ", R.version$platform, "\n", sep = "")
    cat("Bioconductor: ", ifelse(is.na(bioc_version), "unavailable", bioc_version), "\n\n", sep = "")
    cat("Operating-system information (non-identifying fields only)\n")
    print(safe_system)
    cat("\nSession information\n")
    print(sessionInfo())
  },
  file = system_info_file
)

packages <- c(
  "data.table", "ggplot2", "dplyr", "tidyr", "patchwork", "ragg", "png",
  "pheatmap", "ppcor", "NMF", "WGCNA", "msigdbr", "Seurat",
  "SeuratObject", "Matrix", "harmony", "R.utils", "remotes",
  "GEOquery", "Biobase", "GSVA", "limma", "edgeR", "AnnotationDbi",
  "org.Hs.eg.db", "hgu133plus2.db", "hugene10sttranscriptcluster.db",
  "illuminaHumanv2.db", "clusterProfiler", "biomaRt", "zellkonverter",
  "SingleCellExperiment", "SummarizedExperiment", "preprocessCore", "impute",
  "GO.db", "MCPcounter"
)

version_or_na <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

description_field <- function(pkg, field) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  d <- tryCatch(utils::packageDescription(pkg), error = function(e) NULL)
  if (is.null(d) || is.null(d[[field]]) || !nzchar(as.character(d[[field]]))) {
    return(NA_character_)
  }
  as.character(d[[field]])
}

package_table <- data.frame(
  Package = packages,
  Version = vapply(packages, version_or_na, character(1)),
  Available = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
  RemoteUsername = vapply(packages, description_field, character(1), field = "RemoteUsername"),
  RemoteRepo = vapply(packages, description_field, character(1), field = "RemoteRepo"),
  RemoteSha = vapply(packages, description_field, character(1), field = "RemoteSha"),
  stringsAsFactors = FALSE
)

write.csv(
  package_table,
  file.path(OUTDIR, "release_package_versions.csv"),
  row.names = FALSE,
  na = ""
)

validation_source <- file.path("results", "validation", "reference_output_validation.csv")
validation_copy <- file.path(OUTDIR, "release_reference_output_validation.csv")
validation_summary <- "Validation report unavailable"

if (file.exists(validation_source)) {
  validation <- read.csv(validation_source, stringsAsFactors = FALSE, check.names = FALSE)
  file.copy(validation_source, validation_copy, overwrite = TRUE)
  if ("Status" %in% names(validation)) {
    counts <- table(validation$Status)
    get_n <- function(x) if (x %in% names(counts)) unname(counts[[x]]) else 0L
    validation_summary <- paste0(
      "Validation summary: PASS=", get_n("PASS"),
      "; FAIL=", get_n("FAIL"),
      "; SKIP=", get_n("SKIP")
    )
  }
}

writeLines(
  c(
    paste0("Captured: ", capture_time),
    paste0("Git commit: ", ifelse(is.na(git_commit), "unavailable", git_commit)),
    paste0("R: ", R.version.string),
    paste0("R platform: ", R.version$platform),
    paste0("Bioconductor: ", ifelse(is.na(bioc_version), "unavailable", bioc_version)),
    validation_summary
  ),
  file.path(OUTDIR, "release_verification_metadata.txt")
)

cat("Release environment captured.\n")
cat("  ", system_info_file, "\n", sep = "")
cat("  ", file.path(OUTDIR, "release_package_versions.csv"), "\n", sep = "")
cat("  ", file.path(OUTDIR, "release_verification_metadata.txt"), "\n", sep = "")
if (file.exists(validation_copy)) {
  cat("  ", validation_copy, "\n", sep = "")
}
