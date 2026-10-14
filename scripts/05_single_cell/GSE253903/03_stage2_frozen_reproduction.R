###############################################################################
# GSE253903 — CANONICAL FROZEN MANUSCRIPT REPRODUCTION
#
# Reconstructs the manuscript Stage-2 cell-type localization results from the
# frozen donor × broad-cell-type pseudobulk logCPM matrix and metadata.
#
# This canonical path does NOT require the ~3.7-GB historical Seurat object.
# The full Stage-1/Stage-2 scripts remain available as the raw-data
# reconstruction workflow.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260923)

###############################################################################
# 0. SETTINGS
###############################################################################

EXPR_FILE <- "data/frozen/GSE253903_pseudobulk_logCPM_frozen.csv"
META_FILE <- "data/frozen/GSE253903_pseudobulk_metadata_frozen.csv"
MODULE_FILE <- "Whole_transcriptome_module_assignments.csv"

OUTDIR <- "GSE253903_STAGE2_CELLTYPE_VALIDATION"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MIN_DONORS_PER_CELLTYPE <- 4L

BROWN_PRIMARY_CELLTYPES <- c("Macrophage", "Monocyte")
BLUE_PRIMARY_CELLTYPES  <- c("SMC", "Fibroblast")

###############################################################################
# 1. REQUIRED PACKAGES
###############################################################################

required_pkgs <- c("dplyr", "ggplot2", "pheatmap")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs)) {
  stop(
    "Missing package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script."
  )
}

###############################################################################
# 2. LOCKED INFLAMMASOME SIGNATURE
###############################################################################

inflammasome_65 <- c(
  "ABHD17A","ABHD8","ADCY7","ATAT1","BRCC3","BTK","CARD8","CASR",
  "CCDC92","CD36","CPTP","CSNK1A1","DDX3X","DHX33","EIF2AK2","ELP6",
  "FBXL2","GBP2","GBP5","GKN2","HSPA8","IFI16","IRAK1","IRAK2",
  "IRAK4","IRGM","KCNJ8","KCNK13","KCNK6","LAMP2","LATS1","LATS2",
  "MAPK8","MARCHF5","MARK4","MAVS","MEFV","MYD88","NEK7","NLRC3",
  "NLRP1","NLRP3","NLRP6","OGT","P2RX7","PLCG2","PPP2CA","PRKD1",
  "PTPN22","PYCARD","PYDC1","PYDC2","SIRT2","STMP1","TLR4","TLR6",
  "TLR9","TREM2","TRIM11","TRIM31","TRIM65","USP50","ZDHHC1",
  "ZDHHC12","ZDHHC5"
)

###############################################################################
# 3. HELPERS
###############################################################################

read_modules <- function(path) {

  x <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  clean_names <- gsub(
    "[^a-z0-9]",
    "",
    tolower(names(x))
  )

  gene_idx <- which(
    clean_names %in%
      c("gene","genes","symbol","genesymbol","geneid")
  )

  module_idx <- which(
    clean_names %in%
      c(
        "module",
        "modulecolor",
        "modulecolour",
        "color",
        "moduleassignment"
      )
  )

  if (!length(gene_idx) || !length(module_idx)) {
    stop("Could not identify gene/module columns in module file.")
  }

  gene_col <- names(x)[gene_idx[1]]
  module_col <- names(x)[module_idx[1]]

  genes <- toupper(trimws(as.character(x[[gene_col]])))
  modules <- tolower(trimws(as.character(x[[module_col]])))

  list(
    brown = unique(genes[modules == "brown"]),
    blue  = unique(genes[modules == "blue"])
  )
}

row_z_score <- function(mat) {

  mat <- as.matrix(mat)

  z <- t(
    apply(
      mat,
      1,
      function(v) {

        s <- sd(v, na.rm = TRUE)

        if (is.na(s) || s == 0) {
          return(rep(0, length(v)))
        }

        (v - mean(v, na.rm = TRUE)) / s
      }
    )
  )

  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)

  z
}

score_gene_set <- function(logcpm, genes) {

  genes <- intersect(
    toupper(genes),
    rownames(logcpm)
  )

  if (length(genes) < 2L) {
    return(rep(NA_real_, ncol(logcpm)))
  }

  sub <- logcpm[
    genes,
    ,
    drop = FALSE
  ]

  variable <- apply(
    sub,
    1,
    var,
    na.rm = TRUE
  )

  genes_use <- names(variable)[
    !is.na(variable) &
      variable > 0
  ]

  if (length(genes_use) < 2L) {
    return(rep(NA_real_, ncol(logcpm)))
  }

  z <- row_z_score(
    sub[
      genes_use,
      ,
      drop = FALSE
    ]
  )

  colMeans(
    z,
    na.rm = TRUE
  )
}

paired_target_delta <- function(
  scores,
  score_col,
  target_celltypes,
  comparison_name
) {

  patients <- unique(scores$Sample)

  rows <- list()
  rr <- 1L

  for (pat in patients) {

    tmp <- scores[
      scores$Sample == pat,
      ,
      drop = FALSE
    ]

    target <- tmp[
      tmp$Broad_label %in% target_celltypes,
      score_col,
      drop = TRUE
    ]

    other <- tmp[
      !tmp$Broad_label %in% target_celltypes,
      score_col,
      drop = TRUE
    ]

    target <- target[is.finite(target)]
    other  <- other[is.finite(other)]

    if (!length(target) || !length(other)) {
      next
    }

    rows[[rr]] <- data.frame(
      Comparison = comparison_name,
      Score = score_col,
      Sample = pat,
      Target_mean = mean(target),
      Other_mean = mean(other),
      Delta_target_minus_other =
        mean(target) - mean(other),
      stringsAsFactors = FALSE
    )

    rr <- rr + 1L
  }

  if (!length(rows)) {
    return(data.frame())
  }

  do.call(
    rbind,
    rows
  )
}

save_boxplot <- function(
  df,
  score_col,
  title_text,
  file_stem
) {

  med_order <- df |>
    dplyr::group_by(Broad_label) |>
    dplyr::summarise(
      med = median(
        .data[[score_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(med)
    )

  df$Broad_label <- factor(
    df$Broad_label,
    levels = med_order$Broad_label
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Broad_label,
      y = .data[[score_col]]
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      height = 0,
      alpha = 0.65,
      size = 1.5
    ) +
    ggplot2::theme_bw(
      base_size = 11
    ) +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Pseudobulk gene-set score",
      title = title_text
    )

  ggplot2::ggsave(
    paste0(file_stem, ".png"),
    p,
    width = 10,
    height = 6,
    dpi = 300
  )

  ggplot2::ggsave(
    paste0(file_stem, ".pdf"),
    p,
    width = 10,
    height = 6
  )
}

###############################################################################
# 4. LOAD FROZEN PSEUDOBULK INPUTS
###############################################################################

if (!file.exists(EXPR_FILE)) {
  stop("Frozen expression matrix missing: ", EXPR_FILE)
}

if (!file.exists(META_FILE)) {
  stop("Frozen pseudobulk metadata missing: ", META_FILE)
}

x <- read.csv(
  EXPR_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!identical(names(x)[1], "Gene")) {
  stop("First frozen expression column must be Gene.")
}

logcpm <- as.matrix(
  x[, -1, drop = FALSE]
)

storage.mode(logcpm) <- "numeric"

rownames(logcpm) <- toupper(
  trimws(
    as.character(x$Gene)
  )
)

pb_meta <- read.csv(
  META_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_meta <- c(
  "Pseudobulk",
  "Sample",
  "Broad_label"
)

if (!all(required_meta %in% names(pb_meta))) {
  stop("Frozen metadata is missing required columns.")
}

if (!identical(
  colnames(logcpm),
  pb_meta$Pseudobulk
)) {
  stop(
    "Frozen expression columns and metadata pseudobulks are not aligned."
  )
}

stopifnot(
  nrow(logcpm) == 13506L,
  ncol(logcpm) == 141L,
  nrow(pb_meta) == 141L,
  length(unique(pb_meta$Sample)) == 12L,
  length(unique(pb_meta$Broad_label)) == 14L
)

write.csv(
  pb_meta,
  file.path(
    OUTDIR,
    "02_pseudobulk_metadata.csv"
  ),
  row.names = FALSE
)

###############################################################################
# 5. LOAD LOCKED BROWN / BLUE MODULES
###############################################################################

mods <- read_modules(MODULE_FILE)

brown_genes <- mods$brown
blue_genes  <- mods$blue

stopifnot(
  length(brown_genes) == 4262L,
  length(blue_genes) == 5197L,
  length(inflammasome_65) == 65L
)

###############################################################################
# 6. SCORES / TESTS / RANKING
###############################################################################

source(
  "scripts/05_single_cell/GSE253903/stage2_parts/03_scores_tests_ranking.R"
)

###############################################################################
# 7. FIGURES
###############################################################################

save_boxplot(
  scores,
  "Brown_WGCNA",
  "GSE253903: Brown inflammatory WGCNA module by cell type",
  file.path(
    OUTDIR,
    "Figure_Brown_by_celltype"
  )
)

save_boxplot(
  scores,
  "Blue_WGCNA",
  "GSE253903: Blue fibrous WGCNA module by cell type",
  file.path(
    OUTDIR,
    "Figure_Blue_by_celltype"
  )
)

save_boxplot(
  scores,
  "Inflammasome_65",
  "GSE253903: locked 65-gene inflammasome score by cell type",
  file.path(
    OUTDIR,
    "Figure_Inflammasome65_by_celltype"
  )
)

heat <- as.data.frame(
  score_summary[
    ,
    c(
      "Broad_label",
      "Brown_median",
      "Blue_median",
      "Inflammasome65_median"
    ),
    drop = FALSE
  ]
)

rownames(heat) <- heat$Broad_label
heat$Broad_label <- NULL

png(
  file.path(
    OUTDIR,
    "Figure_pseudobulk_module_heatmap.png"
  ),
  width = 1800,
  height = 1400,
  res = 180
)

pheatmap::pheatmap(
  as.matrix(heat),
  scale = "column",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main =
    "GSE253903 donor-pseudobulk module localization"
)

dev.off()

pdf(
  file.path(
    OUTDIR,
    "Figure_pseudobulk_module_heatmap.pdf"
  ),
  width = 8,
  height = 7
)

pheatmap::pheatmap(
  as.matrix(heat),
  scale = "column",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main =
    "GSE253903 donor-pseudobulk module localization"
)

dev.off()

###############################################################################
# 8. COVERAGE
###############################################################################

coverage <- data.frame(
  GeneSet = c(
    "Brown_WGCNA",
    "Blue_WGCNA",
    "Inflammasome_65"
  ),
  Total_genes = c(
    length(brown_genes),
    length(blue_genes),
    length(inflammasome_65)
  ),
  Present_in_filtered_pseudobulk = c(
    length(
      intersect(
        brown_genes,
        rownames(logcpm)
      )
    ),
    length(
      intersect(
        blue_genes,
        rownames(logcpm)
      )
    ),
    length(
      intersect(
        inflammasome_65,
        rownames(logcpm)
      )
    )
  ),
  stringsAsFactors = FALSE
)

coverage$Coverage_fraction <-
  coverage$Present_in_filtered_pseudobulk /
  coverage$Total_genes

write.csv(
  coverage,
  file.path(
    OUTDIR,
    "08_gene_set_coverage.csv"
  ),
  row.names = FALSE
)

###############################################################################
# 9. SUMMARY
###############################################################################

summary_file <- file.path(
  OUTDIR,
  "FINAL_GSE253903_STAGE2_SUMMARY.txt"
)

sink(summary_file)

cat(
  "GSE253903 STAGE 2 — FROZEN MANUSCRIPT REPRODUCTION\n"
)

cat(
  "===================================================\n\n"
)

cat(
  "Frozen logCPM matrix:",
  EXPR_FILE,
  "\n"
)

cat(
  "Frozen pseudobulk metadata:",
  META_FILE,
  "\n"
)

cat(
  "Genes:",
  nrow(logcpm),
  "\n"
)

cat(
  "Pseudobulks:",
  ncol(logcpm),
  "\n"
)

cat(
  "Donors:",
  length(unique(pb_meta$Sample)),
  "\n"
)

cat(
  "Broad cell types:",
  length(unique(pb_meta$Broad_label)),
  "\n\n"
)

cat("GENE-SET COVERAGE\n")
print(
  coverage,
  row.names = FALSE
)

cat("\n\nCELL-TYPE SCORE SUMMARY\n")
print(
  score_summary,
  row.names = FALSE
)

cat("\n\nPRIMARY DONOR-BALANCED TESTS\n")
print(
  paired_tests,
  row.names = FALSE
)

cat("\n\nLOCALIZATION RANKING\n")
print(
  ranking,
  row.names = FALSE
)

cat(
  "\n\nCanonical reproduction begins from the frozen donor × broad-cell-type pseudobulk logCPM matrix.\n"
)

cat(
  "The full raw-data reconstruction workflow remains available separately and includes Stage-1 Seurat clustering, locked cluster annotation, cell exclusions, pseudobulk aggregation and edgeR normalization.\n"
)

sink()

capture.output(
  sessionInfo(),
  file = file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

cat("\n=============================================\n")
cat("GSE253903 FROZEN REPRODUCTION COMPLETE\n")
cat("=============================================\n")
