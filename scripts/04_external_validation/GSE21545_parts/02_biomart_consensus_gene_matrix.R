###############################################################################
# GSE21545 — FROZEN MANUSCRIPT GPL570 PROBE ANNOTATION
#
# The manuscript analysis used a BioMart-assisted consensus GPL570
# probe-to-gene annotation. That exact final mapping is frozen in the
# repository so manuscript results do not depend on future BioMart
# availability or annotation changes.
#
# The frozen mapping contains one row for every GPL570 probe represented
# in the GSE21545 series matrix:
#
#   data/frozen/GSE21545_GPL570_probe_to_gene_frozen.csv
#
# BioMart is therefore NOT queried during the manuscript reproduction run.
###############################################################################

FROZEN_MAP <- file.path(
  "data",
  "frozen",
  "GSE21545_GPL570_probe_to_gene_frozen.csv"
)

if (!file.exists(FROZEN_MAP)) {
  stop(
    "Frozen GSE21545 GPL570 mapping is missing: ",
    FROZEN_MAP,
    "\nThe manuscript reproduction workflow requires this file."
  )
}

message("Using frozen manuscript GPL570 probe-to-gene mapping: ", FROZEN_MAP)

frozen_map <- read.csv(
  FROZEN_MAP,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_cols <- c("ProbeID", "Final_symbol")

if (!all(required_cols %in% names(frozen_map))) {
  stop(
    "Frozen mapping must contain columns: ",
    paste(required_cols, collapse = ", ")
  )
}

frozen_map$ProbeID <- as.character(frozen_map$ProbeID)

if (nrow(frozen_map) != 54675L) {
  stop(
    "Frozen mapping row-count mismatch: expected 54675, observed ",
    nrow(frozen_map)
  )
}

if (anyDuplicated(frozen_map$ProbeID)) {
  stop("Frozen mapping contains duplicated ProbeID values.")
}

if (any(is.na(frozen_map$ProbeID) | !nzchar(frozen_map$ProbeID))) {
  stop("Frozen mapping contains missing/blank ProbeID values.")
}

###############################################################################
# 2. ALIGN FROZEN MAPPING TO CURRENT GEO PROBE MATRIX
###############################################################################

idx <- match(
  probe_ids,
  frozen_map$ProbeID
)

if (anyNA(idx)) {
  missing_probes <- probe_ids[is.na(idx)]

  stop(
    "Current GSE21545 GEO matrix contains ",
    length(missing_probes),
    " probe(s) absent from the frozen manuscript mapping. Examples: ",
    paste(head(missing_probes, 20L), collapse = ", ")
  )
}

final_symbol <- clean_gene(
  frozen_map$Final_symbol[idx]
)

mapped_n <- sum(
  !is.na(final_symbol) &
    nzchar(final_symbol)
)

unmapped_n <- length(final_symbol) - mapped_n

message("Frozen mapping rows aligned: ", length(final_symbol))
message("Mapped probes: ", mapped_n)
message("Unmapped probes: ", unmapped_n)

if (mapped_n != 46806L) {
  stop(
    "Mapped-probe count mismatch: expected 46806, observed ",
    mapped_n
  )
}

if (unmapped_n != 7869L) {
  stop(
    "Unmapped-probe count mismatch: expected 7869, observed ",
    unmapped_n
  )
}

###############################################################################
# 3. WRITE ANNOTATION AUDIT
#
# Historical filename is retained for compatibility with the existing
# validation workflow. The Annotation_source column makes clear that the
# committed frozen mapping is being used rather than a live BioMart query.
###############################################################################

annotation_audit <- data.frame(
  ProbeID = probe_ids,
  Final_symbol = final_symbol,
  Annotation_source = ifelse(
    !is.na(final_symbol) & nzchar(final_symbol),
    "Frozen_manuscript_BioMart_assisted_mapping",
    "Frozen_manuscript_unmapped"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  annotation_audit,
  file.path(
    OUTDIR,
    "02_GSE21545_GPL570_BioMart_annotation_audit.csv"
  ),
  row.names = FALSE,
  na = ""
)

###############################################################################
# 4. CREATE GENE-LEVEL MATRIX
###############################################################################

valid_probe <- !is.na(final_symbol) &
  nzchar(final_symbol)

expr_probe_annotated <- expr_probe[
  valid_probe,
  ,
  drop = FALSE
]

expr_gene_21545 <- limma::avereps(
  expr_probe_annotated,
  ID = final_symbol[valid_probe]
)

rownames(expr_gene_21545) <- clean_gene(
  rownames(expr_gene_21545)
)

expr_gene_21545 <- expr_gene_21545[
  !is.na(rownames(expr_gene_21545)),
  ,
  drop = FALSE
]

if (nrow(expr_gene_21545) != 21875L) {
  stop(
    "Gene-level matrix row-count mismatch: expected 21875, observed ",
    nrow(expr_gene_21545)
  )
}

if (ncol(expr_gene_21545) != 126L) {
  stop(
    "Gene-level matrix sample-count mismatch: expected 126, observed ",
    ncol(expr_gene_21545)
  )
}

write.csv(
  data.frame(
    Gene = rownames(expr_gene_21545),
    expr_gene_21545,
    check.names = FALSE
  ),
  file.path(
    OUTDIR,
    "03_GSE21545_126plaques_BioMart_gene_level.csv"
  ),
  row.names = FALSE
)

###############################################################################
# 5. PROVENANCE RECORD
###############################################################################

writeLines(
  c(
    "GSE21545 frozen GPL570 manuscript mapping",
    "==========================================",
    "",
    paste("Mapping file:", FROZEN_MAP),
    "Historical origin: final BioMart-assisted GPL570 probe annotation",
    "Frozen rows: 54675",
    "Mapped probes: 46806",
    "Unmapped probes: 7869",
    "Expected gene-level genes: 21875",
    "Expected plaque samples: 126",
    paste0(
      "Frozen mapping SHA256: ",
      "7a842b4b75746418ddd74abc5e64f44e1b30f52a667f4c654108ba2c067af1f6"
    ),
    "",
    "Live BioMart is not queried during manuscript reproduction."
  ),
  file.path(
    OUTDIR,
    "02_GSE21545_frozen_mapping_provenance.txt"
  )
)

message(
  "GSE21545 gene-level matrix reconstructed from frozen manuscript mapping."
)
