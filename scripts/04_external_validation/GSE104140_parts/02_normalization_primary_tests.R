###############################################################################
# 8. VERIFY FROZEN EXPRESSION / METADATA ALIGNMENT
###############################################################################

if (nrow(expr_gene) != 20861L) {
  stop(
    "Expected 20861 genes in frozen GSE104140 matrix; found ",
    nrow(expr_gene),
    "."
  )
}

if (ncol(expr_gene) != 32L) {
  stop(
    "Expected 32 samples in frozen GSE104140 matrix; found ",
    ncol(expr_gene),
    "."
  )
}

if (!identical(colnames(expr_gene), meta$GSM)) {
  stop(
    "Frozen GSE104140 matrix and metadata are not in identical sample order."
  )
}

###############################################################################
# 9. NORMALIZATION STATUS
###############################################################################

# expr_gene is already the historical manuscript TMM-normalized log2-CPM
# gene-level matrix. No re-normalization or contemporary gene annotation is
# performed in the canonical manuscript-reproduction workflow.

cat(
  "\nUsing frozen historical TMM/log2-CPM gene-level matrix.\n"
)

###############################################################################
# 10. FROZEN GENE-SET COVERAGE
###############################################################################

coverage <- do.call(
  rbind,
  lapply(names(GENE_SETS), function(nm) {
    gs <- unique(toupper(GENE_SETS[[nm]]))
    present <- intersect(gs, rownames(expr_gene))
    missing <- setdiff(gs, rownames(expr_gene))

    data.frame(
      GeneSet = nm,
      Total_genes = length(gs),
      Present_genes = length(present),
      Missing_genes = length(missing),
      Percent_present = 100 * length(present) / length(gs),
      Missing_gene_names = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  coverage,
  file.path(OUTDIR, "03_gene_set_coverage.csv"),
  row.names = FALSE
)

cat("\nFrozen gene-set coverage after expression filtering:\n")
print(coverage, row.names = FALSE)

###############################################################################
# 11. SCORING
###############################################################################

gene_z <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  s <- apply(mat, 1, sd, na.rm = TRUE)
  good <- is.finite(s) & s > 0

  z <- sweep(mat[good, , drop = FALSE], 1, mu[good], "-")
  z <- sweep(z, 1, s[good], "/")
  z
}

score_frozen_sets <- function(mat, sets) {
  z <- gene_z(mat)

  ans <- data.frame(
    GSM = colnames(z),
    stringsAsFactors = FALSE
  )

  for (nm in names(sets)) {
    present <- intersect(
      unique(toupper(sets[[nm]])),
      rownames(z)
    )

    if (!length(present)) {
      ans[[nm]] <- NA_real_
    } else {
      ans[[nm]] <- colMeans(
        z[present, , drop = FALSE],
        na.rm = TRUE
      )
    }
  }

  ans$Brown_minus_Blue <- ans$Brown_WGCNA - ans$Blue_WGCNA
  ans
}

###############################################################################
# 12. PRIMARY ANALYSIS — STANDARDIZE WITHIN THE 19 ADVANCED FIBROATHEROMAS
###############################################################################

primary_gsm <- meta$GSM[
  meta$Pathology_group %in% c("Calcified_FA", "Noncalcified_FA")
]

if (length(primary_gsm) != 19L) {
  stop("Expected 19 primary fibroatheroma samples; found ", length(primary_gsm), ".")
}

expr_primary <- expr_gene[, primary_gsm, drop = FALSE]
scores_primary <- score_frozen_sets(expr_primary, GENE_SETS)

scores_primary <- merge(
  meta[, c("GSM", "Title", "Pathology_group")],
  scores_primary,
  by = "GSM",
  sort = FALSE
)

scores_primary <- scores_primary[
  match(primary_gsm, scores_primary$GSM),
  ,
  drop = FALSE
]

scores_primary$Pathology_group <- factor(
  scores_primary$Pathology_group,
  levels = c("Noncalcified_FA", "Calcified_FA")
)

write.csv(
  scores_primary,
  file.path(OUTDIR, "04_PRIMARY_19sample_frozen_programme_scores.csv"),
  row.names = FALSE
)

###############################################################################
# 13. PRIMARY STATISTICAL TESTS
###############################################################################

score_names <- c(
  "Brown_WGCNA",
  "Blue_WGCNA",
  "Inflammasome_65",
  "RUNX2_SPP1"
)

test_one <- function(sc) {
  x_calc <- scores_primary[
    scores_primary$Pathology_group == "Calcified_FA",
    sc
  ]

  x_non <- scores_primary[
    scores_primary$Pathology_group == "Noncalcified_FA",
    sc
  ]

  w <- wilcox.test(
    x_calc,
    x_non,
    alternative = "two.sided",
    exact = FALSE
  )

  tt <- t.test(
    x_calc,
    x_non,
    var.equal = FALSE
  )

  data.frame(
    Score = sc,
    N_calcified = length(x_calc),
    N_noncalcified = length(x_non),
    Median_calcified = median(x_calc, na.rm = TRUE),
    Median_noncalcified = median(x_non, na.rm = TRUE),
    Median_difference_calc_minus_non =
      median(x_calc, na.rm = TRUE) - median(x_non, na.rm = TRUE),
    Mean_calcified = mean(x_calc, na.rm = TRUE),
    Mean_noncalcified = mean(x_non, na.rm = TRUE),
    Wilcoxon_W = unname(w$statistic),
    Wilcoxon_P = w$p.value,
    Welch_t = unname(tt$statistic),
    Welch_P = tt$p.value,
    stringsAsFactors = FALSE
  )
}

primary_tests <- do.call(
  rbind,
  lapply(score_names, test_one)
)

primary_tests$FDR_BH_4scores <- p.adjust(
  primary_tests$Wilcoxon_P,
  method = "BH"
)

primary_tests$Role <- ifelse(
  primary_tests$Score == "RUNX2_SPP1",
  "Secondary calcification-associated marker score",
  "Primary frozen programme"
)

write.csv(
  primary_tests,
  file.path(OUTDIR, "05_PRIMARY_calcified_vs_noncalcified_tests.csv"),
  row.names = FALSE
)

cat("\nPRIMARY CALCIFICATION RESULTS:\n")
print(primary_tests, row.names = FALSE)

###############################################################################
# 14. EXPLORATORY BROWN-minus-BLUE AXIS
###############################################################################

axis_test <- wilcox.test(
  Brown_minus_Blue ~ Pathology_group,
  data = scores_primary,
  exact = FALSE
)

axis_summary <- data.frame(
  Analysis = "Exploratory_Brown_minus_Blue",
  N_calcified = sum(scores_primary$Pathology_group == "Calcified_FA"),
  N_noncalcified = sum(scores_primary$Pathology_group == "Noncalcified_FA"),
  Median_calcified = median(
    scores_primary$Brown_minus_Blue[
      scores_primary$Pathology_group == "Calcified_FA"
    ]
  ),
  Median_noncalcified = median(
    scores_primary$Brown_minus_Blue[
      scores_primary$Pathology_group == "Noncalcified_FA"
    ]
  ),
  Median_difference_calc_minus_non =
    median(
      scores_primary$Brown_minus_Blue[
        scores_primary$Pathology_group == "Calcified_FA"
      ]
    ) -
    median(
      scores_primary$Brown_minus_Blue[
        scores_primary$Pathology_group == "Noncalcified_FA"
      ]
    ),
  Wilcoxon_W = unname(axis_test$statistic),
  P_value = axis_test$p.value,
  stringsAsFactors = FALSE
)

write.csv(
  axis_summary,
  file.path(OUTDIR, "06_EXPLORATORY_Brown_minus_Blue_axis.csv"),
  row.names = FALSE
)
