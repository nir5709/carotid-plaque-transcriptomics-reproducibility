###############################################################################
# 12. MODULE / SIGNATURE SCORES
###############################################################################

scores <- pb_meta
scores$Brown_WGCNA <- score_gene_set(logcpm, brown_genes)
scores$Blue_WGCNA <- score_gene_set(logcpm, blue_genes)
scores$Inflammasome_65 <- score_gene_set(logcpm, inflammasome_65)

write.csv(
  scores,
  file.path(OUTDIR, "03_pseudobulk_gene_set_scores.csv"),
  row.names = FALSE
)

###############################################################################
# 13. CELL-TYPE LOCALIZATION SUMMARY
###############################################################################

score_summary <- scores |>
  dplyr::group_by(Broad_label) |>
  dplyr::summarise(
    N_pseudobulks = dplyr::n(),
    N_donors = dplyr::n_distinct(Sample),
    Brown_median = median(Brown_WGCNA, na.rm = TRUE),
    Brown_mean = mean(Brown_WGCNA, na.rm = TRUE),
    Blue_median = median(Blue_WGCNA, na.rm = TRUE),
    Blue_mean = mean(Blue_WGCNA, na.rm = TRUE),
    Inflammasome65_median = median(Inflammasome_65, na.rm = TRUE),
    Inflammasome65_mean = mean(Inflammasome_65, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(N_donors >= MIN_DONORS_PER_CELLTYPE)

write.csv(
  score_summary,
  file.path(OUTDIR, "04_celltype_score_summary.csv"),
  row.names = FALSE
)

###############################################################################
# 14. PRIMARY DONOR-BALANCED TARGET-vs-OTHER TESTS
###############################################################################

delta_brown <- paired_target_delta(
  scores,
  "Brown_WGCNA",
  BROWN_PRIMARY_CELLTYPES,
  "Brown: Macrophage+Monocyte vs other cell types"
)

delta_blue <- paired_target_delta(
  scores,
  "Blue_WGCNA",
  BLUE_PRIMARY_CELLTYPES,
  "Blue: SMC+Fibroblast vs other cell types"
)

delta_inflam <- paired_target_delta(
  scores,
  "Inflammasome_65",
  BROWN_PRIMARY_CELLTYPES,
  "Inflammasome65: Macrophage+Monocyte vs other cell types"
)

deltas <- dplyr::bind_rows(delta_brown, delta_blue, delta_inflam)

write.csv(
  deltas,
  file.path(OUTDIR, "05_primary_target_vs_other_donor_deltas.csv"),
  row.names = FALSE
)

paired_tests <- do.call(
  rbind,
  lapply(
    split(
      deltas,
      interaction(deltas$Comparison, deltas$Score, drop = TRUE)
    ),
    function(x) {
      if (nrow(x) < 4L) {
        return(data.frame(
          Comparison = x$Comparison[1],
          Score = x$Score[1],
          N_donors = nrow(x),
          Median_delta = median(x$Delta_target_minus_other, na.rm = TRUE),
          Wilcoxon_V = NA_real_,
          P_value = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      tst <- suppressWarnings(
        wilcox.test(
          x$Target_mean,
          x$Other_mean,
          paired = TRUE,
          exact = FALSE
        )
      )
      data.frame(
        Comparison = x$Comparison[1],
        Score = x$Score[1],
        N_donors = nrow(x),
        Median_delta = median(x$Delta_target_minus_other, na.rm = TRUE),
        Wilcoxon_V = unname(tst$statistic),
        P_value = tst$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
)

paired_tests$FDR <- p.adjust(paired_tests$P_value, method = "BH")
write.csv(
  paired_tests,
  file.path(OUTDIR, "06_primary_target_vs_other_paired_tests.csv"),
  row.names = FALSE
)

###############################################################################
# 15. LOCALIZATION RANKING
###############################################################################

ranking_brown <- score_summary |>
  dplyr::select(Broad_label, N_donors, Median_score = Brown_median) |>
  dplyr::mutate(GeneSet = "Brown_WGCNA") |>
  dplyr::arrange(dplyr::desc(Median_score))

ranking_blue <- score_summary |>
  dplyr::select(Broad_label, N_donors, Median_score = Blue_median) |>
  dplyr::mutate(GeneSet = "Blue_WGCNA") |>
  dplyr::arrange(dplyr::desc(Median_score))

ranking_inflam <- score_summary |>
  dplyr::select(Broad_label, N_donors, Median_score = Inflammasome65_median) |>
  dplyr::mutate(GeneSet = "Inflammasome_65") |>
  dplyr::arrange(dplyr::desc(Median_score))

ranking <- dplyr::bind_rows(ranking_brown, ranking_blue, ranking_inflam)
ranking <- ranking |>
  dplyr::group_by(GeneSet) |>
  dplyr::mutate(Rank = rank(-Median_score, ties.method = "min")) |>
  dplyr::ungroup()

write.csv(
  ranking,
  file.path(OUTDIR, "07_module_localization_ranking.csv"),
  row.names = FALSE
)
