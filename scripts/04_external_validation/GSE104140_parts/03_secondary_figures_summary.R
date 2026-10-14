###############################################################################
# 15. SECONDARY CONTEXT — DIT vs ADVANCED FIBROATHEROMA
###############################################################################

secondary_gsm <- meta$GSM[
  meta$Pathology_group %in% c("DIT", "Calcified_FA", "Noncalcified_FA")
]

expr_secondary <- expr_gene[, secondary_gsm, drop = FALSE]
scores_secondary <- score_frozen_sets(expr_secondary, GENE_SETS)

scores_secondary <- merge(
  meta[, c("GSM", "Title", "Pathology_group")],
  scores_secondary,
  by = "GSM",
  sort = FALSE
)

scores_secondary <- scores_secondary[
  match(secondary_gsm, scores_secondary$GSM),
  ,
  drop = FALSE
]

scores_secondary$Stage2 <- ifelse(
  scores_secondary$Pathology_group == "DIT",
  "DIT",
  "Advanced_FA"
)

secondary_tests <- do.call(
  rbind,
  lapply(score_names, function(sc) {
    x_adv <- scores_secondary[scores_secondary$Stage2 == "Advanced_FA", sc]
    x_dit <- scores_secondary[scores_secondary$Stage2 == "DIT", sc]

    w <- wilcox.test(x_adv, x_dit, exact = FALSE)

    data.frame(
      Score = sc,
      N_advanced = length(x_adv),
      N_DIT = length(x_dit),
      Median_advanced = median(x_adv, na.rm = TRUE),
      Median_DIT = median(x_dit, na.rm = TRUE),
      Median_difference_advanced_minus_DIT =
        median(x_adv, na.rm = TRUE) - median(x_dit, na.rm = TRUE),
      P_value = w$p.value,
      stringsAsFactors = FALSE
    )
  })
)

secondary_tests$FDR_BH_4scores <- p.adjust(
  secondary_tests$P_value,
  method = "BH"
)

write.csv(
  scores_secondary,
  file.path(OUTDIR, "08_SECONDARY_27sample_DIT_vs_Advanced_programme_scores.csv"),
  row.names = FALSE
)

write.csv(
  secondary_tests,
  file.path(OUTDIR, "07_SECONDARY_advanced_vs_DIT_tests.csv"),
  row.names = FALSE
)

###############################################################################
# 16. FIGURE
###############################################################################

plot_df <- scores_primary %>%
  dplyr::select(
    GSM,
    Pathology_group,
    Brown_WGCNA,
    Blue_WGCNA,
    Inflammasome_65,
    RUNX2_SPP1
  ) %>%
  pivot_longer(
    cols = c(Brown_WGCNA, Blue_WGCNA, Inflammasome_65, RUNX2_SPP1),
    names_to = "Programme",
    values_to = "Score"
  )

plot_df$Programme <- factor(
  plot_df$Programme,
  levels = c(
    "Brown_WGCNA",
    "Blue_WGCNA",
    "Inflammasome_65",
    "RUNX2_SPP1"
  ),
  labels = c(
    "Brown WGCNA",
    "Blue WGCNA",
    "Inflammasome-65",
    "RUNX2-SPP1"
  )
)

p <- ggplot(
  plot_df,
  aes(x = Pathology_group, y = Score)
) +
  geom_boxplot(width = 0.60, outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  facet_wrap(~ Programme, scales = "free_y", nrow = 1) +
  labs(
    x = NULL,
    y = "Frozen programme score (gene-wise z-score mean)",
    title = "GSE104140: calcified vs non-calcified advanced fibroatheroma"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave(
  file.path(OUTDIR, "Figure_GSE104140_calcified_vs_noncalcified.png"),
  p,
  width = 12,
  height = 4.5,
  dpi = 400
)

ggsave(
  file.path(OUTDIR, "Figure_GSE104140_calcified_vs_noncalcified.pdf"),
  p,
  width = 12,
  height = 4.5
)

###############################################################################
# 17. FINAL SUMMARY
###############################################################################

summary_file <- file.path(
  OUTDIR,
  "FINAL_GSE104140_CALCIFICATION_SUMMARY.txt"
)

sink(summary_file)

cat("GSE104140 CALCIFICATION VALIDATION — VERSION 4\n")
cat("================================================\n\n")

cat("DATA HANDLING\\n")
cat(
  "Canonical manuscript reproduction uses the frozen historical ",
  "GSE104140 TMM/log2-CPM gene-level expression matrix.\\n",
  sep = ""
)
cat(
  "The frozen matrix was generated during the original analysis from ",
  "public NCBI GEO RNA-seq counts after group-blind expression filtering, ",
  "TMM normalization and Entrez-to-HGNC annotation.\\n",
  sep = ""
)
cat(
  "The canonical rerun does not re-query contemporary gene annotation ",
  "resources because those resources can change across releases.\\n",
  sep = ""
)
cat(
  "Historical 32-sample pathology metadata are also frozen so that ",
  "sample-group assignments are deterministic.\\n\\n",
  sep = ""
)

cat("CURRENT GEO PATHOLOGY GROUPS\n")
print(table(meta$Pathology_group))

cat("\nGENE-SET COVERAGE\n")
print(coverage, row.names = FALSE)

cat("\nPRIMARY: CALCIFIED vs NON-CALCIFIED ADVANCED FIBROATHEROMA\n")
print(primary_tests, row.names = FALSE)

cat("\nEXPLORATORY BROWN-minus-BLUE AXIS\n")
print(axis_summary, row.names = FALSE)

cat("\nSECONDARY: ADVANCED FIBROATHEROMA vs DIT\n")
print(secondary_tests, row.names = FALSE)

cat("\nINTERPRETATION RULES\n")
cat(
  "1. Primary comparison = 11 calcified vs 8 non-calcified advanced ",
  "fibroatheromas using current GEO labels.\n",
  sep = ""
)
cat("2. Brown, Blue and Inflammasome-65 are frozen programmes.\n")
cat(
  "3. RUNX2-SPP1 is secondary and is not a complete osteogenic pathway.\n"
)
cat(
  "4. Calcified must not be equated automatically with stable/safe, ",
  "and non-calcified must not be equated with vulnerable/dangerous.\n",
  sep = ""
)
cat(
  "5. GSE104140 contains pathology-defined calcification, not continuous ",
  "CT/Agatston calcium scores.\n",
  sep = ""
)
cat(
  "6. A null primary result is retained; do not redefine gene sets after ",
  "seeing the result.\n",
  sep = ""
)
cat(
  "7. DIT vs advanced fibroatheroma is secondary context only.\n"
)

sink()

###############################################################################
# 18. SESSION INFO
###############################################################################

capture.output(
  sessionInfo(),
  file = file.path(OUTDIR, "sessionInfo.txt")
)

###############################################################################
# 19. FINISH
###############################################################################

cat("\n============================================================\n")
cat("GSE104140 CALCIFICATION VALIDATION COMPLETE\n")
cat("============================================================\n\n")

cat("Main outputs:\n")
cat("1. FINAL_GSE104140_CALCIFICATION_SUMMARY.txt\n")
cat("2. 03_gene_set_coverage.csv\n")
cat("3. 04_PRIMARY_19sample_frozen_programme_scores.csv\n")
cat("4. 05_PRIMARY_calcified_vs_noncalcified_tests.csv\n")
cat("5. 06_EXPLORATORY_Brown_minus_Blue_axis.csv\n")
cat("6. 07_SECONDARY_advanced_vs_DIT_tests.csv\n")
cat("7. 08_SECONDARY_27sample_DIT_vs_Advanced_programme_scores.csv\n")
cat("8. Figure_GSE104140_calcified_vs_noncalcified.png\n")
