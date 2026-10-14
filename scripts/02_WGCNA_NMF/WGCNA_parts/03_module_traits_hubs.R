###############################################################################
# 15. RECOMPUTE AND ORDER MODULE EIGENGENES
###############################################################################

MEs <- WGCNA::moduleEigengenes(
  datExpr,
  colors = moduleColors,
  impute = TRUE,
  excludeGrey = FALSE
)$eigengenes

MEs <- WGCNA::orderMEs(MEs)

###############################################################################
# 16. MODULE SIZE TABLE AND ASSIGNMENTS
###############################################################################

module_size_table <- as.data.frame(table(moduleColors), stringsAsFactors = FALSE)
colnames(module_size_table) <- c("Module", "Genes")
module_size_table <- module_size_table[
  order(module_size_table$Genes, decreasing = TRUE),
  , drop = FALSE
]

write_csv_safe(
  module_size_table,
  file.path(OUTDIR, "03_WGCNA_modules", "Module_sizes.csv")
)

module_assignment <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

write_csv_safe(
  module_assignment,
  file.path(OUTDIR, "03_WGCNA_modules", "Whole_transcriptome_module_assignments.csv")
)

non_grey_modules <- sort(unique(moduleColors[moduleColors != "grey"]))
cat("\nNumber of non-grey modules:\n")
print(length(non_grey_modules))

###############################################################################
# 17. BLOCK DENDROGRAMS + MODULE COLORS
###############################################################################

if (length(net$dendrograms) > 0) {
  for (b in seq_along(net$dendrograms)) {
    block_genes <- net$blockGenes[[b]]
    pdf(
      file.path(
        OUTDIR,
        "03_WGCNA_modules",
        paste0("Block_", b, "_gene_dendrogram_modules.pdf")
      ),
      width = 12,
      height = 7
    )
    WGCNA::plotDendroAndColors(
      net$dendrograms[[b]],
      moduleColors[block_genes],
      "Module colors",
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05,
      main = paste("WGCNA block", b)
    )
    dev.off()
  }
}

###############################################################################
# 18. MODULE EIGENGENE CORRELATION
###############################################################################

ME_cor <- WGCNA::bicor(MEs, maxPOutliers = MAX_P_OUTLIERS)
write_csv_safe(
  ME_cor,
  file.path(OUTDIR, "03_WGCNA_modules", "Module_eigengene_correlations.csv")
)

report_MEs <- colnames(MEs)[colnames(MEs) != "MEgrey"]

if (length(report_MEs) >= 2) {
  eig_cor_report <- ME_cor[report_MEs, report_MEs, drop = FALSE]
  pdf(
    file.path(OUTDIR, "03_WGCNA_modules", "Module_eigengene_correlation_heatmap.pdf"),
    width = max(7, 0.4 * length(report_MEs) + 3),
    height = max(6, 0.4 * length(report_MEs) + 2)
  )
  WGCNA::labeledHeatmap(
    Matrix = eig_cor_report,
    xLabels = sub("^ME", "", colnames(eig_cor_report)),
    yLabels = sub("^ME", "", rownames(eig_cor_report)),
    ySymbols = sub("^ME", "", rownames(eig_cor_report)),
    colorLabels = FALSE,
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = matrix(
      sprintf("%.2f", eig_cor_report),
      nrow = nrow(eig_cor_report),
      dimnames = dimnames(eig_cor_report)
    ),
    zlim = c(-1, 1),
    main = "Module eigengene correlations"
  )
  dev.off()
}

###############################################################################
# 19. MODULE-TRAIT CORRELATIONS
###############################################################################

module_trait_test <- WGCNA::bicorAndPvalue(
  MEs,
  datTraits,
  maxPOutliers = MAX_P_OUTLIERS
)

moduleTraitCor <- module_trait_test$bicor
moduleTraitP <- module_trait_test$p
moduleTraitFDR <- matrix(
  p.adjust(as.vector(moduleTraitP), method = "BH"),
  nrow = nrow(moduleTraitP),
  ncol = ncol(moduleTraitP),
  dimnames = dimnames(moduleTraitP)
)

write_csv_safe(
  moduleTraitCor,
  file.path(OUTDIR, "04_Module_trait", "Module_trait_correlations.csv")
)
write_csv_safe(
  moduleTraitP,
  file.path(OUTDIR, "04_Module_trait", "Module_trait_raw_P.csv")
)
write_csv_safe(
  moduleTraitFDR,
  file.path(OUTDIR, "04_Module_trait", "Module_trait_BH_FDR.csv")
)

###############################################################################
# 20. LONG-FORM MODULE-TRAIT TABLE
###############################################################################

module_trait_long <- expand.grid(
  ModuleEigengene = rownames(moduleTraitCor),
  Trait = colnames(moduleTraitCor),
  stringsAsFactors = FALSE
)
module_trait_long$Module <- sub("^ME", "", module_trait_long$ModuleEigengene)
module_trait_long$Spearman_like_bicor <- as.vector(moduleTraitCor)
module_trait_long$P_value <- as.vector(moduleTraitP)
module_trait_long$FDR <- as.vector(moduleTraitFDR)
module_trait_long$Abs_correlation <- abs(module_trait_long$Spearman_like_bicor)
module_trait_long <- module_trait_long[
  order(module_trait_long$Trait, -module_trait_long$Abs_correlation),
  , drop = FALSE
]

write_csv_safe(
  module_trait_long,
  file.path(OUTDIR, "04_Module_trait", "Module_trait_results_long.csv")
)

###############################################################################
# 21. MODULE-TRAIT HEATMAP
###############################################################################

heat_idx <- rownames(moduleTraitCor) != "MEgrey"
heat_cor <- moduleTraitCor[heat_idx, , drop = FALSE]
heat_fdr <- moduleTraitFDR[heat_idx, , drop = FALSE]

text_matrix <- matrix(
  paste0(
    sprintf("%.2f", heat_cor),
    "\nFDR=",
    format.pval(heat_fdr, digits = 2, eps = 1e-99)
  ),
  nrow = nrow(heat_cor),
  dimnames = dimnames(heat_cor)
)

pdf(
  file.path(OUTDIR, "04_Module_trait", "Module_trait_heatmap.pdf"),
  width = 8,
  height = max(6, 0.45 * nrow(heat_cor) + 2)
)
WGCNA::labeledHeatmap(
  Matrix = heat_cor,
  xLabels = colnames(heat_cor),
  yLabels = sub("^ME", "", rownames(heat_cor)),
  ySymbols = sub("^ME", "", rownames(heat_cor)),
  colorLabels = FALSE,
  colors = WGCNA::blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Unsupervised WGCNA modules vs predefined plaque scores"
)
dev.off()

###############################################################################
# 22. STRONGEST MODULE FOR EACH TRAIT
###############################################################################

non_grey_trait_results <- module_trait_long[
  module_trait_long$Module != "grey",
  , drop = FALSE
]

top_module_by_trait <- do.call(
  rbind,
  lapply(unique(non_grey_trait_results$Trait), function(trait_name) {
    x <- non_grey_trait_results[
      non_grey_trait_results$Trait == trait_name,
      , drop = FALSE
    ]
    x[which.max(x$Abs_correlation), , drop = FALSE]
  })
)
rownames(top_module_by_trait) <- NULL

write_csv_safe(
  top_module_by_trait,
  file.path(OUTDIR, "04_Module_trait", "Strongest_module_for_each_trait.csv")
)

significant_module_trait <- non_grey_trait_results[
  is.finite(non_grey_trait_results$FDR) &
    non_grey_trait_results$FDR < 0.05,
  , drop = FALSE
]

write_csv_safe(
  significant_module_trait,
  file.path(OUTDIR, "04_Module_trait", "FDR_significant_module_trait_associations.csv")
)

###############################################################################
# 23. GENE -> MODULE EIGENGENE CORRELATIONS (kME)
###############################################################################

mm_test <- WGCNA::bicorAndPvalue(
  datExpr,
  MEs,
  maxPOutliers = MAX_P_OUTLIERS
)
MM <- mm_test$bicor
MM_P <- mm_test$p

###############################################################################
# 24. GENE -> TRAIT CORRELATIONS (GENE SIGNIFICANCE)
###############################################################################

gs_test <- WGCNA::bicorAndPvalue(
  datExpr,
  datTraits,
  maxPOutliers = MAX_P_OUTLIERS
)
GS <- gs_test$bicor
GS_P <- gs_test$p

###############################################################################
# 25. EXTRACT EACH GENE'S MEMBERSHIP IN ITS OWN MODULE
###############################################################################

own_me_name <- paste0("ME", moduleColors)
own_me_index <- match(own_me_name, colnames(MEs))
own_kME <- rep(NA_real_, ncol(datExpr))
own_kME_P <- rep(NA_real_, ncol(datExpr))
valid_me <- !is.na(own_me_index)
gene_index <- seq_len(ncol(datExpr))

own_kME[valid_me] <- MM[cbind(gene_index[valid_me], own_me_index[valid_me])]
own_kME_P[valid_me] <- MM_P[cbind(gene_index[valid_me], own_me_index[valid_me])]

gene_results <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  kME_own = own_kME,
  kME_own_P = own_kME_P,
  stringsAsFactors = FALSE
)

###############################################################################
# 26. ADD GENE SIGNIFICANCE FOR EACH TRAIT
###############################################################################

for (trait_name in colnames(datTraits)) {
  gene_results[[paste0("GS_", trait_name)]] <- GS[, trait_name]
  gene_results[[paste0("GS_P_", trait_name)]] <- GS_P[, trait_name]
  gene_results[[paste0("GS_FDR_", trait_name)]] <- p.adjust(
    GS_P[, trait_name],
    method = "BH"
  )
}

write_csv_safe(
  gene_results,
  file.path(
    OUTDIR,
    "05_Hub_genes",
    "All_genes_module_membership_and_gene_significance.csv"
  )
)

###############################################################################
# 27. TOP HUB GENES FOR EVERY NON-GREY MODULE
###############################################################################

hub_tables <- lapply(non_grey_modules, function(mod) {
  x <- gene_results[gene_results$Module == mod, , drop = FALSE]
  x <- x[order(abs(x$kME_own), decreasing = TRUE), , drop = FALSE]
  if (nrow(x) > TOP_HUBS_PER_MODULE) {
    x <- x[seq_len(TOP_HUBS_PER_MODULE), , drop = FALSE]
  }
  x$Hub_rank <- seq_len(nrow(x))
  x
})

top_hubs <- do.call(rbind, hub_tables)
rownames(top_hubs) <- NULL

write_csv_safe(
  top_hubs,
  file.path(OUTDIR, "05_Hub_genes", "Top_hub_genes_per_module.csv")
)

###############################################################################
# 28. HUB TABLE FOR TRAIT-ASSOCIATED MODULES
###############################################################################

interesting_modules <- unique(c(
  top_module_by_trait$Module,
  significant_module_trait$Module
))
interesting_modules <- interesting_modules[interesting_modules != "grey"]
interesting_hubs <- top_hubs[top_hubs$Module %in% interesting_modules, , drop = FALSE]

write_csv_safe(
  interesting_hubs,
  file.path(OUTDIR, "05_Hub_genes", "Hub_genes_in_trait_associated_modules.csv")
)

###############################################################################
# 29. MODULE MEMBERSHIP vs GENE SIGNIFICANCE PLOTS
###############################################################################

for (trait_name in unique(top_module_by_trait$Trait)) {
  mod <- top_module_by_trait$Module[top_module_by_trait$Trait == trait_name][1]
  plot_df <- gene_results[gene_results$Module == mod, , drop = FALSE]
  gs_col <- paste0("GS_", trait_name)
  plot_df$GS_abs_for_plot <- abs(plot_df[[gs_col]])

  p <- ggplot(
    plot_df,
    aes(x = abs(kME_own), y = GS_abs_for_plot)
  ) +
    geom_point(alpha = 0.65) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Module ", mod, ": module membership vs gene significance"),
      subtitle = paste0("Trait: ", trait_name),
      x = "|Module membership (kME)|",
      y = paste0("|Gene significance for ", trait_name, "|")
    )

  stem <- paste0(
    "MM_vs_GS_", safe_file(trait_name), "_module_", safe_file(mod)
  )
  ggsave(
    file.path(OUTDIR, "05_Hub_genes", paste0(stem, ".pdf")),
    p, width = 6, height = 5
  )
  ggsave(
    file.path(OUTDIR, "05_Hub_genes", paste0(stem, ".png")),
    p, width = 6, height = 5, dpi = 300
  )
}
