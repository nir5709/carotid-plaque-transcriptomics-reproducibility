###############################################################################
# 30. MAP WGCNA GENES TO ENTREZ IDs
###############################################################################

all_symbols <- unique(colnames(datExpr))

symbol_map <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = all_symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
)

symbol_map <- symbol_map[
  !is.na(symbol_map$SYMBOL) & !is.na(symbol_map$ENTREZID),
  , drop = FALSE
]
symbol_map <- unique(symbol_map)
universe_entrez <- unique(symbol_map$ENTREZID)

###############################################################################
# 31. ENRICH EACH NON-GREY MODULE
###############################################################################

go_module_list <- list()

for (mod in non_grey_modules) {
  module_symbols <- colnames(datExpr)[moduleColors == mod]
  module_entrez <- unique(
    symbol_map$ENTREZID[symbol_map$SYMBOL %in% module_symbols]
  )

  if (length(module_entrez) < 10) next

  ego <- tryCatch(
    clusterProfiler::enrichGO(
      gene = module_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(ego)) next
  ego_df <- as.data.frame(ego)
  if (nrow(ego_df) == 0) next

  ego_df$Module <- mod
  ego_df$Module_genes_total <- length(module_symbols)
  ego_df$Module_genes_with_Entrez <- length(module_entrez)
  go_module_list[[mod]] <- ego_df
}

if (length(go_module_list) > 0) {
  go_all <- do.call(rbind, go_module_list)
  rownames(go_all) <- NULL

  write_csv_safe(
    go_all,
    file.path(OUTDIR, "06_GO_BP", "GO_BP_all_modules.csv")
  )

  go_sig <- go_all[
    is.finite(go_all$p.adjust) & go_all$p.adjust < GO_FDR_CUTOFF,
    , drop = FALSE
  ]

  write_csv_safe(
    go_sig,
    file.path(OUTDIR, "06_GO_BP", "GO_BP_FDR_lt_0.05.csv")
  )

  #############################################################################
  # 32. GO PLOTS FOR THE STRONGEST MODULE FOR EACH TRAIT
  #############################################################################

  for (trait_name in unique(top_module_by_trait$Trait)) {
    mod <- top_module_by_trait$Module[top_module_by_trait$Trait == trait_name][1]
    go_mod <- go_sig[go_sig$Module == mod, , drop = FALSE]
    if (nrow(go_mod) == 0) next

    go_mod <- go_mod[order(go_mod$p.adjust), , drop = FALSE]
    go_plot_df <- head(go_mod, 10)
    go_plot_df$minus_log10_FDR <- -log10(
      pmax(go_plot_df$p.adjust, .Machine$double.xmin)
    )

    p_go <- ggplot(
      go_plot_df,
      aes(
        x = minus_log10_FDR,
        y = reorder(Description, minus_log10_FDR),
        size = Count
      )
    ) +
      geom_point() +
      theme_classic(base_size = 11) +
      labs(
        title = paste0("GO-BP: module ", mod),
        subtitle = paste0("Strongest WGCNA module for ", trait_name),
        x = "-log10(BH FDR)",
        y = NULL
      )

    stem <- paste0(
      "GO_BP_", safe_file(trait_name), "_module_", safe_file(mod)
    )
    ggsave(
      file.path(OUTDIR, "06_GO_BP", paste0(stem, ".pdf")),
      p_go, width = 8, height = 5.5
    )
    ggsave(
      file.path(OUTDIR, "06_GO_BP", paste0(stem, ".png")),
      p_go, width = 8, height = 5.5, dpi = 300
    )
  }
} else {
  go_all <- data.frame()
  go_sig <- data.frame()
  warning("No GO-BP enrichment results were generated.")
}

###############################################################################
# OPTIONAL POST-HOC OVERLAP WITH THE PREVIOUS NMF MODULES
###############################################################################

if (file.exists(NMF_FILE)) {
  nmf_modules <- read.csv(
    NMF_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_nmf_cols <- c("Gene", "Module")

  if (all(required_nmf_cols %in% colnames(nmf_modules))) {
    nmf_modules$Gene <- clean_gene(nmf_modules$Gene)
    nmf_modules <- nmf_modules[
      nmf_modules$Gene %in% colnames(datExpr),
      , drop = FALSE
    ]
    nmf_module_names <- sort(unique(nmf_modules$Module))
    universe_genes <- colnames(datExpr)
    overlap_list <- list()
    counter <- 1L

    for (wgcna_mod in non_grey_modules) {
      wgcna_genes <- universe_genes[moduleColors == wgcna_mod]

      for (nmf_mod in nmf_module_names) {
        nmf_genes <- unique(nmf_modules$Gene[nmf_modules$Module == nmf_mod])
        a <- length(intersect(wgcna_genes, nmf_genes))
        b <- length(setdiff(wgcna_genes, nmf_genes))
        c <- length(setdiff(nmf_genes, wgcna_genes))
        d <- length(universe_genes) - a - b - c

        contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
        fisher_result <- fisher.test(contingency, alternative = "greater")

        overlap_list[[counter]] <- data.frame(
          WGCNA_Module = wgcna_mod,
          NMF_Module = nmf_mod,
          WGCNA_Module_Size = length(wgcna_genes),
          NMF_Genes_Present = length(nmf_genes),
          Overlap_Genes = a,
          Odds_Ratio = unname(fisher_result$estimate),
          P_value = fisher_result$p.value,
          Overlap_gene_names = paste(
            sort(intersect(wgcna_genes, nmf_genes)),
            collapse = ";"
          ),
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }

    nmf_overlap <- do.call(rbind, overlap_list)
    nmf_overlap$FDR <- p.adjust(nmf_overlap$P_value, method = "BH")
    nmf_overlap <- nmf_overlap[
      order(nmf_overlap$FDR, -nmf_overlap$Overlap_Genes),
      , drop = FALSE
    ]

    write_csv_safe(
      nmf_overlap,
      file.path(
        OUTDIR,
        "07_NMF_overlap",
        "WGCNA_vs_previous_NMF_module_overlap.csv"
      )
    )
  } else {
    warning(
      paste(
        "NMF file exists but does not contain columns:",
        paste(required_nmf_cols, collapse = ", ")
      )
    )
  }
} else {
  message("NMF_gene_modules.csv not found. Optional NMF overlap analysis skipped.")
}

###############################################################################
# 33. SAVE CORE R OBJECTS
###############################################################################

save(
  expr_full,
  datExpr,
  datTraits,
  sft,
  softPower,
  net,
  moduleColors,
  MEs,
  moduleTraitCor,
  moduleTraitP,
  moduleTraitFDR,
  module_trait_long,
  gene_results,
  top_hubs,
  top_module_by_trait,
  significant_module_trait,
  file = file.path(
    OUTDIR,
    "08_Reproducibility",
    "IDEA2_WGCNA_core_objects.RData"
  )
)

###############################################################################
# 34. SESSION INFORMATION
###############################################################################

capture.output(
  sessionInfo(),
  file = file.path(OUTDIR, "08_Reproducibility", "sessionInfo.txt")
)

###############################################################################
# 35. FINAL SUMMARY
###############################################################################

summary_file <- file.path(OUTDIR, "FINAL_IDEA2_SUMMARY.txt")

sink(summary_file)
cat("IDEA 2: WHOLE-TRANSCRIPTOME UNSUPERVISED WGCNA\n")
cat("================================================\n\n")
cat("Expression genes used:", ncol(datExpr), "\n")
cat("Samples used:", nrow(datExpr), "\n")
cat("Network type:", NETWORK_TYPE, "\n")
cat("Correlation:", COR_TYPE, "\n")
cat("Soft-threshold power:", softPower, "\n")
cat("Power selection reason:", power_reason, "\n")
cat("Non-grey modules:", length(non_grey_modules), "\n\n")
cat("MODULE SIZES\n")
print(module_size_table)
cat("\n\nSTRONGEST MODULE FOR EACH TRAIT\n")
print(top_module_by_trait)
cat("\n\nFDR-SIGNIFICANT MODULE-TRAIT ASSOCIATIONS\n")
if (nrow(significant_module_trait) > 0) {
  print(significant_module_trait)
} else {
  cat("None at BH FDR < 0.05.\n")
}
cat("\n\nTOP HUB GENES IN TRAIT-ASSOCIATED MODULES\n")
print(interesting_hubs)
if (exists("nmf_overlap")) {
  cat("\n\nWGCNA vs PREVIOUS NMF OVERLAP -- TOP RESULTS\n")
  print(head(nmf_overlap, 20))
}
sink()

###############################################################################
# 36. FINISHED
###############################################################################

cat("\n============================================================\n")
cat("IDEA 2 WGCNA ANALYSIS COMPLETED\n")
cat("============================================================\n\n")
cat("CHECK THESE FILES FIRST:\n\n")
cat("1. ", file.path(OUTDIR, "02_Soft_threshold", "Soft_threshold_results.csv"), "\n", sep = "")
cat("2. ", file.path(OUTDIR, "03_WGCNA_modules", "Module_sizes.csv"), "\n", sep = "")
cat("3. ", file.path(OUTDIR, "04_Module_trait", "Module_trait_results_long.csv"), "\n", sep = "")
cat("4. ", file.path(OUTDIR, "04_Module_trait", "Strongest_module_for_each_trait.csv"), "\n", sep = "")
cat("5. ", file.path(OUTDIR, "05_Hub_genes", "Hub_genes_in_trait_associated_modules.csv"), "\n", sep = "")
cat("6. ", file.path(OUTDIR, "06_GO_BP", "GO_BP_FDR_lt_0.05.csv"), "\n", sep = "")
cat("7. ", summary_file, "\n", sep = "")
cat("\nDo NOT interpret the modules until the module-trait table, hub genes, GO enrichment, and QC diagnostics have been reviewed together.\n")
