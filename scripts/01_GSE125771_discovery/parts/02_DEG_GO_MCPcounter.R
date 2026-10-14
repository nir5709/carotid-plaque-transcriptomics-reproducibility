###############################################################################
# 2. INDEPENDENT TOP10-vs-BOTTOM10 DIFFERENTIAL EXPRESSION
###############################################################################

run_extreme_limma <- function(signature_name) {
  grouping <- group_tables[[signature_name]]
  selected <- grouping[grouping$Group %in% c("Low", "High"), , drop = FALSE]
  selected <- selected[order(selected$Group), , drop = FALSE]
  group_factor <- factor(as.character(selected$Group), levels = c("Low", "High"))
  design <- model.matrix(~ 0 + group_factor)
  colnames(design) <- levels(group_factor)
  fit <- limma::lmFit(expr_full[, selected$Sample, drop = FALSE], design)
  fit <- limma::contrasts.fit(
    fit,
    limma::makeContrasts(High_vs_Low = High - Low, levels = design)
  )
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tab <- limma::topTable(
    fit,
    coef = "High_vs_Low",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  tab$Gene <- rownames(tab)
  tab$Status <- "Not_significant"
  tab$Status[tab$adj.P.Val < FDR_CUTOFF & tab$logFC >= LOGFC_CUTOFF] <- "Up"
  tab$Status[tab$adj.P.Val < FDR_CUTOFF & tab$logFC <= -LOGFC_CUTOFF] <- "Down"
  list(all = tab, significant = tab[tab$Status != "Not_significant", , drop = FALSE])
}

extreme_results <- setNames(lapply(signatures, run_extreme_limma), signatures)

deg_summary <- do.call(rbind, lapply(signatures, function(s) {
  x <- extreme_results[[s]]$significant
  dir.create(file.path(OUTDIR, "03_DEG", s), recursive = TRUE, showWarnings = FALSE)
  write.csv(extreme_results[[s]]$all,
            file.path(OUTDIR, "03_DEG", s, "all_genes_limma.csv"), row.names = FALSE)
  write.csv(x,
            file.path(OUTDIR, "03_DEG", s, "significant_DEGs.csv"), row.names = FALSE)
  data.frame(
    Signature = s,
    Significant_total = nrow(x),
    Upregulated = sum(x$Status == "Up"),
    Downregulated = sum(x$Status == "Down"),
    stringsAsFactors = FALSE
  )
}))
write.csv(deg_summary, file.path(OUTDIR, "03_DEG", "DEG_summary.csv"), row.names = FALSE)

###############################################################################
# 3. DIRECTION-AWARE DEG OVERLAP
###############################################################################

status_table <- function(s, nm) {
  x <- extreme_results[[s]]$all[, c("Gene", "Status")]
  names(x)[2] <- nm
  x
}

overlap <- Reduce(
  function(x,y) merge(x,y,by="Gene",all=TRUE),
  list(
    status_table("Inflammasome", "Inflammasome_status"),
    status_table("Fibrous", "Fibrous_status"),
    status_table("Osteogenic", "Osteogenic_status")
  )
)
for (nm in names(overlap)[-1]) overlap[[nm]][is.na(overlap[[nm]])] <- "Not_significant"
overlap$Significant_in_all_three <-
  overlap$Inflammasome_status != "Not_significant" &
  overlap$Fibrous_status != "Not_significant" &
  overlap$Osteogenic_status != "Not_significant"
overlap$Inflam_up_Osteogenic_up_Fibrous_down <-
  overlap$Inflammasome_status == "Up" &
  overlap$Osteogenic_status == "Up" &
  overlap$Fibrous_status == "Down"
write.csv(overlap, file.path(OUTDIR, "03_DEG", "direction_aware_three_signature_overlap.csv"), row.names = FALSE)

###############################################################################
# 4. GO BIOLOGICAL PROCESS OVER-REPRESENTATION
###############################################################################

map_to_entrez <- function(genes) {
  z <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(genes),
      keytype = "SYMBOL",
      columns = c("SYMBOL","ENTREZID")
    )
  )
  unique(na.omit(z$ENTREZID))
}

run_go <- function(genes, universe) {
  q <- map_to_entrez(genes)
  u <- map_to_entrez(universe)
  if (length(q) < 5L) return(data.frame())
  ans <- tryCatch(
    clusterProfiler::enrichGO(
      gene = q,
      universe = u,
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
  if (is.null(ans)) data.frame() else as.data.frame(ans)
}

for (s in signatures) {
  sig <- extreme_results[[s]]$significant
  up <- run_go(sig$Gene[sig$Status == "Up"], rownames(expr_full))
  down <- run_go(sig$Gene[sig$Status == "Down"], rownames(expr_full))
  write.csv(up, file.path(OUTDIR, "04_GO_BP", paste0(s, "_upregulated_GO_BP.csv")), row.names = FALSE)
  write.csv(down, file.path(OUTDIR, "04_GO_BP", paste0(s, "_downregulated_GO_BP.csv")), row.names = FALSE)
}

###############################################################################
# 5. MCP-COUNTER CELL-COMPOSITION CONTEXT
###############################################################################

mcp_matrix <- NULL
if (requireNamespace("MCPcounter", quietly = TRUE)) {
  mcp_matrix <- MCPcounter::MCPcounter.estimate(expr_full, featuresType = "HUGO_symbols")
} else {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  try(remotes::install_github("ebecht/MCPcounter", subdir = "Source", upgrade = "never"), silent = TRUE)
  if (requireNamespace("MCPcounter", quietly = TRUE)) {
    mcp_matrix <- MCPcounter::MCPcounter.estimate(expr_full, featuresType = "HUGO_symbols")
  }
}

if (!is.null(mcp_matrix)) {
  write.csv(data.frame(Cell_type=rownames(mcp_matrix), mcp_matrix, check.names=FALSE),
            file.path(OUTDIR, "05_MCPcounter", "MCPcounter_scores.csv"), row.names = FALSE)

  mcp_cor <- do.call(rbind, lapply(rownames(mcp_matrix), function(cell) {
    do.call(rbind, lapply(signatures, function(sig) {
      ct <- safe_cor(as.numeric(mcp_matrix[cell, score_df$Sample]), score_df[[sig]])
      data.frame(Cell_type=cell, Signature=sig, Spearman_rho=ct$rho,
                 P_value=ct$P_value, N=ct$N, stringsAsFactors=FALSE)
    }))
  }))
  mcp_cor$FDR_within_signature <- ave(
    mcp_cor$P_value, mcp_cor$Signature,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  write.csv(mcp_cor, file.path(OUTDIR, "05_MCPcounter", "MCPcounter_signature_correlations.csv"), row.names = FALSE)

  monocytic_row <- grep("monocytic|monocyte|macrophage", rownames(mcp_matrix), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(monocytic_row) && "PTPRC" %in% rownames(expr_full) && requireNamespace("ppcor", quietly = TRUE)) {
    common <- intersect(score_df$Sample, colnames(mcp_matrix))
    pcor_res <- ppcor::pcor.test(
      score_df$Inflammasome[match(common, score_df$Sample)],
      as.numeric(mcp_matrix[monocytic_row, common]),
      cbind(
        as.numeric(expr_full["PTPRC", common]),
        colMeans(mcp_matrix[, common, drop=FALSE], na.rm=TRUE)
      ),
      method = "spearman"
    )
    partial_table <- data.frame(
      Cell_type = monocytic_row,
      Partial_Spearman_rho = unname(pcor_res$estimate),
      P_value = pcor_res$p.value,
      N = pcor_res$n,
      Controlled_for = "PTPRC expression and mean MCP-counter score",
      stringsAsFactors = FALSE
    )
    write.csv(partial_table,
              file.path(OUTDIR, "05_MCPcounter", "MCPcounter_monocytic_partial_correlation.csv"),
              row.names = FALSE)
  }
}
