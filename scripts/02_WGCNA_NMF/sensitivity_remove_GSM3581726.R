###############################################################################
# IDEA 2 FINAL SENSITIVITY CHECK
# REMOVE POTENTIAL OUTLIER GSM3581726
###############################################################################

load(
  "IDEA2_GSE125771_WHOLE_TRANSCRIPTOME_WGCNA/08_Reproducibility/IDEA2_WGCNA_core_objects.RData"
)

outlier_sample <- "GSM3581726"

###############################################################################
# 1. REMOVE ONLY FROM SENSITIVITY DATA
###############################################################################

keep_samples <- rownames(datExpr) != outlier_sample

datExpr_sens <- datExpr[
  keep_samples,
  ,
  drop = FALSE
]

datTraits_sens <- datTraits[
  keep_samples,
  ,
  drop = FALSE
]

cat("Primary samples:", nrow(datExpr), "\n")
cat("Sensitivity samples:", nrow(datExpr_sens), "\n")

###############################################################################
# 2. KEEP THE ORIGINAL WGCNA MODULE ASSIGNMENTS
# Recalculate eigengenes with 39 samples
###############################################################################

MEs_sens <- WGCNA::moduleEigengenes(
  datExpr_sens,
  colors = moduleColors,
  impute = TRUE,
  excludeGrey = FALSE
)$eigengenes

MEs_sens <- WGCNA::orderMEs(
  MEs_sens
)

###############################################################################
# 3. REPEAT MODULE-TRAIT ASSOCIATIONS
###############################################################################

sens_trait <- WGCNA::bicorAndPvalue(
  MEs_sens,
  datTraits_sens,
  maxPOutliers = 0.05
)

sens_cor <- sens_trait$bicor
sens_p   <- sens_trait$p

###############################################################################
# 4. BROWN-BLUE DIRECT CORRELATION
###############################################################################

brown_blue_sens <- WGCNA::bicorAndPvalue(
  MEs_sens[, "MEbrown", drop = FALSE],
  MEs_sens[, "MEblue", drop = FALSE],
  maxPOutliers = 0.05
)

cat("\nBrown vs Blue after removing GSM3581726:\n")
print(brown_blue_sens$bicor)
print(brown_blue_sens$p)

###############################################################################
# 5. EXTRACT THE MAIN ASSOCIATIONS
###############################################################################

main_results_sens <- data.frame(
  Relationship = c(
    "Brown vs Inflammasome",
    "Brown vs Fibrous",
    "Brown vs Osteogenic",
    "Blue vs Inflammasome",
    "Blue vs Fibrous",
    "Blue vs Osteogenic",
    "Brown vs Blue"
  ),
  Correlation = c(
    sens_cor["MEbrown", "Inflammasome"],
    sens_cor["MEbrown", "Fibrous"],
    sens_cor["MEbrown", "Osteogenic"],
    sens_cor["MEblue", "Inflammasome"],
    sens_cor["MEblue", "Fibrous"],
    sens_cor["MEblue", "Osteogenic"],
    brown_blue_sens$bicor[1, 1]
  ),
  P_value = c(
    sens_p["MEbrown", "Inflammasome"],
    sens_p["MEbrown", "Fibrous"],
    sens_p["MEbrown", "Osteogenic"],
    sens_p["MEblue", "Inflammasome"],
    sens_p["MEblue", "Fibrous"],
    sens_p["MEblue", "Osteogenic"],
    brown_blue_sens$p[1, 1]
  )
)

main_results_sens$FDR <- p.adjust(
  main_results_sens$P_value,
  method = "BH"
)

print(main_results_sens)

write.csv(
  main_results_sens,
  "IDEA2_GSE125771_WHOLE_TRANSCRIPTOME_WGCNA/Outlier_sensitivity_GSM3581726.csv",
  row.names = FALSE
)
