# 05 — Single-cell analyses

These workflows provide cell-resolved evaluation while preserving donor-level inference. Individual cells are not treated as independent biological replicates for donor-level tests.

## GSE253903

Two complementary routes are retained.

### Canonical manuscript reconstruction

`GSE253903/03_stage2_frozen_reproduction.R` starts from the frozen donor × broad-cell-type pseudobulk logCPM matrix and metadata under `data/frozen/`. This is the deterministic manuscript path and reproduces the score, donor-delta, paired-test, ranking and coverage outputs without the historical multi-gigabyte Seurat object.

### Current/public-data reconstruction

`GSE253903/01_stage1_QC_clustering.R` followed by `GSE253903/02_stage2_celltype_localization.R` reconstructs the workflow from public single-cell data through donor-aware clustering, locked broad-cell-type annotation, pseudobulk formation and downstream validation. `00_repair_download.R` is only a download-recovery utility.

Primary inference uses donor × broad-cell-type pseudobulks. Brown is evaluated in macrophage/monocyte populations, Blue in SMC/fibroblast populations, and Inflammasome-65 across broad cell types.

## GSE210152

GSE210152 contains plaque CD45+ immune cells. It can localize inflammatory macrophage/MNP states but cannot independently validate the Blue structural/SMC programme.

### Canonical manuscript reconstruction

Run:

```bash
Rscript scripts/05_single_cell/GSE210152/03_frozen_historical_reproduction.R
```

The driver stages frozen historical downstream tables from `data/frozen/` and applies the documented `g5` → `5` candidate-cluster label normalization via `02_g5_cluster_label_repair.R`. It reconstructs all seven fixed manuscript tables exactly at the validated tolerance without committing or rerunning the historical multi-gigabyte Seurat object.

### Current/public-H5AD reconstruction

`GSE210152/01_macrophage_state_analysis.R` reconstructs the public H5AD workflow and includes compatibility updates for contemporary Seurat/Harmony releases. Current package versions can produce a different MNP clustering partition from the historical Windows run. This route is therefore retained as a current-software reconstruction and must not replace the frozen historical manuscript outputs.

The candidate PLIN2/TREM1 inflammatory-LAM-like state is interpreted using marker expression, enrichment and donor representation together; it is not treated as one-to-one identity with the Brown WGCNA module. Symptomatic/asymptomatic differences are descriptive only because the dataset contains four symptomatic and two asymptomatic donors.
