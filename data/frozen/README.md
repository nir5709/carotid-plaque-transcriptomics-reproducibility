# Frozen manuscript-reproduction inputs

This directory contains processed mappings, matrices, metadata and author-derived summary tables retained when exact manuscript reconstruction depends on a historical software/annotation state or when rerunning a large upstream object is unnecessary for numerical regression.

These files are **not substitutes for the original public source data**. The authoritative raw datasets remain the corresponding NCBI GEO accessions documented in `docs/datasets.md`.

## Discovery and composition

- `GSE125771_gene_level_matrix_frozen.csv` — annotated, RMA-normalized, probe-collapsed discovery expression matrix used by the GSE125771 manuscript workflows.
- `GSE125771_CIBERSORT_macrophage_myeloid_composite_correlations_frozen.csv` — author-derived historical CIBERSORT/LM22 summary used for manuscript regression.
- `GSE125771_CIBERSORT_signature_correlations_Plt0.05_samples_frozen.csv` — author-derived historical CIBERSORT/LM22 correlation summary used for manuscript regression.

CIBERSORT source code, LM22 and raw CIBERSORT fraction matrices are not redistributed.

## External bulk datasets

- `GSE21545_GPL570_probe_to_gene_frozen.csv` — final historical GPL570 probe-to-gene mapping used for GSE21545 reconstruction.
- `GSE43292_GPL6244_probe_to_gene_frozen.csv` — final historical probe-to-gene mapping used for GSE43292 reconstruction.
- `GSE163154_GPL6104_probe_to_gene_frozen.csv` — final historical probe-to-gene mapping used for GSE163154 reconstruction.
- `GSE104140_TMM_logCPM_gene_level_frozen.csv` — historical TMM/log2-CPM gene-level matrix used for deterministic GSE104140 reconstruction.
- `GSE104140_sample_metadata_frozen.csv` — frozen pathology-group metadata paired with the GSE104140 matrix.

## Single-cell datasets

- `GSE253903_pseudobulk_logCPM_frozen.csv` — donor × broad-cell-type pseudobulk logCPM matrix used by the canonical GSE253903 manuscript reconstruction.
- `GSE253903_pseudobulk_metadata_frozen.csv` — matching donor/cell-type pseudobulk metadata.
- `GSE210152_*_frozen.csv` and `GSE210152_historical_sessionInfo.txt` — historical downstream tables and environment information used by the deterministic GSE210152 reconstruction.

## Reproducibility principle

Frozen artifacts are retained only when they protect the manuscript evidence chain from documented package, annotation or clustering drift. Current-software reconstruction scripts are kept separately where practical and should not be tuned to force agreement with the frozen manuscript values.
