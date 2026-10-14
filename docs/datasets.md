# Datasets

This project uses seven public human transcriptomic datasets from the NCBI Gene Expression Omnibus (GEO). Raw/source data are not duplicated in this repository.

| GEO accession | Study role | Key analytical use |
|---|---|---|
| GSE125771 | Discovery | Fixed transcriptional scores, differential expression, WGCNA, NMF, MCP-counter and historical CIBERSORT summaries |
| GSE21545 | External bulk comparison | Transfer of fixed signatures and pairwise score relationships |
| GSE43292 | Paired lesion-stage comparison | Within-patient advanced plaque versus distant macroscopically intact/earlier-stage tissue |
| GSE163154 | IPH segment comparison | Segment-level comparison of plaque segments with versus without intraplaque haemorrhage |
| GSE104140 | Lesion-stage/calcification | DIT versus advanced fibroatheroma; within-advanced-fibroatheroma calcification analysis |
| GSE253903 | Single-cell localization | Donor-aware pseudobulk localization across broad cell types |
| GSE210152 | Macrophage-state analysis | Localization of a PLIN2/TREM1-high macrophage-like state in plaque CD45+ immune cells |

## Interpretation notes

- GSE125771 and GSE21545 originate from the broader BiKE resource; participant independence between these datasets is not assumed.
- GSE43292 comparator tissue is distant macroscopically intact/earlier-stage tissue and is not treated as a healthy-control cohort.
- GSE163154 contains tissue segments. The metadata used for the manuscript did not provide verified complete patient-to-segment linkage, so the primary comparison is interpreted at the segment level and is not described as a fully patient-adjusted analysis.
- GSE104140 transcriptomic calcification status is not equated with mechanical plaque stability or vulnerability.
- GSE253903 donor × broad-cell-type pseudobulks are the primary inferential unit; individual cells are not independent biological replicates.
- GSE210152 contains CD45+ immune cells. It can evaluate inflammatory macrophage/MNP-state localization but cannot independently validate the structural smooth-muscle programme. Symptomatic/asymptomatic differences are descriptive because the dataset contains four symptomatic and two asymptomatic donors.

## Data acquisition and frozen manuscript inputs

The authoritative source for each raw dataset is its GEO accession. Where historical annotation or software versions affected exact manuscript values, the repository retains frozen mappings, processed matrices, pseudobulk inputs, or derived result tables under `data/frozen/`. These objects provide deterministic reconstruction of manuscript-facing results.

Current-software reconstructions are also provided where practical, but they are not tuned to force numerical identity with historical results when package or annotation drift occurs.
