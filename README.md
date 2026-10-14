# Carotid Plaque Transcriptomics

Computational workflows and reproducibility materials supporting integrative bulk and single-cell transcriptomic analyses of human carotid atherosclerosis.

## Associated study

**Inflammasome-associated and structural transcriptional programmes in human carotid plaques: an integrative analysis of bulk and single-cell transcriptomes**

This repository provides the code, fixed gene-signature definitions, reproducibility documentation, software-environment information, and selected derived outputs used for the associated study.

## Study overview

The study integrates seven publicly available human transcriptomic datasets to examine relationships between inflammatory and smooth-muscle structural transcriptional programmes in carotid atherosclerotic plaques, their association with lesion stage and intraplaque haemorrhage, and their localization across bulk and single-cell transcriptomic data.

### Datasets

| GEO accession | Role in the study |
|---|---|
| GSE125771 | Discovery cohort: predefined signatures, differential expression, WGCNA, NMF and cell-composition analyses |
| GSE21545 | External fixed-signature comparison |
| GSE43292 | Paired lesion-stage comparison |
| GSE163154 | Segment-level comparison by intraplaque haemorrhage status |
| GSE104140 | Lesion-stage and within-fibroatheroma calcification analyses |
| GSE253903 | Donor-aware broad cell-type localization |
| GSE210152 | Macrophage-state localization |

## Transcriptional programmes and markers

The manuscript evaluates three prespecified biological components:

- **Inflammasome signature** — the fixed 65-gene set represented in the discovery dataset from the 67-gene source set. The two source genes not represented in GSE125771 were `NLRP2B` and `PYDC5`.
- **Fibrous remodelling signature** — `ACTA2`, `CNN1`, `MYH11`, `MYOCD`, `TAGLN`, `LMOD1` and `ELN`.
- **Osteogenic markers** — `RUNX2` and `SPP1`.

The locked machine-readable gene sets used by the analysis are provided in [`signatures/`](signatures/). File names are retained as computational identifiers so that the analysis scripts remain reproducible.

## Repository structure

```text
carotid-plaque-transcriptomics-reproducibility/
├── README.md
├── CITATION.cff
├── data/
│   └── frozen/
├── docs/
├── environment/
├── results/
├── scripts/
│   ├── 01_GSE125771_discovery/
│   ├── 02_WGCNA_NMF/
│   ├── 03_cell_composition/
│   ├── 04_external_validation/
│   └── 05_single_cell/
└── signatures/
```

## Reproducibility

Raw GEO data are not duplicated in this repository. Source datasets should be obtained directly from the corresponding NCBI Gene Expression Omnibus accessions listed above.

The repository contains manuscript-facing analysis workflows together with locked signature definitions, documentation, and selected frozen intermediate or derived outputs required to reproduce and verify reported results. Where exact reconstruction depends on historical software or annotation environments, frozen mappings or derived reference outputs are retained as the canonical regression path and are distinguished from current-software reruns.

The manuscript-facing workflow has undergone clean-clone reproducibility and regression testing. The locked regression validator reports **97 PASS, 0 FAIL and 0 SKIP** for the configured checks.

## Third-party resources

The study includes results produced using CIBERSORT/LM22. This repository does not redistribute CIBERSORT source code, packaged implementations, LM22, or other restricted third-party resources. See [`docs/THIRD_PARTY_RESOURCES.md`](docs/THIRD_PARTY_RESOURCES.md) for details.

## Data availability

All primary transcriptomic datasets are publicly available from the NCBI Gene Expression Omnibus under the accession numbers listed above. This repository contains analysis code and derived reproducibility materials rather than redistributed raw source data.

## Authors

Nihar G. Garg; Manan Pancholy; Samir Pancholy; Jagat S. Chauhan; Tejas M. Patel.

Affiliations and machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## Citation

Please cite the associated manuscript and the archived repository release. A persistent Zenodo DOI will be added to this repository after the first archival release.

## License

No license is asserted here for third-party data, software, reference resources, or materials governed by their original terms. A repository-level software license will be added only after the authors/institution approve the intended reuse terms.
