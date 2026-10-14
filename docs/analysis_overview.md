# Analysis overview

This repository follows the analytical sequence reported in the manuscript.

## 1. Discovery in GSE125771

The discovery cohort contains 40 human carotid plaque transcriptomes. Three predefined biological components are evaluated:

- **Inflammasome signature** — fixed 65-gene set represented from the 67-gene source set;
- **Fibrous remodelling signature** — `ACTA2`, `CNN1`, `MYH11`, `MYOCD`, `TAGLN`, `LMOD1`, and `ELN`;
- **Osteogenic markers** — `RUNX2` and `SPP1`.

Discovery analyses include continuous score relationships, independent high-versus-low score contrasts, differential expression and functional annotation.

## 2. WGCNA and NMF

Whole-transcriptome WGCNA identifies co-expression modules in the discovery cohort. NMF is applied to inflammatory-contrast differentially expressed genes as a complementary unsupervised analysis. Cross-method overlap is used to evaluate recurrent inflammatory and structural transcriptional structure.

## 3. Bulk cell-composition context

MCP-counter and historical CIBERSORT/LM22 analyses provide complementary expression-derived estimates of cellular composition. These results are interpreted as composition-associated transcriptomic patterns rather than histological cell counts or direct evidence of biochemical inflammasome activation.

The repository includes the author-written MCP-counter workflow. For CIBERSORT/LM22, only small author-derived historical result summaries required for manuscript regression are retained; the CIBERSORT implementation and LM22 matrix are not redistributed. See [`THIRD_PARTY_RESOURCES.md`](THIRD_PARTY_RESOURCES.md).

## 4. External bulk evaluation

Fixed discovery signatures/modules are transferred to external datasets without re-optimizing their membership:

- **GSE21545** — external fixed-signature relationships;
- **GSE43292** — paired lesion-stage comparison;
- **GSE163154** — segment-level comparison by intraplaque-haemorrhage status;
- **GSE104140** — lesion-stage and within-fibroatheroma calcification analyses.

## 5. Cell-resolved evaluation

- **GSE253903** — donor-aware pseudobulk localization across broad plaque cell types;
- **GSE210152** — macrophage/mononuclear-phagocyte state localization, marker enrichment, and donor-aware pseudobulk evaluation.

Individual cells are not treated as independent biological replicates for donor-level inference.

## 6. Evidence integration

The study integrates discovery, external bulk and cell-resolved evidence to evaluate recurrence and cellular localization of transcriptional programmes. The analyses do not establish a causal temporal sequence of plaque progression, direct biochemical inflammasome activation, or a validated quantitative measure of tissue calcium burden.

The repository targets analytical and numerical reproducibility. Final publication-figure layout and graphical assembly are outside its computational scope.
