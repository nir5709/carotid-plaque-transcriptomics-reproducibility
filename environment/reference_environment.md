# Historical analysis-environment reference

This file records package versions recovered from retained historical analysis/session material. These versions are provenance for manuscript reconstruction; they are not a claim that every workflow used one identical environment.

## Retained bulk-analysis reference

| Component | Version |
|---|---|
| R | 4.4.2 |
| NMF | 0.28 |
| GSVA | 2.0.7 |
| limma | 3.62.2 |
| msigdbr | 26.1.0 |
| GO.db | 3.20 |

Locked signature membership is stored under `signatures/`, so live MSigDB/package matching is not required for manuscript score definitions.

## GSE210152 historical environment

The recovered historical GSE210152 session used R 4.4.2 on Windows with Seurat 5.5.0, SeuratObject 5.4.0, harmony 1.2.4 and edgeR 4.4.2, together with the other packages recorded in `data/frozen/GSE210152_historical_sessionInfo.txt`.

## Historical CIBERSORT/LM22 provenance

The manuscript included a historical CIBERSORT/LM22 analysis. Because CIBERSORT/LM22 is a third-party resource with separate terms of use, the public archive does not distribute its implementation, LM22 matrix, raw cell-fraction output or an executable wrapper. The manuscript-facing correlation/test summaries from the historical run are retained as small author-derived tables under `data/frozen/` and are regression-checked directly.

A later platform rerun showed small numerical drift while preserving direction and interpretation; this is why the frozen historical summaries remain the manuscript reference. See `../docs/THIRD_PARTY_RESOURCES.md`.

Package-version matching alone is not treated as proof of equivalence. The canonical release must also pass `scripts/99_validate_reference_outputs.R` against the frozen manuscript targets.
