# 03 — Cell-composition analyses

This directory contains the author-written MCP-counter analysis used for composition-aware interpretation of GSE125771.

## MCP-counter

`01_MCPcounter_composition_adjustment.R` is part of the manuscript evidence chain. It evaluates monocytic-lineage abundance in relation to the locked transcriptional programmes and includes the documented partial-correlation sensitivity analyses.

Expression-derived composition estimates are interpreted as transcriptomic composition estimates rather than histological cell counts or direct measures of cell-intrinsic inflammasome activation.

## CIBERSORT / LM22

The manuscript also reports CIBERSORT/LM22 analyses. CIBERSORT and LM22 are third-party resources with separate terms of use. To keep this archive conservative, this repository does **not** redistribute CIBERSORT source code, packaged implementations, LM22, raw cell-fraction matrices, or executable wrappers/installers for those resources.

For manuscript reproducibility, only small author-derived historical summary tables are retained under `data/frozen/` and checked by `scripts/99_validate_reference_outputs.R`. These summaries contain the reported correlation/test results and do not reproduce the CIBERSORT implementation or LM22 matrix.

Users wishing to perform a new CIBERSORT analysis should obtain authorized access/resources directly from the provider and follow the applicable terms. See `docs/THIRD_PARTY_RESOURCES.md`.
