# Clean rerun protocol

This protocol is intended for release verification from a fresh clone of the repository.

## 1. Clone and enter the repository

```bash
git clone https://github.com/nir5709/carotid-plaque-transcriptomics-reproducibility.git
cd carotid-plaque-transcriptomics-reproducibility
```

## 2. Confirm frozen inputs

Before running the canonical manuscript reconstruction, confirm that the required files documented in `data/frozen/README.md` are present. Raw GEO source data are not required for the canonical frozen reconstruction paths where a frozen processed input is explicitly provided.

## 3. Prepare the R environment

```bash
Rscript scripts/00_install_dependencies.R
```

Package installation can vary by operating system. The retained historical versions are documented in `environment/reference_environment.md`; exact manuscript regression is protected by locked inputs and numerical validation rather than by assuming that package-version matching alone guarantees equivalence.

## 4. Run the manuscript workflows

Run from the repository root in this order:

```bash
Rscript scripts/01_GSE125771_discovery/01_discovery_core.R
Rscript scripts/02_WGCNA_NMF/01_whole_transcriptome_WGCNA.R
Rscript scripts/03_cell_composition/01_MCPcounter_composition_adjustment.R
Rscript scripts/04_external_validation/GSE21545_fixed_signature_validation.R
Rscript scripts/04_external_validation/GSE43292_paired_stage_validation.R
Rscript scripts/04_external_validation/GSE163154_IPH_segment_validation.R
Rscript scripts/04_external_validation/GSE104140_lesion_calcification_validation.R
Rscript scripts/05_single_cell/GSE253903/03_stage2_frozen_reproduction.R
Rscript scripts/05_single_cell/GSE210152/03_frozen_historical_reproduction.R
```

The current/public-data single-cell reconstruction routes are documented separately in `scripts/05_single_cell/README.md` and are not required for the canonical frozen regression path.

## 5. Validate manuscript numerical outputs

```bash
Rscript scripts/99_validate_reference_outputs.R
```

A release candidate is numerically acceptable only when all configured checks pass and no required check is skipped. The validator writes:

`results/validation/reference_output_validation.csv`

Expected values and tolerances are locked to the manuscript analysis and should not be altered to force a passing result.

## 6. Capture release environment

After successful validation:

```bash
Rscript scripts/98_capture_release_environment.R
```

This records non-identifying operating-system information, R/Bioconductor versions, package versions, the Git commit and a copy/summary of the validator output under `environment/`.

## 7. Final release checks

Before creating a versioned GitHub/Zenodo release, confirm that no credentials, patient-identifiable information, machine-specific absolute paths, raw restricted third-party resources or unintended large intermediate objects have been added.
