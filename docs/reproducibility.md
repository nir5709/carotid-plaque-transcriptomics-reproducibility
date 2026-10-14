# Reproducibility policy

This repository supports transparent reconstruction of the computational analyses and manuscript-facing numerical results.

## Repository contents

The repository includes:

- author-written analysis scripts used for reported results;
- fixed transcriptional signature definitions;
- selected derived tables needed to reproduce reported statistics;
- frozen historical mappings, matrices or intermediate tables where exact manuscript reconstruction depends on version-sensitive software or annotations;
- software/session information;
- documentation of preprocessing, annotation and statistical decisions;
- a manuscript regression validator with locked reference values and tolerances.

The repository does not intentionally include:

- raw GEO data that can be retrieved from the original accession;
- patient-identifiable or restricted information;
- local credentials, tokens, `.Renviron` or `.env` files;
- absolute machine-specific paths;
- unnecessary large intermediate objects when smaller frozen tables are sufficient;
- third-party resources whose redistribution rights are not controlled by the authors;
- final journal-layout or publication-figure assembly code.

## Historical and current-software reconstruction

Some manuscript values depend on historical software or annotation states. For these cases the canonical regression path begins from a documented frozen mapping, matrix, pseudobulk object, or derived table. Current-software reconstructions are retained separately where useful for transparency and executability, but they are not substituted for the historical manuscript path when version drift changes numerical results.

## File-path convention

Manuscript-facing scripts use paths relative to the repository root. Machine-specific absolute paths should not be introduced into release code.

## Statistical reproducibility

For reported analyses, the workflow documents or encodes the relevant biological unit, statistical test, paired/unpaired structure, multiplicity-correction family, fixed score/gene-set definition, inclusion/exclusion criteria and software/version provenance.

## Third-party resources

CIBERSORT/LM22 is the principal restricted third-party resource in this project. The repository does not redistribute the CIBERSORT implementation, LM22 matrix, raw fraction matrices or automated download/install wrappers. Only author-derived historical summary tables required for manuscript regression are retained. See [`THIRD_PARTY_RESOURCES.md`](THIRD_PARTY_RESOURCES.md).

## Release validation

Before a versioned archival release, the canonical workflows should be run from a clean clone, followed by:

```r
source("scripts/99_validate_reference_outputs.R")
```

The validator writes `results/validation/reference_output_validation.csv` and stops on a numerical mismatch. Release-environment information can then be captured with `scripts/98_capture_release_environment.R`.

Reference values and tolerances are part of the reproducibility contract and should not be changed merely to obtain a passing validation result.
