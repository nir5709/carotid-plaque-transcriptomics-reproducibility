# Figure scope

This repository is designed to reproduce the computational evidence underlying the manuscript. Analysis scripts generate the numerical tables and, where useful, diagnostic or dataset-level plots needed to verify the reported results.

Final multi-panel manuscript figure assembly, typography, panel lettering, graphical layout and journal-specific artwork preparation are intentionally outside the computational reproducibility scope. Those presentation steps do not alter the underlying statistical results.

Where a manuscript panel depends on a reported statistic, the corresponding numerical source should be traceable to the relevant analysis output and to `scripts/99_validate_reference_outputs.R` where the value is part of the locked regression set.
