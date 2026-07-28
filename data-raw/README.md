# data-raw

Scripts for datasets this package can generate itself.

- `sim_decoupling.R` — simulated example, built from `sim_trends()` and
  `sim_noise_pair()` with fixed seeds. Reproducible from a clean checkout:
  `source("data-raw/sim_decoupling.R")`.

`data/morro_bay.rda` has no script here. It is field data, produced by the
analysis repository
([`lomad-analysis`](https://github.com/ruizt/lomad-analysis),
`mb-analysis/export_example.R`) and copied in as a finished `.rda`. A script
that only read one file and re-saved it would add indirection, not provenance;
the processing chain is documented in `?morro_bay` and lives with the code that
performs it.

To refresh it: run `mb-analysis/export_example.R` in the analysis repo, then
copy the resulting `morro_bay.rda` into `data/`.
