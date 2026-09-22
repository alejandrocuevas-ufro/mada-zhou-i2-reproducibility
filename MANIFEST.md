# Reproducibility manifest

## Inputs

- `data/published_reanalyses_2x2.csv` — clean aggregate 2×2 tables for published reanalyses.
- `data/metadta_telomerase_reconstruction.csv` — public telomerase example.
- `data/data_provenance.csv` — provenance and extraction method.

## Analysis code

- `scripts/zhou_i2_utils.R` — common I² calculations and diagnostics.
- `scripts/00_environment.R` — version check, session information, and source-code capture.
- `scripts/01_validate_public_examples.R` — AuditC6/AuditC validation.
- `scripts/02_validate_metadta.R` — independent formula comparison.
- `scripts/03_reanalyse_published.R` — published applied examples.
- `scripts/04_simulation.R` — Monte Carlo experiment and Figure 1.
- `run_all.R` — sequential execution.

## Expected generated outputs

After `Rscript run_all.R`, `results/` should contain:

- `sessionInfo.txt`
- `summary.reitsma_mada_0.5.12.R`
- `public_example_validation.csv`
- `public_example_model_diagnostics.csv`
- `metadta_formula_validation.csv`
- `published_reanalyses_summary.csv`
- `published_reanalyses_all_variants.csv`
- `published_reproduction_check.csv`
- `simulation_replicates.csv`
- `simulation_summary.csv`
- `Figure1_simulation.png`
- `Figure1_simulation.pdf`

## Audit note

- `AUDITED_CODE_SNIPPET.md` — minimal code fragment showing the audited assignments and reconstructed counterparts.
