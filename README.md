# Zhou–Dendukuri bivariate I² in `mada`: reproducibility files

This repository accompanies a methodological audit of the Zhou–Dendukuri bivariate I² reported by `summary.reitsma()` in R package `mada` 0.5.12.

The repository has four aims:

1. reproduce the I² currently returned by `mada`;
2. reconstruct the statistic using the between-study covariance matrix and reference-status denominators specified in the Zhou–Dendukuri formulation;
3. reproduce and recalculate published applications for which study-level 2×2 tables could be recovered; and
4. reproduce the Monte Carlo experiment evaluating the behaviour of the statistic as the number of studies increases while population heterogeneity remains fixed.

## Repository structure

```text
.
├── README.md
├── README_ES.md
├── run_all.R
├── data/
│   ├── published_reanalyses_2x2.csv
│   ├── metadta_telomerase_reconstruction.csv
│   ├── data_provenance.csv
│   └── DATA_DICTIONARY.md
├── scripts/
│   ├── 00_environment.R
│   ├── 01_validate_public_examples.R
│   ├── 02_validate_metadta.R
│   ├── 03_reanalyse_published.R
│   ├── 04_simulation.R
│   └── zhou_i2_utils.R
├── expected/
│   ├── published_i2_reference.csv
│   └── simulation_summary_reference.csv
└── results/
```

## Software

The manuscript audit targets **`mada` 0.5.12**. The scripts deliberately stop if another `mada` version is installed, because a later package version may change the code under study.

A convenient way to install the audited version is:

```r
install.packages("remotes")
remotes::install_version("mada", version = "0.5.12", repos = "https://cloud.r-project.org")
```

The simulation also uses `MASS`, which is a recommended R package in standard R installations.

## Complete reproduction

From the repository root, run:

```bash
Rscript run_all.R
```

The run will:

- record `sessionInfo()` and the exact installed `summary.reitsma()` source;
- reproduce the AuditC6 and AuditC calculations;
- reproduce the algebra used in the published `metadta` telomerase example;
- fit the published 2×2 datasets and reproduce the I² values printed in the source articles;
- calculate the two intermediate corrections and the full Zhou–Dendukuri reconstruction;
- run 500 Monte Carlo replicates for each K = 6, 10, 20, 40, 80, and 160; and
- generate the simulation figure in PNG and PDF formats.

Outputs are written to `results/`.

## Published reanalyses

`data/published_reanalyses_2x2.csv` contains only public, study-level aggregate counts used in the manuscript reanalyses. Provenance is documented in `data/data_provenance.csv`.

The analyses are:

- Dragota 2026: six studies of the systemic immune-inflammation index for thyroid cancer;
- Manetas-Stavrakakis 2023: 18 PPG arms and 32 single-lead ECG arms for atrial fibrillation detection;
- Pinto-Villalba 2026: 14 transcranial ultrasonography arms, including the TCCS ischaemic and haemorrhagic subgroup analyses.

The script first reproduces the reported `mada` I². Reproduction is considered successful when the difference from the published one-decimal value is no greater than 0.15 percentage points. The corrected I² is then calculated from the same fitted Reitsma object.

## Continuity corrections

For analyses containing zero cells, the model is fitted with:

```r
correction = 0.5
correction.control = "all"
```

This reproduces the relevant published analyses. The native `mada` I² calculation uses `fit$freqdata`, which retains the original 2×2 counts; the audit script follows this behaviour when reproducing `summary.reitsma()`.

## Monte Carlo experiment

The simulation seed is fixed at:

```r
set.seed(20260918)
```

Population parameters are stored directly in `scripts/04_simulation.R`. The population Zhou–Dendukuri I² is approximately 90.56%. Only K changes across scenarios; the generating means, between-study covariance matrix, and diseased/non-diseased sample sizes remain fixed.

`expected/simulation_summary_reference.csv` records the final reference values from the validated R run. The archived `results/` directory contains the corresponding complete outputs.

## Interpretation of the near-zero I² example

The Pinto-Villalba TCCS ischaemic analysis is retained because it illustrates a separate property of the original Zhou–Dendukuri statistic. Its random-effects variances are non-zero, but the estimated correlation is close to −1, making the determinant of the between-study covariance matrix almost zero. This example should not be interpreted as evidence of the `mada` implementation discrepancy itself.

## Data status

No individual participant data are included. All reanalyses use aggregate 2×2 counts available in the cited publications or reconstructed directly from published event/total information. See `DATA_NOTICE.md` and `data/data_provenance.csv`.

## Citation and archival version

Before manuscript submission, create a versioned public release (for example, GitHub + Zenodo) and replace the repository placeholder in the manuscript with the permanent DOI. The archived release should include the `results/sessionInfo.txt` generated by the final R run.

## Audited source fragment

`AUDITED_CODE_SNIPPET.md` records the minimal `summary.reitsma()` fragment relevant to the two implementation discrepancies and the corresponding reconstruction used in this repository.
