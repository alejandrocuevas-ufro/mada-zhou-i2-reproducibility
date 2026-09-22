# Zhou–Dendukuri bivariate I² in `mada`: reproducibility files

This repository accompanies a methodological audit of the Zhou–Dendukuri bivariate I² reported by `summary.reitsma()` in R package `mada` 0.5.12.

**Archived release v1.0.0:** https://doi.org/10.5281/zenodo.22887677

The repository has four aims:

1. reproduce the I² currently returned by `mada`;
2. reconstruct the statistic using the between-study covariance matrix and reference-status denominators specified in the Zhou–Dendukuri formulation;
3. reproduce and recalculate published applications for which study-level 2×2 tables could be recovered; and
4. reproduce the Monte Carlo experiment evaluating the behaviour of the statistic as the number of studies increases while population heterogeneity remains fixed.

## Complete reproduction

The manuscript audit targets **`mada` 0.5.12**. Install that version if necessary:

```r
install.packages("remotes")
remotes::install_version("mada", version = "0.5.12", repos = "https://cloud.r-project.org")
```

From the repository root:

```bash
Rscript run_all.R
```

The workflow records `sessionInfo()` and the installed `summary.reitsma()` source, validates AuditC6/AuditC and the published `metadta` example, reproduces published applications, applies the reconstructed Zhou–Dendukuri calculation, runs 500 Monte Carlo replicates for each K = 6, 10, 20, 40, 80, and 160, and generates the simulation figure.

## Key files

- `AUDITED_CODE_SNIPPET.md` — minimal source fragment underlying the audit.
- `data/published_reanalyses_2x2.csv` — aggregate 2×2 data used in the published reanalyses.
- `data/data_provenance.csv` — provenance of each dataset.
- `scripts/zhou_i2_utils.R` — exact/native and reconstructed I² calculations with numerical boundary safeguards.
- `scripts/03_reanalyse_published.R` — reanalysis of published applications.
- `scripts/04_simulation.R` — Monte Carlo experiment and Figure 1.
- `expected/` — validated reference results.
- `results/` — selected final outputs from the validated R run.

The full replicate-level simulation table and figure files are generated deterministically by `scripts/04_simulation.R` from seed `20260918`; they need not be stored to reproduce the analysis.

## Published reanalyses

The analyses are:

- Dragota 2026: six studies of the systemic immune-inflammation index for thyroid cancer;
- Manetas-Stavrakakis 2023: 18 PPG arms and 32 single-lead ECG arms for atrial fibrillation detection;
- Pinto-Villalba 2026: 14 transcranial ultrasonography arms, including TCCS ischaemic and haemorrhagic subgroup analyses.

The script first reproduces the reported `mada` I². Reproduction is considered successful when the difference from the published one-decimal value is no greater than 0.15 percentage points. The reconstructed I² is then calculated from the same fitted Reitsma object.

## Continuity corrections

For analyses containing zero cells, the model is fitted with:

```r
correction = 0.5
correction.control = "all"
```

The native `mada` I² calculation uses `fit$freqdata`, which retains the original 2×2 counts; the audit script follows this behaviour when reproducing `summary.reitsma()`.

## Monte Carlo experiment

The simulation seed is fixed at:

```r
set.seed(20260918)
```

The population Zhou–Dendukuri I² is 90.56%. Only K changes across scenarios; the generating means, between-study covariance matrix, and diseased/non-diseased sample sizes remain fixed.

## Near-zero I² example

The Pinto-Villalba TCCS ischaemic analysis illustrates a separate property of the original Zhou–Dendukuri statistic. Its random-effects variances are non-zero, but the estimated correlation is close to −1, making the determinant of the between-study covariance matrix almost zero. This is not evidence of the `mada` implementation discrepancy itself.

## Data status

No individual participant data are included. All reanalyses use aggregate 2×2 counts available in the cited publications or reconstructed directly from published event/total information.

## Repository and archived release

GitHub: https://github.com/alejandrocuevas-ufro/mada-zhou-i2-reproducibility

Archived release v1.0.0: https://doi.org/10.5281/zenodo.22887677
