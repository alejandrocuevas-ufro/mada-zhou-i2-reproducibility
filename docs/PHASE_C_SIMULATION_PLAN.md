# Phase C: Monte Carlo design for determinant collapse with known rho

## Central question

How often does the fitted Zhou–Dendukuri bivariate I² approach zero because a fitted covariance matrix becomes nearly singular when the true random-effect covariance is **not** singular?

**Distinguish:** (i) the intrinsic determinant collapse when the population correlation is nearly ±1, (ii) fitting-boundary artifacts under moderate true correlation, and (iii) the precision and failure rate of the estimator. This simulation tests a specific bivariate normal–normal implementation in `mada` 0.5.12, not all diagnostic-test meta-analysis models.

## Generator

For study `i`, generate two independent standard normal scores `z1_i,z2_i`; set

```
eta_Se,i  = mu_Se  + sd_Se  * z1_i
eta_FPR,i = mu_FPR + sd_FPR * [rho * z1_i + sqrt(1-rho^2) * z2_i]
TP_i ~ Binomial(nD, plogis(eta_Se,i))
FP_i ~ Binomial(nND, plogis(eta_FPR,i))
FN_i = nD - TP_i
TN_i = nND - FP_i
```

For numerical stability calculate `sqrt(1-rho^2)` as `sqrt((1-|rho|)(1+|rho|))`. With `rho=±1` the true covariance is singular. Studies are independently generated; binomial counts within a study are conditionally independent given the paired latent effects.

Population parameters are identical to the earlier study of K dependence:

- `mu_Se = 2.0997`, `mu_FPR = -1.2637`
- `sd_Se = 1.1746`, `sd_FPR = 0.6381`
- `rho_true = -1, -.9999, -.99, -.9, -.5, 0, .5, .9, .99, .9999, 1`
- `K = 6, 10, 20, 40, 80, 160`
- `nD/nND = 50/50` (small) or `150/500` (original)
- 132 independent scenarios; **100 replication datasets/scenario for the full run** (13,200 primary REML fits), with additional ML fits for `K=6,20,80` at `rho=0,+/-.99` in both sample schemes (1,800 more ML fits).
- Seed `20261015`. PR smoke: 3 replications/scenario.

## Estimands and endpoints

The **population** Zhou index follows the Zhou–Dendukuri formula using known marginal between-study variances and expected within-study logit variances in the corresponding disease-status groups. Note that for fixed positive marginal variances and sample sizes the population index converges to zero when `|rho_true|→1`. This is a property of the target parameter, not an estimation error.

The **fitted** Zhou index uses the previously reconstructed parameters from `fit$Psi`, with the corrected reference-status denominators, **not** the package's native implementation that uses fixed-effect SEs. The estimator is REML by default, using continuity correction +0.5 with correction.control="all", matching the audit protocol.

For each independent replication record:

1. Status (successful, fitting error, extraction error) and warning count; failures remain denominators, not silently discarded.
2. `rho_hat`, `tau2_Se_hat`, `tau2_FPR_hat`, `I2_hat`; marginal heterogeneity diagnostic fractions and generalized-eigenvalue directions.
3. Fitted boundary `|rho_hat|>=0.999` and extreme boundary `>=0.999999` (descriptive thresholds, not tests).
4. Low fitted `I2_hat<0.05` (under 5%).
5. **Discordance indicator**: fitted `I2_hat < 0.05` **when population** `I2_ZD >= 0.25` (at least 25%). This targets substantial negative estimation discrepancy, not a formal false-positive rate or clinically validated heterogeneity classification.
6. Mean bias and absolute error of fitted I² with respect to the **population** target, conditional on success.
7. Scenario success fraction, zero-cell fraction and descriptive 2.5%, 50% and 97.5% empirical quantiles.

Evaluate discordance especially under true `|rho|<=0.9`; display population and fitted I² side-by-side. For `|rho_true|≈1`, a small population index is expected and not an error in estimating the chosen target. Also examine `K` dependence and sample-size precision, not only overall pooled averages. Conditional-on-success estimates can be biased by fit failures; report attempts and successes.

ML sensitivity uses exactly the same simulated counts as REML for each paired scenario-replicate. The REML and ML fits of a matched replicate are therefore paired, not independent simulation runs.

## Reproduction

Run `Rscript scripts/09_test_phase_c.R` after installing `mada` version 0.5.12. Use `PHASE_C_REPS=100` for the full run. GitHub Actions workflow `.github/workflows/r-phase-c.yml` runs 3 replicates/scenario on PR and 100/scenario on pushes to `main` (or manually configurable via Actions UI).

Output files in `phase_C_outputs/`:

- `phaseC_design_population.csv`: exact population inputs/targets.
- `phaseC_replicate_results.csv`: fitted and failed model-level records.
- `phaseC_scenario_summary.csv`: scenario-level metrics and denominator counts.
- `phaseC_nonextreme_discordance.csv`: key subset to evaluate estimation-induced low I².
- `phaseC_boundary_frequency.png`: boundary probability by `rho`, `K` and sample size.
- `phaseC_discordant_I2.png`: low-I² discordance probability.
- `sessionInfo.txt`: computation environment, seed and package version.

## Limitations and scientific claims

This uses a correctly specified bivariate **normal–normal approximation** for fitting, while count data arise from a binomial-normal hierarchical generative process. With low denominators and zero cells, discretisation and continuity correction may affect behaviour. The results do not replace direct fitting of the underlying binomial-normal GLMM. There is no adjustment for multiple test thresholds because the endpoints are descriptive, but comparisons should not be interpreted as formal significance tests. The design fixes marginal between-study variances, means and per-study denominators, and thus does not characterize all realistic test-accuracy meta-analyses. Expansion to alternative true variance configurations and a direct binomial-normal model belongs to a later simulation extension.

## Background

Zhou Y & Dendukuri N (2014), *Statistics in Medicine*, DOI 10.1002/sim.6115. Original `mada` archive and code audit in this repository.
