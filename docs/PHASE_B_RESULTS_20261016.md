# Phase B results — empirical stability and profile of the fitted correlation

**Validated:** 2026-10-16 (UTC); GitHub Actions, R release, mada **0.5.12**.
**Case:** Pinto-Villalba et al. (2026), TCCS ischaemic subgroup, eight archived analysis arms.
**Analysis code:** scripts/07_phase_b_empirical_stability.R, scripts/07_test_phase_b.R and scripts/08_phase_b_profile_likelihood.R.
**Seed:** 20261015; **400 bootstrap draws per arm/family strategy AND per ML/REML estimator**.
**Validation run (all checks passed):** https://github.com/alejandrocuevas-ufro/mada-zhou-i2-reproducibility/actions/runs/38835153868
**Results artifact:** https://github.com/alejandrocuevas-ufro/mada-zhou-i2-reproducibility/actions/runs/38835153868/artifacts/12019505101

The archived dataset and prior reproduction results have NOT been modified.

## 1. Baseline model replication

| Fit | tau2 Se | tau2 FPR | rho(Se,FPR), rounded | corrected Zhou I2 |
|---|---:|---:|---:|---:|
| REML, all, +0.5 | 0.10265573 | 2.442952 | about -1 | 0.00256359% |
| ML, all, +0.5 | 0.08006207 | 2.068790 | about -1 | 0.00034558% |

Both model fits succeeded with no warnings in the baseline. The REML solution reproduces the previously archived correlation and corrected I2 to numerical tolerance. The raw Pearson correlation among the eight study-level logit sensitivity and logit FPR pairs, with 0.5 added to all cells, is about -0.250; that observable descriptive quantity is NOT the latent between-study rho.

## 2. Sensitivity analyses

With 0.25, 0.5 or 1 as continuity correction, applied to all arms or only arms with a zero cell, there were 12 successful ML/REML sensitivity fits. Estimated rho remained nearly -1 (the plain console print rounds it to -1).

Corrected Zhou I2 remained tiny:
- REML: 0.001941% to 0.010878% across all corrections/modes.
- ML: 0.000346% to 0.002851% across all corrections/modes.

There were 16 leave-one-arm-out fits and 14 leave-one-publication-family-out fits (including combined deletion of Postert 1998 and Postert 1999), all successful. Estimated correlations remained near -1, and corrected Zhou I2 ranged from roughly 0.0000032% to 0.015224% across these exclusions. These findings show stability of the POINT ESTIMATE to these chosen perturbations, not strong identification of a population correlation.

The two Postert papers originated from the same Bochum research group, but participant-level independence/overlap has not been verified. Treating them as one publication family is a sensitivity check, not an assertion that they analysed overlapping patients.

## 3. Nonparametric bootstrap (full validation, B=400 per cell)

Threshold for a near-boundary estimate: absolute fitted rho >= 0.999. Fractions below refer to successful refits. Quantiles of the index are in **percentage points of I2**, NOT proportions.

| Resampling unit | Estimator | Successful / attempted | Near-boundary | Median corrected I2 | 97.5% quantile corrected I2 |
|---|---|---:|---:|---:|---:|
| Eight arms | REML | 400 / 400 | 380 / 400 = 95.00% | 0.002915% | 40.4055% |
| Eight arms | ML | 400 / 400 | 382 / 400 = 95.50% | 0.001677% | 30.4344% |
| Seven publication families | REML | 400 / 400 | 387 / 400 = 96.75% | 0.004314% | 0.05544% |
| Seven publication families | ML | 400 / 400 | 391 / 400 = 97.75% | 0.002123% | 0.02325% |

The stricter near-perfect threshold |rho| >= 0.999999 occurred, respectively, in 93.0%, 93.5%, 94.25%, and 95.25% of successful refits.

The 97.5% percentile of rho itself was around **+0.999997** for arm bootstrap under either estimator, revealing a tail of sign-reversed, extreme correlation fits. The 97.5% percentile of rho was roughly -0.0397 (REML) and -0.0513 (ML) under family resampling. These tails deserve study and must not be masked by reporting only median rho near -1.

**Critical caveat:** resampling only eight arms with replacement leaves a median of FIVE distinct arms per resample. These bootstrap distributions are small-sample diagnostics, not calibrated 95% confidence intervals, and have visible boundary masses and asymmetric tails. All refits completing is not the same as reliable parameter inference. Do not use these percentiles to make claims about coverage, type-I error or the true rho.

## 4. Fixed-rho profile likelihood

For each fixed rho, the two marginal between-study variances were optimised and both pooled means were analytically profiled from the multivariate normal-normal likelihood using the study-specific transformed within-study variances. This matches the Gaussian likelihood geometry of the mada 0.5.12 Reitsma fit, up to terms constant in rho; optimized values were checked against the package fit.

| Fixed rho | Relative ML deviance | Relative REML deviance |
|---:|---:|---:|
| -1.00 | 0.0000 | 0.0000 |
| -0.95 | 0.1408 | 0.1118 |
| -0.80 | 0.5160 | 0.4310 |
| -0.50 | 1.1587 | 1.0440 |
| 0.00 | 1.8872 | 2.0253 |
| +0.50 | 1.8880 | 2.0505 |
| +1.00 | 1.8889 | 2.0512 |

**Key observation:** a maximum at the negative correlation boundary coexists with relatively small relative deviance for markedly less negative rho. Under ML, when rho is fixed at 0, the optimised between-study sensitivity variance decreases to approximately 9.4e-8; when rho is fixed positive, it reaches the optimisation lower bound (about 1.52e-8). The correlation is not well defined as a marginal variance tends to zero.

**Do NOT apply a standard chi-square likelihood ratio cutoff here as a calibrated 95% confidence bound.** The covariance parameter space and singular boundaries violate the usual regular interior assumptions. These relative differences diagnose a potential identification problem; a rigorously calibrated inferential conclusion will need dedicated methodology, such as a parametric-bootstrap calibration or an alternative model.

## 5. Synthesis and what the evidence does (and does not) support

These empirical analyses support two complementary claims:

1. The near-singular covariance and near-zero corrected Zhou–Dendukuri index are stable as **fitted outputs** of the two-stage Reitsma model across ML/REML, correction choices and leave-one-out perturbations.
2. This stability should not be conflated with precision of the latent correlation: the near-boundary maximum, the relatively flat likelihood profile, alternative solutions with small sensitivity variance, and bootstrap tails point to weak or nonregular identification.

They do NOT establish that the true random-effects correlation equals -1. They do NOT establish that the fitted between-study variances represent accurately estimated population heterogeneity, nor that an alternative proposed index is a calibrated replacement.

The original implementation discrepancy in mada and the structural determinant problem analysed in Phase A remain separate issues.

## 6. Follow-on priorities

- Source verification for possible patient/cohort overlap between Postert 1998 and 1999 (PMIDs 9587989 and 10471430), and scrutiny of ultrasound-examination rather than participant counting.
- Fit a **binomial-normal** bivariate hierarchical model using the raw counts as a methodological sensitivity; do not present this as the same likelihood as the existing normal-normal Reitsma approximation.
- Calibrate constrained rho inference or use parametric simulations under multiple true correlations and study sizes (Phase C).
- Evaluate the eigenvalue-based diagnostic alongside the original Zhou index, with uncertainty and effective-rank behavior.

## Sources

- Pinto-Villalba et al. (2026): https://doi.org/10.3389/fneur.2026.1726046
- Postert et al. (1998): https://pubmed.ncbi.nlm.nih.gov/9587989/
- Postert et al. (1999): https://pubmed.ncbi.nlm.nih.gov/10471430/
- Hamza, Reitsma & Stijnen (2008): https://doi.org/10.1177/0272989X08323917
