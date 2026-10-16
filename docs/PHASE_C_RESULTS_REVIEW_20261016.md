# Phase C — Review of the completed known-rho Monte Carlo

**Review date:** 16 October 2026 (UTC). **Data:** GitHub Actions [run 38837447182](https://github.com/alejandrocuevas-ufro/mada-zhou-i2-reproducibility/actions/runs/38837447182) and its [full numeric artifact](https://github.com/alejandrocuevas-ufro/mada-zhou-i2-reproducibility/actions/runs/38837447182/artifacts/12020836327). **Code:** pinned `mada` 0.5.12, seed 20261015.

## Computational and consistency checks

- 132 prespecified population scenarios (11 correlations, 6 study counts, 2 precision regimes).
- 100 generated datasets per scenario and 13,200 REML fits; 1,800 paired ML fits in a targeted subset; all **15,000** fits produced finite model results, with **0 captured warnings** and no captured fit or extraction errors.

> **Correction (Phase D1, 2026-10-16).** The preceding sentence must not be read
> as a quality check. `status == "ok"` means only that `reitsma()` returned an
> object, and `mada`/`mvmeta` emit no warning when the fitted correlation
> reaches the singular boundary — which is why 0 warnings were recorded across
> 15,000 fits even in scenarios where 79 % of fits sit exactly at |ρ̂| = 1. The
> Phase D1 audit additionally shows `converged == TRUE` in 100 % of replicates
> in every scenario examined, boundary solutions included. Any report of these
> runs should state the fraction of boundary solutions next to the success
> count. See [`PHASE_D1_AUDIT.md`](PHASE_D1_AUDIT.md) §3.
- Independently checked all **150** method-specific scenario summaries against their 100 underlying replicate records. Fitted I² medians and key event probabilities agree exactly.
- Independently recomputed all 132 theoretical population Zhou–Dendukuri I² targets with the original means, random-effects variances, and within-study formulas; all agreed numerically.
- Each targeted ML replicate has the corresponding REML replicate on the same scenario and replicate ID.
- **Qualification:** the existing script's `status=ok` verifies that a model object and valid derived quantities were returned, but does *not yet independently inspect all low-level optimizer convergence codes or Hessian identifiability*. Absence of warnings is not equivalent to proof of reliable inference.

## Key observed comparisons

All percentages below use the corrected Zhou–Dendukuri reconstruction based on `fit$Psi` and the reference-status sample denominators. The two regimes are 50 cases / 50 non-cases versus 150 cases / 500 non-cases.

| True rho | K | Precision | Population I² | Fraction with estimated I² <5% |
|---:|---:|---|---:|---:|
| 0 | 6 | 50/50 | 77.11% | 41/100 |
| 0 | 10 | 50/50 | 77.11% | 13/100 |
| 0 | 20 | 50/50 | 77.11% | 0/100 |
| 0 | 6 | 150/500 | 94.86% | 3/100 |
| 0 | 20 | 150/500 | 94.86% | 0/100 |
| +0.99 | 6 | 50/50 | 32.21% | 87/100 |
| +0.99 | 160 | 50/50 | 32.21% | 88/100 |
| +0.99 | 6 | 150/500 | 72.24% | 68/100 |
| +0.99 | 160 | 150/500 | 72.24% | 6/100 |

The `rho_true=0`, `K=6`, 50/50 result is a direct counterexample to interpreting the near-zero fitted statistic as proof of a nearly singular *population* between-study covariance. Low precision, small K, and true correlation close to the boundary substantially affect how often such errors occur.

The persistent 88% value for `rho_true=+0.99`, `K=160`, 50/50 warrants a *model-based explanation*, not an assumption of an algorithm defect. Study-level binomial counts are generated from a binomial-normal hierarchical model but analysed by a two-stage normal-normal approximation with continuity correction. In finite samples per study, this approximation may differ materially from the data-generating likelihood even with many studies.

## Visualization defect found and corrected

**The original two PNG charts in the completed run must not be used in a manuscript.** Their plotting function used `subset(sums, method=="reml" & scheme==scheme)`. Under R's non-standard evaluation, `scheme==scheme` compared the dataframe's column to itself, so **both panels mixed both precision regimes**. It affected plot filtering and titles, NOT the Monte Carlo generated records or CSV summary calculations.

The permanent R correction uses explicit column indexing: `sums[sums$method=="reml" & sums$scheme==scheme, , drop=FALSE]` and asserts that each panel has the expected 66 scenario results with a single sample-size regime. It also displays **NA rather than 0** for the severe-discordance diagnostic where population I² is below 25%, because the discordance event is inapplicable there.

Corrected plots may be regenerated from the archived scenario CSV without redoing any statistical estimation. After integrating this code correction, future GitHub Actions reruns will generate the corrected PNGs.

## Interpretation and scientific status

**Supported:** determinant-based heterogeneity is highly sensitive to near-boundary correlations. In the specified small-K, low-precision binomial-normal generating scenarios, fitted near-singular covariance estimates and near-zero corrected Zhou I² can appear even if the population rho is zero and the corresponding population I² is large. Adding more studies and increasing within-study precision generally decreases such discordance when rho is moderate, but the behaviour near true |rho|=0.99 is more complex.

**Not established:** whether these results generalize to binomial-normal GLMM fitting, other levels of marginal heterogeneity, or all distributions of study sizes. The `I²<5%` and `population I²≥25%` thresholds are descriptive; the frequencies are NOT type-I error or clinically calibrated false-negative rates. With 100 replications per scenario, uncertainty is substantial (for a 41% observed fraction, its Monte Carlo standard error is approximately 4.9 percentage points); publication-level precision calls for at least 500–1000 replications per priority scenario. Numerical convergence and parameter identifiability need separate checks.

## Reproducible next steps

1. Expand priority scenarios to at least 1000 replications with stable random-number streams and uncertainty intervals for Monte Carlo proportions.
2. Record optimizer convergence and Hessian/curvature diagnostics where exposed by `mada` / `mvmeta`, and distinguish valid fit-object return from reliable model convergence.
3. Fit a binomial-normal hierarchical model directly to the 2×2 counts for key cases, contrasting it with the normal-normal `mada` approximation.
4. Report the joint estimated `rho`, marginal random-effects variances, both generalized eigenvalues, and Zhou I² rather than Zhou I² alone.

**Conclusion:** the numeric results are internally consistent and promising, but the *original full-run figures were incorrect* and have been corrected in code. The simulation provides methodological evidence, not yet a fully externally validated general estimator-performance study.
