# Phase B — empirical stability of the extreme Reitsma correlation

## Research question

Does the correlation near -1, and thus the near-zero determinant-based bivariate I², persist under changes to the estimation method, continuity correction, individual-study contributions and resampling scheme in the published TCCS ischaemic analysis?

The primary example is the eight-arm Pinto-Villalba et al. (2026) TCCS ischaemic dataset, reconstructed previously. The archived REML result from mada 0.5.12 has estimated rho(Se,FPR) approximately -0.9999999993, with positive marginal random-effects variances. It represents a boundary estimate from a fitted model, **not proof of a true population correlation of -1**.

## Unit-of-analysis caution

There are eight analysis arms and two arms with zero cells (Postert 1998 and Bar 2010). Postert 1998 and Postert 1999 are separate publications, but independence of their participant samples has NOT been verified from the original full papers.

The primary analysis follows the eight archived rows. An exploratory publication-family grouping treats both Postert rows as one cluster (seven clusters). This neither asserts actual overlap nor proves that the two observations are dependent; this remains a study-source verification task.

## Prespecified analyses

1. Primary reproduction: fit intercept-only logit(Se), logit(FPR) with REML and the baseline continuity correction (0.5 added to all arms if a zero cell exists). Validate rho and corrected Zhou–Dendukuri I² against the archive.
2. Estimator sensitivity: compare REML and maximum likelihood (ML). Both are normal-normal two-stage Reitsma models; this is NOT a comparison with a binomial-normal GLMM.
3. Continuity correction: 0.25, 0.5, 1, combined with "all" and "single" rules, under both ML and REML. Corrected Zhou internal expected variance denominators use the ORIGINAL counts, consistently with our previous audit.
4. Leave-one-arm-out: eight exclusions per estimator; examine changes in correlation, marginal variances, determinant and I².
5. Leave-one-publication-family-out: seven exclusions per estimator, including removal of both Postert rows together.
6. Exploratory nonparametric bootstrap: resample complete 2x2 rows with replacement (eight arm draws) or publication families (seven clusters, potentially including both Postert rows together). Run each resample with both REML and ML, reusing the SAME resampled dataset to compare estimators.
7. Diagnostics: failed fits, warnings, convergence codes when exposed, number of distinct units, marginal random-effects variances, rho, distance to |rho|=1, corrected Zhou I², marginal variance ratios, strongest/weakest generalized-eigenvalue fractions, and experimental determinant-of-total comparator.

## Predefined interpretation rules

The near-boundary flag is abs(rho) >= 0.999, with a more extreme secondary threshold abs(rho) >= 0.999999. These are numerical DESCRIPTORS, not formal hypothesis tests.

Reported 2.5%, 50% and 97.5% bootstrap quantiles summarize only successful refits. They are NOT validated 95% confidence intervals at the covariance boundary. Failure rates and denominators are retained explicitly. A large fraction of boundary estimates can occur with weak identifiability and does not prove the true correlation equals -1.

At k=8, resampling studies with replacement creates duplicates and sometimes few unique independent units. Such refits are diagnostic rather than reliable standalone inference; a likelihood profile and fully binomial-normal GLMM may be needed to investigate true correlation uncertainty.

## Reproduce in GitHub Actions

Workflow: Actions > R | Phase B Reitsma stability > Run workflow, after this workflow is merged into the default branch. Manual runs default to 400 bootstrap draws **for each strategy** (arm and publication family) and each estimation method (REML and ML), for a total of 1600 bootstrap refits. The pull request executes a smoke run of 30 draws per strategy and method (120 refits). User can select 0 to 2000 draws; random seed 20261015.

Local run from repo root, with mada 0.5.12 installed:

    PHASE_B_BOOTSTRAP=100 Rscript scripts/07_test_phase_b.R

The workflow uploads phase_B_outputs as a downloadable artifact, including: phaseB_study_table.csv, phaseB_baseline.csv, phaseB_correction_sensitivity.csv, phaseB_leave_one_out.csv, phaseB_bootstrap_replicates.csv, phaseB_bootstrap_summary.csv, phaseB_loo_correlation.png, phaseB_bootstrap_correlation.png and sessionInfo.txt.

## Scientific interpretation

Persistence across ML/REML, corrections and study deletion would indicate empirical stability of the NEAR-BOUNDARY FIT in these perturbations, not that the underlying population covariance is singular.

Large fluctuations under limited perturbations, nonconvergence, or a large proportion of boundary refits suggest instability or weak identification. Failure to converge is a substantive output rather than a row to silently discard.

The diagnostic is deliberately separate from the earlier mada implementation discrepancy. The original input data and archived numerical results are unchanged.

## References

- Pinto-Villalba et al. (2026): https://doi.org/10.3389/fneur.2026.1726046
- Hamza, Reitsma & Stijnen (2008): https://doi.org/10.1177/0272989X08323917
- mada 0.5.12: https://cran.r-project.org/package=mada
