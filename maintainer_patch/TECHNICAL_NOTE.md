# Technical note: Zhou–Dendukuri bivariate I² in mada 0.5.12

## Summary

The issue concerns the post-estimation Zhou–Dendukuri I² reported by `summary.reitsma()`. The Reitsma model fit itself is not changed.

Two differences were identified between the published Zhou–Dendukuri formulation and the code in `mada` 0.5.12.

### 1. Between-study variance terms

The current implementation uses the squared standard errors of the pooled fixed-effect coefficients:

```r
tau_a <- (tabfixed[1,2])^2
tau_b <- (tabfixed[2,2])^2
cor_sq <- (corRandom[1,2])^2
tau <- (1-cor_sq)*tau_a*tau_b
```

These are estimates of `Var(mu_hat)`, not the between-study variances. The fitted Reitsma object already contains the between-study variance-covariance matrix in `Psi`.

For the Zhou–Dendukuri determinant, the proposed patch uses:

```r
tau_a <- Psi[1,1]
tau_b <- Psi[2,2]
tau <- det(Psi)
```

The determinant is clipped to zero only when a negative value is within a floating-point tolerance compatible with a positive-semidefinite boundary.

### 2. Reference-status denominators

The current implementation uses:

```r
Nia <- 1/(object$freqdata$TP + object$freqdata$FP)
Nib <- 1/(object$freqdata$TN + object$freqdata$FN)
```

Sensitivity is defined among participants with the reference condition and specificity among those without it. The corresponding denominators are therefore:

```r
Nia <- 1/(object$freqdata$TP + object$freqdata$FN)
Nib <- 1/(object$freqdata$TN + object$freqdata$FP)
```

## Independent validation

Using the same fitted `reitsma` objects:

| Dataset | mada 0.5.12 | Psi only | Denominators only | Full reconstruction |
|---|---:|---:|---:|---:|
| AuditC6 | 52.5234% | 82.2689% | 49.5947% | 80.4936% |
| AuditC | 40.3835% | 86.0467% | 36.7300% | 84.0890% |

The dominant numerical difference in these examples arises from substituting `Var(mu_hat)` for the between-study variance components. The denominator correction is smaller and can either increase or decrease the result.

The reconstructed calculation was also checked against the algebra used by `metadta`, which uses the between-study covariance matrix and the reference-status group sizes.

## Why the substitution depends on the number of studies

Under regular conditions, the variance of a pooled effect estimator decreases approximately as

```text
Var(mu_hat) = O(1/K).
```

If two such variances are inserted into the determinant term, then

```text
|T| = O(1/K^2)
```

and its square root decreases approximately as `O(1/K)`, even when the underlying between-study covariance matrix is held fixed.

A Monte Carlo experiment with fixed population heterogeneity gave:

| K | mada mean I² | reconstructed mean I² | population I² |
|---:|---:|---:|---:|
| 6 | 52.25% | 72.24% | 90.56% |
| 10 | 50.77% | 83.30% | 90.56% |
| 20 | 39.03% | 88.11% | 90.56% |
| 40 | 25.97% | 89.76% | 90.56% |
| 80 | 15.32% | 90.14% | 90.56% |
| 160 | 8.42% | 90.30% | 90.56% |

The simulation is not needed to justify the code change, but it demonstrates the practical consequence of using pooled-estimator variances in place of the between-study components.

## Published examples

Five published estimates that could be reconstructed from aggregate 2×2 data were reproduced before applying the corrected calculation:

| Analysis | mada reproduced | reconstructed | change |
|---|---:|---:|---:|
| SII for thyroid cancer | 71.77% | 90.46% | +18.69 pp |
| AI-PPG for atrial fibrillation | 12.53% | 58.09% | +45.56 pp |
| Single-lead ECG for atrial fibrillation | 9.23% | 53.19% | +43.96 pp |
| Transcranial ultrasound, overall | 14.41% | 50.57% | +36.16 pp |
| Ultrasound for haemorrhagic stroke | 13.12% | 26.49% | +13.37 pp |

These examples demonstrate consequences of the post-estimation calculation only. They do not imply that the underlying Reitsma fits or pooled sensitivity/specificity estimates are incorrect.

## Numerical boundary

When the fitted random-effects correlation is extremely close to ±1, `det(Psi)` may be mathematically zero but become a tiny negative number through floating-point arithmetic. This occurred in a validated TCCS example with an estimated correlation of approximately -0.9999999993.

The proposed safeguard changes only tiny negative determinants within numerical tolerance to zero. A genuinely negative determinant remains an error.

## Scope of the proposed patch

The patch is deliberately narrow. It changes only the Zhou–Dendukuri I² calculation.

Two related scope questions are documented but **not changed automatically**:

1. `summary.reitsma()` can also be called for meta-regression models. The current Zhou calculation uses the first two fixed-effect rows, whereas the published statistic is most straightforwardly defined for the intercept-only bivariate model.
2. `reitsma()` permits non-logit transformations through `alphasens` and `alphafpr`, whereas the exponential expression used in the Zhou calculation corresponds to the logit formulation.

These points may merit an explicit restriction, warning, or separate derivation, but they should be decided by the maintainers rather than folded into the core correction without further validation.

## Reproducibility

Full code, aggregate 2×2 data, software versions, simulation seed, and validated outputs are archived at:

https://doi.org/10.5281/zenodo.22887677
