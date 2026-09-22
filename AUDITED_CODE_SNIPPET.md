# Audited `mada` code fragment

Audit target: `mada` 0.5.12, function `summary.reitsma()`.

The implementation assigns the squared standard errors of the pooled fixed-effect coefficients to the quantities used as between-study variances:

```r
tau_a <- (tabfixed[1,2])^2
tau_b <- (tabfixed[2,2])^2
cor_sq <- (corRandom[1,2])^2
tau <- (1-cor_sq)*tau_a*tau_b
```

In the fitted `reitsma` object, the between-study covariance matrix is `Psi`. The reconstruction of the Zhou-Dendukuri determinant therefore uses:

```r
tau_a <- Psi[1,1]
tau_b <- Psi[2,2]
cor_sq <- (corRandom[1,2])^2
tau <- (1-cor_sq)*tau_a*tau_b
```

The second discrepancy concerns the denominators used for the expected within-study variances. The audited implementation uses:

```r
Nia <- 1/(object$freqdata$TP + object$freqdata$FP)
Nib <- 1/(object$freqdata$TN + object$freqdata$FN)
```

whereas sensitivity and specificity are based on the reference-status groups:

```r
Nia <- 1/(object$freqdata$TP + object$freqdata$FN)
Nib <- 1/(object$freqdata$TN + object$freqdata$FP)
```

The repository scripts reproduce the native `mada` result before applying either substitution. Numerical safeguards in `scripts/zhou_i2_utils.R` only clip correlations that exceed `[-1,1]` within floating-point tolerance at the positive-semidefinite boundary.
