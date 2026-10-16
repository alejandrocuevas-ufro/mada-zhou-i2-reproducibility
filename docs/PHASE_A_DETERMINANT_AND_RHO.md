# Phase A — Determinant, singularity and extreme correlations

**Status:** reproducible mathematical characterisation; not a Monte Carlo study or
a claim about inferential performance. This is a separate investigation from
the implementation discrepancies in `mada` 0.5.12.

**Scope:** the two-dimensional Zhou–Dendukuri determinant-based statistic,
with the random-effects covariance `T` (positive semidefinite) and the
expected within-study variance matrix `W` (positive definite). In the
original bivariate DTA parametrisation of the existing scripts, `W` is
diagonal, and the coordinates are logit(sensitivity) and logit(FPR).

## 1. Definition and zero-set theorem

```text
I_ZD(T,W) = sqrt(det(T)) / [sqrt(det(T)) + sqrt(det(W))]
```

**Proposition A1.** For any symmetric positive-semidefinite 2x2 `T`
and symmetric positive-definite 2x2 `W`:

```text
I_ZD(T,W) = 0   if and only if   det(T) = 0
           if and only if        rank(T) < 2.
```

**Proof.** Both square roots are nonnegative; `det(W)>0`, so the
denominator is strictly positive. The ratio is therefore zero precisely
when `det(T)=0`. For a PSD 2x2 matrix the determinant is the product
of its nonnegative eigenvalues, so it vanishes exactly when the matrix
has rank 0 or 1. QED.

**Interpretation:** a rank-one between-study covariance may encode
substantial variation in one direction. A zero index means that the
*two-dimensional generalized variance* vanishes; it does not prove that
all marginal between-study variances vanish. When the scientific target
is inherently two-dimensional volume, that behaviour is mathematically
coherent rather than necessarily a defect. It becomes potentially
misleading if reported as "no heterogeneity between studies".

## 2. Correlation limit and rate (with fixed variances)

With marginal between-study variances `v1>0, v2>0` and `|rho|<=1`:

```text
T = [ v1                   rho*sqrt(v1*v2) ]
    [ rho*sqrt(v1*v2)      v2              ]
det(T) = v1*v2*(1-rho^2)
```

Let `delta = 1-|rho|`, `c = sqrt(v1*v2)/sqrt(det(W))`,
and `h = sqrt(2*delta-delta^2)`. Then:

```text
I_ZD = c*h / (1+c*h)
I_ZD ~ c*sqrt(2*delta)               as delta -> 0+
```

This is a square-root approach to zero. For `r=|rho|` in `(0,1)`:

```text
dI_ZD/dr = -c*r / [sqrt(1-r^2)*(1+c*sqrt(1-r^2))^2]
```

Hence the index is strictly decreasing in `|rho|` for fixed positive
marginal variances and `W`, and the magnitude of this slope diverges
as `|rho| -> 1`. The index is even in `rho`: positive and negative
correlations of equal magnitude produce the same value for fixed
marginal parameters.

In the DTA model, replacing FPR with specificity changes the *sign*
of the estimated between-study correlation, not its absolute value or
this determinant.

## 3. Rank-one counterexample

Set:

```text
T(c) = c * [1 -1; -1 1],   c > 0
W = diag(0.2,0.2).
```

Then `det(T(c))=0`, `I_ZD=0` for **all** `c>0`, but
`tr(T(c))=2c` and the largest eigenvalue of `T(c)` is `2c`,
both unbounded as `c->infinity`. These are deterministic algebraic
counterexamples, not fitted DTA models.

## 4. Spectral decomposition and interpretation

Define the *whitened random-effect matrix*
`G = W^(-1/2) T W^(-1/2)`. Let its nonnegative eigenvalues be
`lambda1 <= lambda2`. Then

```text
det(T)/det(W) = lambda1*lambda2
I_ZD = sqrt(lambda1*lambda2)/(1+sqrt(lambda1*lambda2))
```

Thus one vanishing direction forces `I_ZD` to zero even with a
nonzero or very large eigenvalue in the other direction. For
diagnostic descriptions, compare each direction's
`lambda_j/(1+lambda_j)`, which is **not** the published
Zhou–Dendukuri index.

Exploratory comparator (NOT the JWR I_R^2 estimator):

```text
V(T,W) = 1 - sqrt(det(W)/det(W+T))
       = 1 - 1/sqrt((1+lambda1)*(1+lambda2)).
```

This is a determinant-of-total *descriptive comparator* whose values
can stay positive for rank-one `T`. It is not a validated
replacement for the published index and should not be called the
"proportion of heterogeneity" without additional theory.

The unstandardized trace ratio
`tr(T)/(tr(T)+tr(W))` is shown only as an illustrative,
**coordinate-dependent** comparator; it is not an affine-invariant
multivariate `I^2`.

## 5. Invariance theorem and scope of the DTA estimator

For any invertible real 2x2 matrix `A`, let
`T* = A T A^T`, `W* = A W A^T`. Because

```text
det(A T A^T) = det(A)^2 det(T)
det(A W A^T) = det(A)^2 det(W),
```

the *matrix determinant ratio* in `I_ZD` is unchanged.
The generalized eigenvalues of the pair `(T,W)` are also unchanged,
as are the two determinant-based descriptive functions above.

**Important qualifier.** Under a general mixing of outcomes, the
transformed `W*` may be non-diagonal. The original
Zhou–Dendukuri DTA plug-in uses diagonal expected within-study
variances in its chosen coordinates. Consequently, the theorem applies
to the algebraic matrix functional with *both matrices transformed*;
it does not establish reparametrisation invariance of a workflow that
recomputes only diagonal variances after mixing outcomes.

## 6. Numerical precision and statistical identifiability

For a correlation matrix with `|rho|<1` the 2-norm condition number is

```text
kappa(R) = (1+|rho|)/(1-|rho|).
```

As `|rho|->1`, the matrix becomes ill-conditioned and determinant
calculations can be sensitive to floating-point errors. Use
`(1-|rho|)*(1+|rho|)`, a PSD/eigenvalue check and a recorded numeric
tolerance for diagnostics; never silently interpret a negative
determinant of a truly indefinite covariance estimate as valid.
Clipping a tiny negative eigenvalue within tolerance is a *numerical
policy*, not a statistical solution.

If an optimiser estimates `|rho|=1`, this can indicate a true
population boundary or limited information/unstable estimation.
The algebraic result alone cannot distinguish these explanations.
Estimation behaviour requires Phase B and simulation requires Phase C.

## 7. Archived empirical case (descriptive plug-in illustration)

The Pinto-Villalba TCCS ischaemic reconstruction already archived in
`results/published_reanalyses_all_variants.csv` gives approximately:

| Quantity | Value |
|---|---:|
| n analysis arms | 8 |
| Between-study var(logit sensitivity) | 0.102655735 |
| Between-study var(logit FPR) | 2.442952258 |
| Correlation (sensitivity / FPR) | -0.999999999317 |
| Reconstructed Zhou–Dendukuri I2 | 0.00256359% |
| Approximate `kappa` of the correlation matrix | 2.93 billion |
| Largest generalized-direction diagnostic | 63.86% |
| Total-determinant exploratory comparator | 39.88% |

**Caveats:** these are fitted *plug-in* parameters, not known population
heterogeneity. The empirical near-boundary value is not evidence that
the true correlation equals -1. The quoted alternative percentages
measure different descriptive constructs, not competing estimates of
a single unique quantity.

## 8. How to reproduce

From the repository root (base R only):

```bash
Rscript scripts/06_test_phase_a.R
```

The test sources `scripts/06_phase_a_mathematical_properties.R`.
Saved outputs are under `phase_A_outputs/`: numeric grid, rank-one
scaling, transformation experiment, archived case, high-resolution
PNG/PDF plot, and `sessionInfo.txt`. GitHub Actions' lightweight
workflow also runs this test and uploads the outputs.

## 9. Claims supported versus open questions

**Mathematically established:** the zero-set theorem, square-root
boundary limit, monotone dependence on `|rho|` at fixed variances,
evenness in `rho`, rank-one counterexample, spectral explanation and
matrix-congruence invariance.

**Not yet established:** how frequently fitted Reitsma models have
boundary correlations; whether model instability or genuine rank-one
variation causes an individual empirical near-zero estimate; the
calibration, power or inferential utility of candidate complementary
indices; and whether a particular alternative improves decision-making.

## References

- Zhou Y, Dendukuri N. Statistics for quantifying heterogeneity in
  univariate and bivariate meta-analyses of binary data: The case of
  meta-analyses of diagnostic accuracy. *Statistics in Medicine*.
  2014;33:2701–2717. https://doi.org/10.1002/sim.6115
- Jackson D, White IR, Riley RD. Quantifying the impact of between-study
  heterogeneity in multivariate meta-analyses. *Statistics in Medicine*.
  2012;31:3805–3820. https://doi.org/10.1002/sim.5453
