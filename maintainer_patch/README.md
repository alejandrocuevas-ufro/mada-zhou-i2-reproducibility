# Draft contribution to the mada maintainers

This directory packages the proposed correction to the Zhou–Dendukuri bivariate I² calculation in `mada` 0.5.12.

It is intentionally separated from the main reproducibility workflow so that the maintainers can review the software change independently from the manuscript-oriented analyses.

## Contents

- `mada_0.5.12_zhou_i2.patch` — minimal proposed source-code patch for `R/reitsma.R`.
- `tests/zhou_i2-regression.R` — base-R regression tests suitable for the package's `tests/` directory.
- `run_public_validation.R` — self-contained validation using AuditC6 and the packaged AuditC dataset.
- `run_boundary_example.R` — numerical-boundary example with an estimated random-effects correlation close to -1.
- `TECHNICAL_NOTE.md` — rationale, scope, validated numerical consequences, and points requiring maintainer judgement.
- `EXPECTED_RESULTS.md` — validated outputs from the reproducibility archive.
- `MAINTAINER_REPLY_TEMPLATE.md` — short reply to use if the maintainers invite a patch or pull request.
- `TRANSMISSION_GUIDE.md` — suggested way to transmit the contribution.

## Proposed core change

The proposed patch changes only the post-estimation Zhou–Dendukuri I² calculation in `summary.reitsma()`.

It does **not** alter:

- the Reitsma model fit;
- pooled sensitivity or specificity;
- the estimated between-study covariance matrix `Psi`;
- the Holling I² calculations.

The proposed changes are:

1. use the diagonal elements of `Psi` as the between-study variance components;
2. use `TP + FN` and `TN + FP` as the reference-status group sizes;
3. compute the determinant directly from `Psi`;
4. protect only tiny negative determinants attributable to floating-point error at the positive-semidefinite boundary.

## Reproducibility archive

Full analyses, aggregate data, simulation code, and validated outputs:

https://doi.org/10.5281/zenodo.22887677

Development repository:

https://github.com/alejandrocuevas-ufro/mada-zhou-i2-reproducibility
