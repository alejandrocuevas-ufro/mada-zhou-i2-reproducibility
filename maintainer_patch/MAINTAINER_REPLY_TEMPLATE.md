Subject: Re: Zhou–Dendukuri bivariate I² in mada

Dear Dr Doebler,

Thank you for confirming this.

I have prepared a small proposed patch against the `summary.reitsma()` implementation in `mada` 0.5.12, together with base-R regression tests and two short validation scripts.

The proposed change is restricted to the post-estimation Zhou–Dendukuri I² calculation. It uses the diagonal elements of `Psi` as the between-study variances, the reference-status denominators `TP+FN` and `TN+FP`, and a numerical safeguard for a determinant that is slightly negative only because of floating-point roundoff at the positive-semidefinite boundary. It does not change the Reitsma fit or the pooled estimates.

The draft contribution is available here:

https://github.com/alejandrocuevas-ufro/mada-zhou-i2-reproducibility/tree/maintainer-patch-draft/maintainer_patch

The complete reproducibility archive, including the simulation and published reanalyses, is available at:

https://doi.org/10.5281/zenodo.22887677

I would be happy to adapt the patch to your preferred workflow, including R-Forge, a source tarball patch, or another repository if that is more convenient.

Best wishes,

Alejandro Cuevas Villegas
Universidad de La Frontera
Temuco, Chile
