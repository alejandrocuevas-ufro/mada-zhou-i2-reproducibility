# How to transmit the proposed correction

## Recommended sequence

### 1. Wait for confirmation of the interpretation

Do not push changes to an upstream repository before the maintainer confirms that a contribution would be useful. The CRAN GitHub mirror is not the development upstream and should not receive a pull request.

### 2. Reply in the existing email thread

Keep the discussion in the same thread. Send a concise message with:

- the link to this draft patch directory;
- the permanent Zenodo DOI;
- a one-sentence statement that the Reitsma fit is unchanged;
- an offer to adapt the patch to the maintainer's preferred workflow.

Avoid attaching the full manuscript unless requested.

### 3. Supply the contribution in the format requested

Typical maintainer workflows are:

- **R-Forge/SVN**: send a unified `.patch` file or commit against a checked-out working copy.
- **Git repository**: fork/branch and open a pull request.
- **Email-only maintenance**: send the unified patch and test file as attachments.
- **Maintainer-applied change**: provide the exact replacement lines plus expected regression-test outputs.

The files in this directory support all four workflows.

### 4. Keep the software correction and the academic evidence distinct

For code review, lead with:

- exact source lines;
- proposed replacement;
- AuditC6/AuditC regression tests;
- boundary safeguard.

For scientific documentation, provide separately:

- the five published reanalyses;
- the K-dependence simulation;
- the permanent reproducibility archive.

This makes the patch easy to review without hiding its practical consequences.

### 5. Ask how they want attribution handled

If the correction is accepted, appropriate outcomes may include:

- acknowledgement in `NEWS` or release notes;
- a contributor entry if the maintainer considers the contribution substantial;
- citation/link to the reproducibility archive;
- acknowledgement of the report in the package issue tracker.

Do not prescribe an authorship or contributor role; let the maintainers apply their project's conventions.

## Before any upstream submission

Run against the maintainer's current development version, not only CRAN 0.5.12:

```bash
R CMD build mada
R CMD check --as-cran mada_*.tar.gz
```

and execute the regression script:

```bash
Rscript tests/zhou_i2-regression.R
```

The draft patch should be rebased if the development source has changed since 0.5.12.
