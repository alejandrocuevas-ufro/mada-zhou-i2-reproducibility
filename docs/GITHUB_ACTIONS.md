# Running this reproducibility project on GitHub Actions

The workflows in `.github/workflows/` execute R on GitHub-hosted Ubuntu runners; no local R installation is required.

## Automated lightweight validation

**R | determinant boundary checks** (`r-boundary-check.yml`) runs on pull requests and pushes to `main` and can also be started using **Actions > R | determinant boundary checks > Run workflow**. It installs R with `r-lib/actions/setup-r@v2`, requires no additional R packages, runs `scripts/05_test_determinant_boundary.R`, and saves its CSV and PNG outputs as the `determinant-boundary-results` artifact. Tests cover the limit at `rho = +/-1`, the `rho`-sign symmetry, the monotone determinant collapse for fixed marginal variances, and consistency with the archived Pinto-Villalba TCCS ischaemic example.

## Full archived audit and Monte Carlo simulation

**R | full mada 0.5.12 reproduction** (`r-full-reproduction.yml`) is **manual only** to avoid accidentally executing 3,000 REML fits after every commit.

1. Open the repository **Actions** tab.
2. Select **R | full mada 0.5.12 reproduction**.
3. Click **Run workflow**, choose `main`, and confirm.
4. When finished, open the run summary and download the `full-mada-reproduction` artifact, which contains `results/` and the determinant diagnostics.

The full workflow installs precisely `mada` 0.5.12 from CRAN archives using `remotes::install_version()`, runs the existing `Rscript run_all.R`, and executes the determinant checks. Its outputs include fitted-model audit tables, the original Monte Carlo figures and replicate data, and `sessionInfo()`. It has a 300-minute job limit. Workflow runs use read-only repository permissions and do not commit generated results to `main`.

**Notes:** GitHub offers manual `workflow_dispatch` runs from the default branch after a workflow has been merged there. The first automated boundary check is expected on this pull request if repository Actions are enabled. A passing lightweight check does *not* imply that the full archived `mada` run has been verified; the latter must be run separately. Uploaded artifacts are retained for 30 days, subject to repository settings.

To reproduce the fast check locally: `Rscript scripts/05_test_determinant_boundary.R`.
To reproduce the full run locally after installing dependencies: `Rscript run_all.R`.
