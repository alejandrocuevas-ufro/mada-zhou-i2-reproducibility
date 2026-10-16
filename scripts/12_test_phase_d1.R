# Validation for the Phase D1 audit machinery.
# The Hessian and profile diagnostics are only meaningful if neg2ll_chol() is
# the same objective mada maximises.  These checks establish that against
# mada's own fitted output, on datasets with and without zero cells.
# Run from repository root:  Rscript scripts/12_test_phase_d1.R
suppressPackageStartupMessages({library(mada); library(numDeriv)})
source("scripts/zhou_i2_utils.R")
source("scripts/10_rng_streams.R")
source("scripts/12_phase_d1_audit_fns.R")

mu <- c(2.0997, -1.2637); SD <- c(1.1746, 0.6381)
checks <- list()
add <- function(name, pass, detail = "") checks[[length(checks) + 1L]] <<-
  data.frame(check = name, pass = pass, detail = detail, stringsAsFactors = FALSE)

cases <- list(list(K = 20L, rho = 0.0, rep = 3L),
              list(K = 6L,  rho = 0.0, rep = 11L),
              list(K = 60L, rho = 0.9, rep = 5L))

for (ci in seq_along(cases)) {
  cs <- cases[[ci]]
  pool <- draw_latent_pool(1L, cs$rep, 160L)
  counts <- counts_from_pool(pool, cs$K, cs$rho, 50L, 50L, mu, SD)
  fit <- suppressWarnings(mada::reitsma(counts, method = "reml",
                                        correction = .5, correction.control = "all"))
  tr <- mada_transform(counts, 0.5, "all")

  # (1) Gradient of our objective at mada's optimum must be ~0.
  gr <- numDeriv::grad(neg2ll_chol, fit$par, Y = tr$Y, S = tr$S, reml = TRUE)
  scale <- max(1, abs(neg2ll_chol(fit$par, tr$Y, tr$S, TRUE)))
  add(sprintf("case%d_gradient_zero_at_mada_optimum", ci),
      max(abs(gr)) / scale < 1e-4,
      sprintf("max|grad|/scale = %.3e", max(abs(gr)) / scale))

  # (2) Free optimisation of our objective must recover mada's Psi.
  o <- stats::optim(fit$par, neg2ll_chol, Y = tr$Y, S = tr$S, reml = TRUE,
                    method = "Nelder-Mead",
                    control = list(maxit = 4000, reltol = 1e-12))
  L <- matrix(c(o$par[1], o$par[2], 0, o$par[3]), 2, 2)
  Psi_ours <- L %*% t(L)
  add(sprintf("case%d_recovers_mada_Psi", ci),
      max(abs(Psi_ours - fit$Psi)) < 1e-4 * max(1, max(abs(fit$Psi))),
      sprintf("max abs diff = %.3e", max(abs(Psi_ours - fit$Psi))))

  # (3) Our objective must not find a materially better optimum than mada's.
  add(sprintf("case%d_no_better_optimum_found", ci),
      o$value >= neg2ll_chol(fit$par, tr$Y, tr$S, TRUE) - 1e-6,
      sprintf("ours %.6f vs mada-par %.6f", o$value,
              neg2ll_chol(fit$par, tr$Y, tr$S, TRUE)))

  # (4) The GLS mean implied by our transform must match mada's coefficients.
  v1 <- fit$par[1]^2; c12 <- fit$par[1] * fit$par[2]
  v2 <- fit$par[2]^2 + fit$par[3]^2
  a <- tr$S[, 1] + v1; b <- tr$S[, 2] + v2; d <- a * b - c12^2
  i11 <- b/d; i12 <- -c12/d; i22 <- a/d
  A <- matrix(c(sum(i11), sum(i12), sum(i12), sum(i22)), 2, 2)
  beta <- solve(A, c(sum(i11*tr$Y[,1] + i12*tr$Y[,2]),
                     sum(i12*tr$Y[,1] + i22*tr$Y[,2])))
  add(sprintf("case%d_coefficients_match", ci),
      max(abs(beta - as.numeric(fit$coefficients))) < 1e-6,
      sprintf("max abs diff = %.3e", max(abs(beta - as.numeric(fit$coefficients)))))
}

# (5) The "single" branch of the transform must differ from "all" exactly when
# some but not all tables contain a zero.
pool <- draw_latent_pool(1L, 2L, 160L)
cz <- counts_from_pool(pool, 40L, 0.0, 50L, 50L, mu, SD)
if (any(as.matrix(cz[, c("TP","FN","FP","TN")]) == 0)) {
  ta <- mada_transform(cz, 0.5, "all"); ts <- mada_transform(cz, 0.5, "single")
  add("transform_all_differs_from_single", !isTRUE(all.equal(ta$Y, ts$Y)),
      "dataset contains at least one zero cell")
} else {
  add("transform_all_differs_from_single", TRUE, "no zero cells in probe dataset")
}

res <- do.call(rbind, checks)
dir.create("results", showWarnings = FALSE)
write.csv(res, "results/D1_likelihood_validation.csv", row.names = FALSE)
print(res, row.names = FALSE)
if (!all(res$pass)) stop("Phase D1 likelihood validation failed.")
cat("\nAll Phase D1 likelihood validation checks passed.\n")
