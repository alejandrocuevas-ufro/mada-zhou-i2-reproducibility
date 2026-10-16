# Phase D1, part 4: oracle arm with known within-study variances.
#
# THE QUESTION
# Boundary solutions could come from (a) the bivariate covariance estimator
# itself, which has a stationary point at |rho| = 1 for every dataset, or
# (b) the binomial-to-normal approximation and its continuity correction.
# Phase C cannot separate them because both act at once.
#
# THE CONTRAST
# The oracle arm removes (b) entirely: study-level effects are drawn from the
# same bivariate normal using the SAME latent z1,z2 as the mada arm, but the
# observed outcomes are drawn directly from the normal-normal model with the
# within-study covariance known and supplied to the estimator.  No counts, no
# logit transform, no continuity correction.  The estimator is the identical
# REML objective validated in scripts/12_test_phase_d1.R.
#
# If boundary attraction persists here, it belongs to the covariance
# estimator.  If it largely disappears, it belongs to the approximation.
#
# Run from repository root:  Rscript scripts/14_phase_d1_oracle.R
source("scripts/10_rng_streams.R")
source("scripts/12_phase_d1_audit_fns.R")

N_REP <- as.integer(Sys.getenv("D1_ORACLE_REPS", "400"))
detect_cores <- function(default = 4L) {
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(n) || is.na(n) || n < 1L) n <- suppressWarnings(parallel::detectCores())
  if (!is.finite(n) || is.na(n) || n < 1L) n <- default
  max(1L, min(8L, as.integer(n)))
}
N_CORES <- detect_cores()
outdir <- "results"; dir.create(outdir, showWarnings = FALSE)
mu <- c(2.0997, -1.2637); SD <- c(1.1746, 0.6381)

# Multi-start REML fit of the bivariate normal-normal model with known S.
fit_known_S <- function(Y, S) {
  starts <- list(c(1, 0, 0.5), c(0.5, 0.2, 0.2), c(1.5, -0.3, 0.8),
                 c(0.8, 0.5, 0.05), c(1.2, -0.6, 0.3))
  best <- NULL
  for (st in starts) {
    o <- tryCatch(stats::optim(st, neg2ll_chol, Y = Y, S = S, reml = TRUE,
                               method = "Nelder-Mead",
                               control = list(maxit = 3000, reltol = 1e-12)),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && (is.null(best) || o$value < best$value))
      best <- o
  }
  best
}

oracle_one <- function(scenario, rep) {
  K <- scenario$K; rho <- scenario$rho_true
  pool <- draw_latent_pool(1L, rep, 160L)          # same latent studies as mada arm
  idx <- seq_len(K)
  theta1 <- mu[1] + SD[1] * pool$z1[idx]
  theta2 <- mu[2] + SD[2] * (rho * pool$z1[idx] +
                               sqrt((1 - abs(rho)) * (1 + abs(rho))) * pool$z2[idx])
  # Known within-study variances.  They MUST be independent of theta: if s_i
  # were computed from the study's own latent p_i (as a first version of this
  # script did), conditioning on s_i would be informative about theta_i, the
  # normal-normal model would be misspecified, and tau2 would be attenuated
  # even with perfect data - at K = 4000 that error produced a ratio of 0.80
  # instead of 1.  Drawing the variance-generating proportions from an
  # independent copy of the same marginal distribution keeps the realistic
  # spread of within-study precision while restoring correct specification.
  set.seed(rng_seed_for(40L, rep) + K, kind = RNG_KIND)
  q1 <- plogis(mu[1] + SD[1] * rnorm(K))
  q2 <- plogis(mu[2] + SD[2] * rnorm(K))
  s1 <- 1 / (scenario$nD * q1 * (1 - q1))
  s2 <- 1 / (scenario$nND * q2 * (1 - q2))
  # Measurement noise on its own stream so the latent studies stay shared.
  set.seed(rng_seed_for(30L, rep) + K, kind = RNG_KIND)
  Y <- cbind(theta1 + rnorm(K, sd = sqrt(s1)), theta2 + rnorm(K, sd = sqrt(s2)))
  S <- cbind(s1, s2)

  out <- data.frame(scenario_label = scenario$label, arm = "oracle_normal_known_S",
                    rho_true = rho, K = K, replicate = rep, status = "ok",
                    rho_hat = NA_real_, exact_boundary = NA, tau2_se = NA_real_,
                    tau2_fpr = NA_real_, I2 = NA_real_, stringsAsFactors = FALSE)
  o <- fit_known_S(Y, S)
  if (is.null(o)) { out$status <- "fit_error"; return(out) }
  L11 <- o$par[1]; L21 <- o$par[2]; L22 <- o$par[3]
  Psi <- matrix(c(L11^2, L11 * L21, L11 * L21, L21^2 + L22^2), 2, 2)
  rho_hat <- Psi[1, 2] / sqrt(Psi[1, 1] * Psi[2, 2])
  rho_hat <- max(-1, min(1, rho_hat))
  out$rho_hat <- rho_hat
  out$exact_boundary <- abs(rho_hat) >= 1 - 1e-12
  out$tau2_se <- Psi[1, 1]; out$tau2_fpr <- Psi[2, 2]
  # Zhou index with the SAME known within-study variances (mean of 1/precision
  # is replaced by the mean of the known variances, which is the matching
  # quantity when S is known and diagonal).
  detT <- max(0, Psi[1, 1] * Psi[2, 2] - Psi[1, 2]^2)
  detW <- mean(s1) * mean(s2)
  out$I2 <- sqrt(detT) / (sqrt(detT) + sqrt(detW))
  out
}

scenarios <- data.frame(
  label = c("rho0_K6", "rho0_K20", "rho0.9_K160", "rho0.99_K160"),
  rho_true = c(0, 0, 0.9, 0.99), K = c(6L, 20L, 160L, 160L),
  nD = 50L, nND = 50L, stringsAsFactors = FALSE)

jobs <- expand.grid(s = seq_len(nrow(scenarios)), rep = seq_len(N_REP))
cat("Oracle arm:", nrow(jobs), "fits on", N_CORES, "cores\n")
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(jobs)), function(i)
  oracle_one(scenarios[jobs$s[i], ], jobs$rep[i]), mc.cores = N_CORES)
res <- do.call(rbind, res)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
write.csv(res, file.path(outdir, "D1_oracle_arm.csv"), row.names = FALSE)

ok <- res[res$status == "ok", ]
mcse_p <- function(p, n) sqrt(p * (1 - p) / n)
summ <- do.call(rbind, lapply(scenarios$label, function(lb) {
  x <- ok[ok$scenario_label == lb, ]
  data.frame(scenario = lb, n = nrow(x),
             frac_exact_boundary = mean(x$exact_boundary),
             mcse_boundary = mcse_p(mean(x$exact_boundary), nrow(x)),
             frac_I2_below_5pct = mean(x$I2 < 0.05),
             median_I2 = median(x$I2),
             ratio_tau2_se = median(x$tau2_se) / SD[1]^2,
             ratio_tau2_fpr = median(x$tau2_fpr) / SD[2]^2,
             median_rho_hat = median(x$rho_hat), row.names = NULL)
}))
write.csv(summ, file.path(outdir, "D1_oracle_summary.csv"), row.names = FALSE)
cat("\nD1 ORACLE ARM (normal-normal data, within-study variances known)\n")
print(summ, row.names = FALSE, digits = 4)
cat("\nCompare frac_exact_boundary against the mada arm in",
    "results/D1_convergence_summary.csv\n")
