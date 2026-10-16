# Phase D1, part 2: convergence and identifiability audit.
#
# Phase C reported "15,000 fits, 0 failures".  That statement only means
# reitsma() returned an object.  This script asks whether those objects are
# well-identified, by recording for every replicate:
#
#   * mada/mvmeta's own convergence flag and iteration count;
#   * the fitted Cholesky factor.  mvmeta parameterises Psi = L L' with
#     L = [[L11,0],[L21,L22]] and optimises L UNCONSTRAINED, so |rho| = 1 is
#     the interior event L22 = 0 rather than a box constraint.  |L22| is the
#     natural distance to singularity and is recorded directly;
#   * eigenvalues and condition number of the estimated Psi;
#   * eigenvalues and condition number of the numerical Hessian of the REML
#     objective at the optimum - a near-zero eigenvalue is weak identification
#     regardless of what the convergence flag says;
#   * the range of the profile objective over rho, i.e. how much the data
#     actually say about the correlation.
#
# Run from repository root:  Rscript scripts/12_phase_d1_convergence_audit.R
suppressPackageStartupMessages({library(mada); library(numDeriv)})
stopifnot(as.character(packageVersion("mada")) == "0.5.12")
source("scripts/zhou_i2_utils.R")
source("scripts/10_rng_streams.R")

N_REP <- as.integer(Sys.getenv("D1_AUDIT_REPS", "400"))
detect_cores <- function(default = 4L) {
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(n) || is.na(n) || n < 1L) n <- suppressWarnings(parallel::detectCores())
  if (!is.finite(n) || is.na(n) || n < 1L) n <- default
  max(1L, min(8L, as.integer(n)))
}
N_CORES <- detect_cores()
outdir <- "results"; dir.create(outdir, showWarnings = FALSE)

mu <- c(2.0997, -1.2637); SD <- c(1.1746, 0.6381)

source("scripts/12_phase_d1_audit_fns.R")

audit_one <- function(scenario, rep) {
  pool <- draw_latent_pool(scenario$scheme_id, rep, 160L)
  counts <- counts_from_pool(pool, scenario$K, scenario$rho_true,
                             scenario$nD, scenario$nND, mu, SD)
  base <- data.frame(scenario_label = scenario$label, rho_true = scenario$rho_true,
                     K = scenario$K, scheme = scenario$scheme, replicate = rep,
                     seed = pool$seed, n_zero_arms = sum(rowSums(counts == 0) > 0),
                     correction_triggered = any(as.matrix(
                       counts[, c("TP","FN","FP","TN")]) == 0),
                     status = "ok", converged = NA, niter = NA_integer_,
                     L11 = NA_real_, L21 = NA_real_, L22 = NA_real_,
                     abs_L22 = NA_real_, rho_hat = NA_real_,
                     exact_boundary = NA, tau2_se = NA_real_, tau2_fpr = NA_real_,
                     psi_eig_min = NA_real_, psi_cond = NA_real_, I2 = NA_real_,
                     hess_eig_min = NA_real_, hess_cond = NA_real_,
                     profile_range = NA_real_,
                     sing_curvature = NA_real_, sing_is_local_min = NA,
                     sing_excess_deviance = NA_real_, stringsAsFactors = FALSE)

  fit <- tryCatch(suppressWarnings(
    mada::reitsma(counts, method = "reml", correction = .5,
                  correction.control = "all")), error = function(e) e)
  if (inherits(fit, "error")) { base$status <- "fit_error"; return(base) }

  base$converged <- isTRUE(fit$converged)
  base$niter <- as.integer(fit$niter)
  base$L11 <- fit$par[1]; base$L21 <- fit$par[2]; base$L22 <- fit$par[3]
  base$abs_L22 <- abs(fit$par[3])
  rho <- safe_rho_from_Psi(fit$Psi)
  base$rho_hat <- rho
  base$exact_boundary <- abs(rho) >= 1 - 1e-12
  base$tau2_se <- fit$Psi[1, 1]; base$tau2_fpr <- fit$Psi[2, 2]
  ev <- eigen(fit$Psi, symmetric = TRUE, only.values = TRUE)$values
  base$psi_eig_min <- min(ev)
  base$psi_cond <- if (min(ev) > 0) max(ev) / min(ev) else Inf
  cc <- tryCatch(zhou_components(fit, use_Psi = TRUE, correct_denominators = TRUE),
                 error = function(e) NULL)
  if (!is.null(cc)) base$I2 <- cc$I2

  tr <- mada_transform(counts, 0.5, "all")
  h <- tryCatch(numDeriv::hessian(neg2ll_chol, fit$par, Y = tr$Y, S = tr$S,
                                  reml = TRUE), error = function(e) NULL)
  if (!is.null(h) && all(is.finite(h))) {
    he <- eigen((h + t(h)) / 2, symmetric = TRUE, only.values = TRUE)$values
    base$hess_eig_min <- min(he)
    base$hess_cond <- if (min(abs(he)) > 0) max(abs(he)) / min(abs(he)) else Inf
  }
  pr <- profile_over_rho(tr$Y, tr$S, reml = TRUE)
  if (any(is.finite(pr))) base$profile_range <- diff(range(pr, na.rm = TRUE))

  sd_ <- singular_point_diagnostics(
    tr$Y, tr$S, reml = TRUE,
    unconstrained_value = neg2ll_chol(fit$par, tr$Y, tr$S, TRUE))
  base$sing_curvature <- sd_$curvature
  base$sing_is_local_min <- sd_$is_local_min
  base$sing_excess_deviance <- sd_$excess
  base
}

scenarios <- data.frame(
  label = c("rho0_K6", "rho0_K20", "rho0.9_K160", "rho0.99_K160", "rho0.99_K6"),
  rho_true = c(0, 0, 0.9, 0.99, 0.99),
  K = c(6L, 20L, 160L, 160L, 6L),
  scheme = "low_precision", scheme_id = 1L, nD = 50L, nND = 50L,
  stringsAsFactors = FALSE)

jobs <- expand.grid(s = seq_len(nrow(scenarios)), rep = seq_len(N_REP))
cat("Convergence audit:", nrow(jobs), "fits on", N_CORES, "cores\n")
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(jobs)), function(i)
  audit_one(scenarios[jobs$s[i], ], jobs$rep[i]), mc.cores = N_CORES)
res <- do.call(rbind, res)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
write.csv(res, file.path(outdir, "D1_convergence_audit.csv"), row.names = FALSE)

ok <- res[res$status == "ok", ]
summ <- do.call(rbind, lapply(scenarios$label, function(lb) {
  x <- ok[ok$scenario_label == lb, ]
  data.frame(scenario = lb, n = nrow(x),
             frac_converged_flag = mean(x$converged),
             frac_exact_boundary = mean(x$exact_boundary),
             frac_converged_AND_boundary = mean(x$converged & x$exact_boundary),
             median_abs_L22 = median(x$abs_L22),
             frac_abs_L22_below_1e3 = mean(x$abs_L22 < 1e-3),
             median_hess_eig_min = median(x$hess_eig_min, na.rm = TRUE),
             frac_hess_nonpositive = mean(x$hess_eig_min <= 0, na.rm = TRUE),
             median_hess_cond = median(x$hess_cond, na.rm = TRUE),
             median_profile_range = median(x$profile_range, na.rm = TRUE),
             frac_correction_triggered = mean(x$correction_triggered),
             frac_I2_below_5pct = mean(x$I2 < 0.05, na.rm = TRUE),
             frac_singular_is_local_min = mean(x$sing_is_local_min, na.rm = TRUE),
             median_singular_excess_deviance = median(x$sing_excess_deviance, na.rm = TRUE),
             row.names = NULL)
}))
write.csv(summ, file.path(outdir, "D1_convergence_summary.csv"), row.names = FALSE)
cat("\nD1 CONVERGENCE AND IDENTIFIABILITY SUMMARY\n")
print(summ, row.names = FALSE, digits = 4)
cat("\nA convergence flag of TRUE alongside a boundary solution and a near-zero",
    "\nHessian eigenvalue is the signature the Phase C status column could not see.\n")
