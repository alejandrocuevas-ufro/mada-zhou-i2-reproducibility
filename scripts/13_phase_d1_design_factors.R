# Phase D1, part 3: design factors Phase C held fixed.
#
# Two choices were constants in Phase C and are confounded with its findings:
#
#  (1) CONTINUITY CORRECTION.  reitsma(correction.control = "all") adds 0.5 to
#      EVERY table as soon as a single zero cell appears anywhere in the
#      dataset.  The probability of at least one zero cell grows with K, so in
#      the 50/50 scheme the whole dataset is perturbed in 38% of replicates at
#      K = 6 but in 100% at K = 160.  That is why the attenuation of the
#      between-study variances does not improve with K.  "single" corrects
#      only the offending tables and is the natural comparator.
#
#  (2) EQUAL STUDY SIZES.  Phase C gave every study the same n.  The Zhou
#      index builds its within-study term from mean(1/n), which in a real
#      review is dominated by the smallest study.  Log-normal sizes around the
#      same median test whether the findings survive realistic size spread.
#
# Common random numbers are used throughout, so each factor is compared on the
# same latent studies and the same uniforms.
#
# Run from repository root:  Rscript scripts/13_phase_d1_design_factors.R
suppressPackageStartupMessages(library(mada))
stopifnot(as.character(packageVersion("mada")) == "0.5.12")
source("scripts/zhou_i2_utils.R")
source("scripts/10_rng_streams.R")
source("scripts/12_phase_d1_audit_fns.R")

N_REP <- as.integer(Sys.getenv("D1_FACTOR_REPS", "400"))
detect_cores <- function(default = 4L) {
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(n) || is.na(n) || n < 1L) n <- suppressWarnings(parallel::detectCores())
  if (!is.finite(n) || is.na(n) || n < 1L) n <- default
  max(1L, min(8L, as.integer(n)))
}
N_CORES <- detect_cores()
outdir <- "results"; dir.create(outdir, showWarnings = FALSE)
mu <- c(2.0997, -1.2637); SD <- c(1.1746, 0.6381)

scenarios <- data.frame(
  label = c("rho0_K6", "rho0_K20", "rho0.9_K160", "rho0.99_K160"),
  rho_true = c(0, 0, 0.9, 0.99), K = c(6L, 20L, 160L, 160L),
  stringsAsFactors = FALSE)
factors <- expand.grid(correction_control = c("all", "single"),
                       size_pattern = c("equal", "lognormal"),
                       stringsAsFactors = FALSE)

run_cell <- function(scenario, fac, rep) {
  pool <- draw_latent_pool(1L, rep, 160L)
  if (fac$size_pattern == "equal") {
    nD <- 50L; nND <- 50L
  } else {
    sz <- draw_study_sizes(1L, rep, 160L, median_nD = 50, median_nND = 50)
    nD <- sz$nD; nND <- sz$nND
  }
  counts <- counts_from_pool(pool, scenario$K, scenario$rho_true, nD, nND, mu, SD)
  out <- data.frame(scenario_label = scenario$label, rho_true = scenario$rho_true,
                    K = scenario$K, correction_control = fac$correction_control,
                    size_pattern = fac$size_pattern, replicate = rep,
                    seed = pool$seed,
                    n_zero_arms = sum(rowSums(counts == 0) > 0),
                    status = "ok", rho_hat = NA_real_, exact_boundary = NA,
                    tau2_se = NA_real_, tau2_fpr = NA_real_, I2 = NA_real_,
                    stringsAsFactors = FALSE)
  fit <- tryCatch(suppressWarnings(
    mada::reitsma(counts, method = "reml", correction = 0.5,
                  correction.control = fac$correction_control)),
    error = function(e) e)
  if (inherits(fit, "error")) { out$status <- "fit_error"; return(out) }
  cc <- tryCatch(zhou_components(fit, use_Psi = TRUE, correct_denominators = TRUE),
                 error = function(e) NULL)
  if (is.null(cc)) { out$status <- "extraction_error"; return(out) }
  out$rho_hat <- safe_rho_from_Psi(fit$Psi)
  out$exact_boundary <- abs(out$rho_hat) >= 1 - 1e-12
  out$tau2_se <- fit$Psi[1, 1]; out$tau2_fpr <- fit$Psi[2, 2]; out$I2 <- cc$I2
  out
}

jobs <- expand.grid(s = seq_len(nrow(scenarios)), f = seq_len(nrow(factors)),
                    rep = seq_len(N_REP))
cat("Design factor experiment:", nrow(jobs), "fits on", N_CORES, "cores\n")
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(jobs)), function(i)
  run_cell(scenarios[jobs$s[i], ], factors[jobs$f[i], ], jobs$rep[i]),
  mc.cores = N_CORES)
res <- do.call(rbind, res)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
write.csv(res, file.path(outdir, "D1_design_factors.csv"), row.names = FALSE)

tau2_gen <- SD^2
ok <- res[res$status == "ok", ]
mcse_p <- function(p, n) sqrt(p * (1 - p) / n)
summ <- do.call(rbind, lapply(split(ok, list(ok$scenario_label,
                                             ok$correction_control,
                                             ok$size_pattern), drop = TRUE),
  function(x) data.frame(
    scenario = x$scenario_label[1], correction = x$correction_control[1],
    sizes = x$size_pattern[1], n = nrow(x),
    frac_exact_boundary = mean(x$exact_boundary),
    mcse_boundary = mcse_p(mean(x$exact_boundary), nrow(x)),
    frac_I2_below_5pct = mean(x$I2 < 0.05),
    mcse_I2_low = mcse_p(mean(x$I2 < 0.05), nrow(x)),
    median_I2 = median(x$I2),
    ratio_tau2_se = median(x$tau2_se) / tau2_gen[1],
    ratio_tau2_fpr = median(x$tau2_fpr) / tau2_gen[2],
    frac_any_zero = mean(x$n_zero_arms > 0), row.names = NULL)))
summ <- summ[order(summ$scenario, summ$sizes, summ$correction), ]
write.csv(summ, file.path(outdir, "D1_design_factor_summary.csv"), row.names = FALSE)
cat("\nD1 DESIGN FACTOR SUMMARY\n")
print(summ, row.names = FALSE, digits = 4)
