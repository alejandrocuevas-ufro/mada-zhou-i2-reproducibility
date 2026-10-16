# Phase D1, part 1: pseudo-true values of the fitting procedure.
#
# WHY THIS EXISTS
# Phase C compared the fitted Zhou-Dendukuri I2 against a POPULATION target
# computed from the generating parameters.  That comparison is only a bias if
# the estimator converges to the generating parameters.  It does not: the
# normal-normal approximation combined with a +0.5 continuity correction on
# every table converges to a different point.  This script estimates that
# point - the pseudo-true value - by running the identical fitting procedure
# on datasets large enough (K = 4000) for Monte Carlo error to be negligible.
#
# The pseudo-true value lets Phase D separate
#     approximation bias = pseudo-true - population
#     estimation bias    = fitted      - pseudo-true
# which the single "bias" column of Phase C conflated.
#
# Run from repository root:  Rscript scripts/11_phase_d1_pseudotrue.R
suppressPackageStartupMessages(library(mada))
stopifnot(as.character(packageVersion("mada")) == "0.5.12")
source("scripts/zhou_i2_utils.R")
source("scripts/10_rng_streams.R")

K_LARGE <- as.integer(Sys.getenv("D1_PSEUDOTRUE_K", "4000"))
N_REP <- as.integer(Sys.getenv("D1_PSEUDOTRUE_REPS", "5"))
detect_cores <- function(default = 4L) {
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(n) || is.na(n) || n < 1L)
    n <- suppressWarnings(parallel::detectCores())
  if (!is.finite(n) || is.na(n) || n < 1L) n <- default
  max(1L, min(8L, as.integer(n)))
}
N_CORES <- detect_cores()
outdir <- "results"; dir.create(outdir, showWarnings = FALSE)

mu <- c(2.0997, -1.2637)
SD <- c(1.1746, 0.6381)
tau2_gen <- SD^2
rho_grid <- c(-1, -.9999, -.99, -.9, -.5, 0, .5, .9, .99, .9999, 1)
schemes <- data.frame(scheme = c("low_precision", "original_precision"),
                      scheme_id = c(1L, 2L), nD = c(50L, 150L),
                      nND = c(50L, 500L), stringsAsFactors = FALSE)

# Population target exactly as Phase C defined it, for comparability.
population_I2 <- function(rho, nD, nND) {
  w1 <- (exp(mu[1] + tau2_gen[1]/2) + exp(-mu[1] + tau2_gen[1]/2) + 2) / nD
  th2 <- -mu[2]
  w2 <- (exp(th2 + tau2_gen[2]/2) + exp(-th2 + tau2_gen[2]/2) + 2) / nND
  rootT <- prod(SD) * sqrt((1 - abs(rho)) * (1 + abs(rho)))
  c(pop_I2 = rootT / (sqrt(w1 * w2) + rootT), pop_w_se = w1, pop_w_fpr = w2)
}

design <- do.call(rbind, lapply(seq_len(nrow(schemes)), function(j)
  data.frame(schemes[j, ], rho_true = rho_grid, stringsAsFactors = FALSE, row.names = NULL)))
design$cell_id <- seq_len(nrow(design))

one_fit <- function(cell, rep) {
  # A distinct base seed (not a large replicate index, which would overflow
  # the integer seed map) keeps these large-K draws disjoint from the ordinary
  # replicate stream used by the audit scripts.
  pool <- draw_latent_pool(cell$scheme_id + 10L, rep, K_LARGE,
                           base_seed = RNG_BASE_SEED + 777L)
  d <- counts_from_pool(pool, K_LARGE, cell$rho_true, cell$nD, cell$nND, mu, SD)
  f <- tryCatch(mada::reitsma(d, method = "reml", correction = .5,
                              correction.control = "all"),
                error = function(e) e)
  if (inherits(f, "error"))
    return(data.frame(ok = FALSE, tau2_se = NA_real_, tau2_fpr = NA_real_,
                      rho = NA_real_, I2 = NA_real_, mu_se = NA_real_,
                      mu_fpr = NA_real_))
  cc <- zhou_components(f, use_Psi = TRUE, correct_denominators = TRUE)
  data.frame(ok = TRUE, tau2_se = f$Psi[1, 1], tau2_fpr = f$Psi[2, 2],
             rho = safe_rho_from_Psi(f$Psi), I2 = cc$I2,
             mu_se = as.numeric(f$coefficients[1]),
             mu_fpr = as.numeric(f$coefficients[2]))
}

jobs <- expand.grid(cell_id = design$cell_id, rep = seq_len(N_REP))
cat("Pseudo-true estimation:", nrow(jobs), "fits at K =", K_LARGE,
    "on", N_CORES, "cores\n")
t0 <- Sys.time()
raw <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  cell <- design[design$cell_id == jobs$cell_id[i], ]
  cbind(cell_id = jobs$cell_id[i], rep = jobs$rep[i],
        one_fit(cell, jobs$rep[i]))
}, mc.cores = N_CORES)
raw <- do.call(rbind, raw)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    "min\n")

mcse <- function(x) if (sum(is.finite(x)) > 1L)
  sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))) else NA_real_
agg <- do.call(rbind, lapply(design$cell_id, function(cid) {
  x <- raw[raw$cell_id == cid & raw$ok, ]
  cell <- design[design$cell_id == cid, ]
  pop <- population_I2(cell$rho_true, cell$nD, cell$nND)
  data.frame(cell, n_ok = nrow(x), K_large = K_LARGE,
             pop_I2 = unname(pop["pop_I2"]),
             tau2_se_gen = tau2_gen[1], tau2_fpr_gen = tau2_gen[2],
             pseudotrue_tau2_se = mean(x$tau2_se),
             pseudotrue_tau2_fpr = mean(x$tau2_fpr),
             pseudotrue_rho = mean(x$rho),
             pseudotrue_I2 = mean(x$I2),
             mcse_tau2_se = mcse(x$tau2_se),
             mcse_I2 = mcse(x$I2),
             ratio_tau2_se = mean(x$tau2_se) / tau2_gen[1],
             ratio_tau2_fpr = mean(x$tau2_fpr) / tau2_gen[2],
             approximation_bias_I2 = mean(x$I2) - unname(pop["pop_I2"]),
             row.names = NULL)
}))
write.csv(raw, file.path(outdir, "D1_pseudotrue_raw_fits.csv"), row.names = FALSE)
write.csv(agg, file.path(outdir, "D1_pseudotrue_values.csv"), row.names = FALSE)

cat("\nPSEUDO-TRUE VALUES (K =", K_LARGE, ",", N_REP, "large datasets per cell)\n")
print(agg[, c("scheme", "rho_true", "pop_I2", "pseudotrue_I2",
              "approximation_bias_I2", "ratio_tau2_se", "ratio_tau2_fpr",
              "pseudotrue_rho", "mcse_I2")],
      row.names = FALSE, digits = 4)
cat("\nIf ratio_tau2_* departs from 1 at this K, the fitting procedure is not",
    "\nconsistent for the generating variances and the Phase C 'bias' column",
    "\nmixes approximation bias with estimation bias.\n")
