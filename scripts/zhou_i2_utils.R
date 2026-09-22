# Utility functions for auditing the Zhou-Dendukuri bivariate I2
# in mada::summary.reitsma().
# Scope: intercept-only Reitsma model with logit sensitivity and logit FPR.
#
# Numerical note:
# Random-effects covariance estimates can lie on the boundary |rho| = 1.
# Floating-point arithmetic may then return values such as
# rho = 1.0000000000000002 and a tiny negative determinant.  The helper
# functions below clamp only deviations that are within numerical tolerance;
# larger violations stop execution.

invlogit <- function(x) 1 / (1 + exp(-x))

check_reitsma_scope <- function(fit) {
  if (!inherits(fit, "reitsma")) stop("fit must inherit from 'reitsma'.")
  if (!isTRUE(all.equal(fit$alphasens, 1)) || !isTRUE(all.equal(fit$alphafpr, 1))) {
    stop("Only logit-logit Reitsma models (alpha=1) are supported.")
  }
  if (length(fit$coefficients) != 2L) {
    stop("Only the intercept-only Reitsma model is supported.")
  }
}

safe_rho_from_Psi <- function(Psi, tol = sqrt(.Machine$double.eps)) {
  rho_raw <- Psi[1, 2] / sqrt(Psi[1, 1] * Psi[2, 2])
  if (!is.finite(rho_raw)) stop("Non-finite random-effects correlation derived from Psi.")
  if (abs(rho_raw) > 1 + tol) {
    stop(sprintf("Estimated correlation outside [-1,1] beyond numerical tolerance: %.17g", rho_raw))
  }
  max(-1, min(1, rho_raw))
}

safe_nonnegative_determinant <- function(x,
                                         scale = 1,
                                         tol = sqrt(.Machine$double.eps)) {
  if (!is.finite(x)) stop("Non-finite determinant encountered.")
  threshold <- tol * max(1, abs(scale))
  if (x < 0 && abs(x) <= threshold) return(0)
  if (x < 0) {
    stop(sprintf("Negative determinant beyond numerical tolerance: %.17g", x))
  }
  x
}

zhou_components <- function(fit,
                            use_Psi = FALSE,
                            correct_denominators = FALSE) {
  check_reitsma_scope(fit)

  sm <- summary(fit)
  beta <- as.numeric(fit$coefficients)
  theta_se <- beta[1]
  theta_sp <- -beta[2]

  Psi <- fit$Psi
  rho <- safe_rho_from_Psi(Psi)

  if (use_Psi) {
    tau2_se <- Psi[1, 1]
    tau2_sp <- Psi[2, 2]
  } else {
    tau2_se <- sm$coefficients[1, 2]^2
    tau2_sp <- sm$coefficients[2, 2]^2
  }

  d <- fit$freqdata

  if (correct_denominators) {
    inv_n_se <- 1 / (d$TP + d$FN)
    inv_n_sp <- 1 / (d$TN + d$FP)
  } else {
    inv_n_se <- 1 / (d$TP + d$FP)
    inv_n_sp <- 1 / (d$TN + d$FN)
  }

  Evar_se <- (
    exp(theta_se + tau2_se / 2) +
      exp(-theta_se + tau2_se / 2) + 2
  ) * mean(inv_n_se)

  Evar_sp <- (
    exp(theta_sp + tau2_sp / 2) +
      exp(-theta_sp + tau2_sp / 2) + 2
  ) * mean(inv_n_sp)

  det_between_raw <- (1 - rho^2) * tau2_se * tau2_sp
  det_between <- safe_nonnegative_determinant(
    det_between_raw,
    scale = tau2_se * tau2_sp
  )

  det_within <- Evar_se * Evar_sp
  det_within <- safe_nonnegative_determinant(
    det_within,
    scale = Evar_se * Evar_sp
  )

  sqrt_between <- sqrt(det_between)
  sqrt_within <- sqrt(det_within)

  denom <- sqrt_within + sqrt_between
  I2 <- if (denom == 0) 0 else sqrt_between / denom

  data.frame(
    theta_Se = theta_se,
    theta_Sp = theta_sp,
    tau2_Se = tau2_se,
    tau2_FPR = tau2_sp,
    rho = rho,
    Evar_Se = Evar_se,
    Evar_Sp = Evar_sp,
    det_between = det_between,
    det_expected_within = det_within,
    I2 = I2,
    I2_percent = 100 * I2,
    row.names = NULL
  )
}

audit_zhou <- function(fit, dataset_name) {
  native <- unname(summary(fit)$i2$Zhou[1])

  a <- zhou_components(fit, use_Psi = FALSE, correct_denominators = FALSE)
  b <- zhou_components(fit, use_Psi = TRUE,  correct_denominators = FALSE)
  c <- zhou_components(fit, use_Psi = FALSE, correct_denominators = TRUE)
  d <- zhou_components(fit, use_Psi = TRUE,  correct_denominators = TRUE)

  if (!isTRUE(all.equal(native, a$I2, tolerance = 1e-10))) {
    warning(sprintf("Exact reconstruction did not match native mada for %s", dataset_name))
  }

  out <- rbind(
    data.frame(dataset = dataset_name, implementation = "mada native / exact reconstruction", a),
    data.frame(dataset = dataset_name, implementation = "Psi only; mada denominators retained", b),
    data.frame(dataset = dataset_name, implementation = "mada SE^2; corrected denominators only", c),
    data.frame(dataset = dataset_name, implementation = "full Zhou-Dendukuri reconstruction", d)
  )
  rownames(out) <- NULL
  out
}

model_diagnostics <- function(fit, dataset_name) {
  sm <- summary(fit)
  Psi <- fit$Psi

  data.frame(
    dataset = dataset_name,
    k = nrow(fit$freqdata),
    pooled_sensitivity = invlogit(as.numeric(fit$coefficients[1])),
    pooled_specificity = 1 - invlogit(as.numeric(fit$coefficients[2])),
    mu_logit_Se = as.numeric(fit$coefficients[1]),
    mu_logit_FPR = as.numeric(fit$coefficients[2]),
    var_mu_Se = sm$coefficients[1, 2]^2,
    var_mu_FPR = sm$coefficients[2, 2]^2,
    tau2_Se = Psi[1, 1],
    tau2_FPR = Psi[2, 2],
    rho_Se_FPR = safe_rho_from_Psi(Psi),
    row.names = NULL
  )
}

as_mada_data <- function(d) {
  data.frame(TP = d$TP, FN = d$FN, FP = d$FP, TN = d$TN)
}
