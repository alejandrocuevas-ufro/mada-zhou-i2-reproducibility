## Regression tests for the Zhou-Dendukuri bivariate I2 calculation.
## Intended for mada/tests/zhou_i2-regression.R.
##
## These tests use base R only, avoiding a new testthat dependency.

library(mada)

tol <- 1e-7

expect_close <- function(observed, expected, tolerance = tol, label = "") {
  if (!isTRUE(all.equal(as.numeric(observed), as.numeric(expected),
                        tolerance = tolerance))) {
    stop(sprintf("%s: observed %.12f, expected %.12f",
                 label, observed, expected))
  }
}

safe_det <- function(Psi, tolerance = sqrt(.Machine$double.eps)) {
  d <- det(Psi)
  scale <- prod(diag(Psi))
  threshold <- tolerance * max(1, abs(scale))
  if (d < 0 && abs(d) <= threshold) d <- 0
  if (d < 0) stop("Negative determinant beyond numerical tolerance.")
  d
}

zhou_reference <- function(fit) {
  stopifnot(inherits(fit, "reitsma"))
  stopifnot(length(fit$coefficients) == 2L)
  stopifnot(isTRUE(all.equal(fit$alphasens, 1)))
  stopifnot(isTRUE(all.equal(fit$alphafpr, 1)))

  theta_a <- as.numeric(fit$coefficients[1])
  theta_b <- -as.numeric(fit$coefficients[2])

  tau_a <- fit$Psi[1, 1]
  tau_b <- fit$Psi[2, 2]

  d <- fit$freqdata
  Nia <- 1 / (d$TP + d$FN)
  Nib <- 1 / (d$TN + d$FP)

  E_a <- (exp(tau_a / 2 + theta_a) +
          exp(tau_a / 2 - theta_a) + 2) * mean(Nia)
  E_b <- (exp(tau_b / 2 + theta_b) +
          exp(tau_b / 2 - theta_b) + 2) * mean(Nib)

  det_T <- safe_det(fit$Psi)
  det_E <- E_a * E_b

  sqrt(det_T) / (sqrt(det_E) + sqrt(det_T))
}

## AuditC6 is reproduced from the mada vignette.
AuditC6 <- data.frame(
  TP = c(47, 126, 19, 36, 130, 84),
  FN = c(9, 51, 10, 3, 19, 2),
  FP = c(101, 272, 12, 78, 211, 68),
  TN = c(738, 1543, 192, 276, 959, 89)
)

fit6 <- reitsma(AuditC6)
coef6_before <- fit6$coefficients
Psi6_before <- fit6$Psi
s6 <- summary(fit6)

expect_close(
  s6$i2$Zhou,
  0.804936101328665,
  label = "AuditC6 corrected Zhou-Dendukuri I2"
)
expect_close(
  s6$i2$Zhou,
  zhou_reference(fit6),
  label = "AuditC6 independent reference"
)

## Calling summary must not alter the fitted model.
stopifnot(identical(fit6$coefficients, coef6_before))
stopifnot(identical(fit6$Psi, Psi6_before))

## AuditC is distributed with mada.
data("AuditC")
fitC <- reitsma(AuditC)
coefC_before <- fitC$coefficients
PsiC_before <- fitC$Psi
sC <- summary(fitC)

expect_close(
  sC$i2$Zhou,
  0.840890186812083,
  label = "AuditC corrected Zhou-Dendukuri I2"
)
expect_close(
  sC$i2$Zhou,
  zhou_reference(fitC),
  label = "AuditC independent reference"
)

stopifnot(identical(fitC$coefficients, coefC_before))
stopifnot(identical(fitC$Psi, PsiC_before))

## Numerical-boundary example: TCCS ischaemic subgroup.
## The fitted correlation is extremely close to -1. The statistic should
## remain finite and close to zero rather than becoming NaN.
tccs <- data.frame(
  TP = c(20, 25, 11, 44, 13, 22, 9, 11),
  FN = c(5, 11, 0, 23, 5, 0, 1, 3),
  FP = c(2, 2, 0, 3, 8, 2, 1, 1),
  TN = c(18, 10, 9, 4, 3, 7, 75, 71)
)

fit_boundary <- reitsma(
  tccs,
  correction = 0.5,
  correction.control = "all"
)
i2_boundary <- summary(fit_boundary)$i2$Zhou

stopifnot(is.finite(i2_boundary))
stopifnot(i2_boundary >= 0)
stopifnot(i2_boundary < 1e-4)

message("All Zhou-Dendukuri I2 regression tests passed.")
