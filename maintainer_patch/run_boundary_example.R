## Boundary example for the Zhou-Dendukuri determinant.
##
## The purpose is not to redefine the statistic. It checks that an estimated
## covariance matrix on the |rho| = 1 boundary does not produce NaN solely
## because floating-point arithmetic returns a tiny negative determinant.

library(mada)

tccs <- data.frame(
  TP = c(20, 25, 11, 44, 13, 22, 9, 11),
  FN = c(5, 11, 0, 23, 5, 0, 1, 3),
  FP = c(2, 2, 0, 3, 8, 2, 1, 1),
  TN = c(18, 10, 9, 4, 3, 7, 75, 71)
)

fit <- reitsma(
  tccs,
  correction = 0.5,
  correction.control = "all"
)

rho <- fit$Psi[1,2] / sqrt(fit$Psi[1,1] * fit$Psi[2,2])

cat("tau2 sensitivity:", fit$Psi[1,1], "\n")
cat("tau2 FPR:", fit$Psi[2,2], "\n")
cat("rho:", format(rho, digits = 17), "\n")
cat("det(Psi):", format(det(fit$Psi), digits = 17), "\n")
cat("Zhou-Dendukuri I2 (%):",
    format(100 * summary(fit)$i2$Zhou, digits = 12), "\n")

stopifnot(is.finite(summary(fit)$i2$Zhou))
