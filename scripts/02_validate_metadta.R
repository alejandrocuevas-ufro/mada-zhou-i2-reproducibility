# Reproduce the Zhou-Dendukuri calculation described in the metadta tutorial.
# Published model parameters are rounded, so this is an algebraic/formula check.

d <- read.csv("data/metadta_telomerase_reconstruction.csv", check.names = FALSE)

theta_se <- 1.19
theta_sp <- 2.34
tau2_se <- 0.18
tau2_sp <- 3.32

n_se <- d$TP + d$FN
n_sp <- d$TN + d$FP

Evar_se <- (
  exp(tau2_se / 2 + theta_se) +
    exp(tau2_se / 2 - theta_se) + 2
) * mean(1 / n_se)

Evar_sp <- (
  exp(tau2_sp / 2 + theta_sp) +
    exp(tau2_sp / 2 - theta_sp) + 2
) * mean(1 / n_sp)

I2_se <- tau2_se / (Evar_se + tau2_se)
I2_sp <- tau2_sp / (Evar_sp + tau2_sp)

I2_biv_published <- 0.02 / 100
det_E <- Evar_se * Evar_sp
sqrt_det_T <- I2_biv_published / (1 - I2_biv_published) * sqrt(det_E)
det_T <- sqrt_det_T^2
rho2 <- 1 - det_T / (tau2_se * tau2_sp)
rho_implied <- -sqrt(rho2)

out <- data.frame(
  Evar_Se = Evar_se,
  Evar_Sp = Evar_sp,
  I2_Se_percent = 100 * I2_se,
  I2_Sp_percent = 100 * I2_sp,
  published_bivariate_I2_percent = 100 * I2_biv_published,
  implied_det_between = det_T,
  implied_rho = rho_implied
)
write.csv(out, "results/metadta_formula_validation.csv", row.names = FALSE)
cat("metadta formula validation completed.\n")
