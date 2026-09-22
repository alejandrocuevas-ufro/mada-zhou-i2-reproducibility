# Monte Carlo demonstration of K-dependence in mada 0.5.12 Zhou I2.

library(mada)
library(MASS)
source("scripts/zhou_i2_utils.R")

set.seed(20260918)
K_values <- c(6, 10, 20, 40, 80, 160)
R <- 500

mu <- c(2.0997, -1.2637)
sd_se <- 1.1746
sd_fpr <- 0.6381
rho <- 0.8543

Psi_population <- matrix(c(
  sd_se^2, rho * sd_se * sd_fpr,
  rho * sd_se * sd_fpr, sd_fpr^2
), 2, 2, byrow = TRUE)

nD <- 150
nND <- 500

simulate_once <- function(K) {
  eta <- MASS::mvrnorm(K, mu = mu, Sigma = Psi_population)
  Se <- invlogit(eta[, 1])
  FPR <- invlogit(eta[, 2])

  TP <- rbinom(K, nD, Se)
  FN <- nD - TP
  FP <- rbinom(K, nND, FPR)
  TN <- nND - FP

  fit <- reitsma(
    data.frame(TP = TP, FN = FN, FP = FP, TN = TN),
    method = "reml",
    correction = 0.5,
    correction.control = "single"
  )

  a <- zhou_components(fit, use_Psi = FALSE, correct_denominators = FALSE)
  b <- zhou_components(fit, use_Psi = TRUE, correct_denominators = TRUE)
  sm <- summary(fit)

  c(
    mada = a$I2,
    corrected = b$I2,
    var_mu_Se = sm$coefficients[1, 2]^2,
    var_mu_FPR = sm$coefficients[2, 2]^2,
    tau2_Se = fit$Psi[1, 1],
    tau2_FPR = fit$Psi[2, 2]
  )
}

res_list <- list()
idx <- 1L
for (K in K_values) {
  for (r in seq_len(R)) {
    z <- simulate_once(K)
    res_list[[idx]] <- data.frame(
      K = K,
      replicate = r,
      mada = z["mada"],
      corrected = z["corrected"],
      var_mu_Se = z["var_mu_Se"],
      var_mu_FPR = z["var_mu_FPR"],
      tau2_Se = z["tau2_Se"],
      tau2_FPR = z["tau2_FPR"]
    )
    idx <- idx + 1L
  }
}
res <- do.call(rbind, res_list)
rownames(res) <- NULL

summary_table <- do.call(rbind, lapply(K_values, function(K) {
  x <- res[res$K == K, ]
  data.frame(
    K = K,
    mada_mean_percent = 100 * mean(x$mada),
    corrected_mean_percent = 100 * mean(x$corrected),
    mada_median_percent = 100 * median(x$mada),
    corrected_median_percent = 100 * median(x$corrected),
    mean_var_mu_Se = mean(x$var_mu_Se),
    mean_var_mu_FPR = mean(x$var_mu_FPR),
    mean_tau2_Se = mean(x$tau2_Se),
    mean_tau2_FPR = mean(x$tau2_FPR)
  )
}))

theta_se <- mu[1]
theta_sp <- -mu[2]
tau2_se <- Psi_population[1, 1]
tau2_sp <- Psi_population[2, 2]

ESe <- (exp(tau2_se / 2 + theta_se) + exp(tau2_se / 2 - theta_se) + 2) / nD
ESp <- (exp(tau2_sp / 2 + theta_sp) + exp(tau2_sp / 2 - theta_sp) + 2) / nND
detT <- (1 - rho^2) * tau2_se * tau2_sp
I2_population <- sqrt(detT) / (sqrt(ESe * ESp) + sqrt(detT))

summary_table$population_I2_percent <- 100 * I2_population

write.csv(res, "results/simulation_replicates.csv", row.names = FALSE)
write.csv(summary_table, "results/simulation_summary.csv", row.names = FALSE)

png("results/Figure1_simulation.png", width = 2400, height = 1650, res = 300)
plot(
  summary_table$K, summary_table$corrected_mean_percent,
  type = "b", log = "x", ylim = c(0, 100),
  xlab = "Number of studies (K)", ylab = expression(I^2~"(%)"),
  xaxt = "n", pch = 16
)
axis(1, at = K_values, labels = K_values)
lines(summary_table$K, summary_table$mada_mean_percent, type = "b", pch = 1, lty = 2)
abline(h = 100 * I2_population, lty = 3)
legend(
  "bottomleft",
  legend = c("Zhou-Dendukuri reconstruction", "mada 0.5.12 implementation", "Population I2"),
  lty = c(1, 2, 3), pch = c(16, 1, NA), bty = "n"
)
dev.off()

pdf("results/Figure1_simulation.pdf", width = 8, height = 5.5)
plot(
  summary_table$K, summary_table$corrected_mean_percent,
  type = "b", log = "x", ylim = c(0, 100),
  xlab = "Number of studies (K)", ylab = expression(I^2~"(%)"),
  xaxt = "n", pch = 16
)
axis(1, at = K_values, labels = K_values)
lines(summary_table$K, summary_table$mada_mean_percent, type = "b", pch = 1, lty = 2)
abline(h = 100 * I2_population, lty = 3)
legend(
  "bottomleft",
  legend = c("Zhou-Dendukuri reconstruction", "mada 0.5.12 implementation", "Population I2"),
  lty = c(1, 2, 3), pch = c(16, 1, NA), bty = "n"
)
dev.off()

print(summary_table, digits = 6)
cat("Population Zhou-Dendukuri I2:", 100 * I2_population, "%\n")
cat("Simulation completed.\n")
