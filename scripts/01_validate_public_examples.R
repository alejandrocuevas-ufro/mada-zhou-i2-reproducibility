# Validate the implementation discrepancy using AuditC6 and AuditC.

library(mada)
source("scripts/zhou_i2_utils.R")

AuditC6 <- data.frame(
  TP = c(47, 126, 19, 36, 130, 84),
  FN = c(9, 51, 10, 3, 19, 2),
  FP = c(101, 272, 12, 78, 211, 68),
  TN = c(738, 1543, 192, 276, 959, 89)
)

fit_AuditC6 <- reitsma(AuditC6)
res_AuditC6 <- audit_zhou(fit_AuditC6, "AuditC6")

data("AuditC", package = "mada")
fit_AuditC <- reitsma(AuditC)
res_AuditC <- audit_zhou(fit_AuditC, "AuditC")

results <- rbind(res_AuditC6, res_AuditC)
diagnostics <- rbind(
  model_diagnostics(fit_AuditC6, "AuditC6"),
  model_diagnostics(fit_AuditC, "AuditC")
)

write.csv(results, "results/public_example_validation.csv", row.names = FALSE)
write.csv(diagnostics, "results/public_example_model_diagnostics.csv", row.names = FALSE)

cat("Public example validation completed.\n")
