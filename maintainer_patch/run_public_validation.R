## Public, self-contained validation for the proposed Zhou-Dendukuri I2 patch.
##
## Run after installing the patched mada source:
##   Rscript maintainer_patch/run_public_validation.R

library(mada)

AuditC6 <- data.frame(
  TP = c(47, 126, 19, 36, 130, 84),
  FN = c(9, 51, 10, 3, 19, 2),
  FP = c(101, 272, 12, 78, 211, 68),
  TN = c(738, 1543, 192, 276, 959, 89)
)

data("AuditC")

fit6 <- reitsma(AuditC6)
fitC <- reitsma(AuditC)

out <- data.frame(
  dataset = c("AuditC6", "AuditC"),
  expected_mada_0_5_12_percent = c(
    52.5234082261936,
    40.3835194520243
  ),
  expected_patched_percent = c(
    80.4936101328665,
    84.0890186812083
  ),
  observed_percent = 100 * c(
    summary(fit6)$i2$Zhou,
    summary(fitC)$i2$Zhou
  )
)

print(out, row.names = FALSE)

stopifnot(
  isTRUE(all.equal(out$observed_percent,
                   out$expected_patched_percent,
                   tolerance = 1e-7))
)

message("Public validation passed.")
