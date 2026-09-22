# Reanalyse published applications of mada/Zhou-Dendukuri using the same
# Reitsma fit and replacing only the post-estimation I2 calculation.

library(mada)
source("scripts/zhou_i2_utils.R")

d <- read.csv(
  "data/published_reanalyses_2x2.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

subset_analysis <- function(id) {
  if (id == "dragota_primary") return(d[d$article == "Dragota 2026", ])
  if (id == "manetas_ppg") return(d[d$article == "Manetas-Stavrakakis 2023" & d$group == "PPG", ])
  if (id == "manetas_ecg") return(d[d$article == "Manetas-Stavrakakis 2023" & d$group == "ECG de una derivación", ])
  if (id == "pinto_global") return(d[d$article == "Pinto-Villalba 2026", ])
  if (id == "pinto_tccs_ischemic") {
    x <- d[d$article == "Pinto-Villalba 2026", ]
    return(x[grepl("isqu", x$group, ignore.case = TRUE) & grepl("TCCS", x$group), ])
  }
  if (id == "pinto_hemorrhagic") {
    x <- d[d$article == "Pinto-Villalba 2026", ]
    return(x[grepl("hemorr", x$group, ignore.case = TRUE), ])
  }
  stop("Unknown analysis id: ", id)
}

published_I2 <- c(
  dragota_primary = 71.8,
  manetas_ppg = 12.5,
  manetas_ecg = 9.2,
  pinto_global = 14.4,
  pinto_tccs_ischemic = 0.0,
  pinto_hemorrhagic = 13.1
)

analysis_ids <- names(published_I2)
all_variants <- list()
summary_rows <- list()

for (id in analysis_ids) {
  x <- subset_analysis(id)

  fit <- reitsma(
    as_mada_data(x),
    method = "reml",
    correction = 0.5,
    correction.control = "all"
  )

  aud <- audit_zhou(fit, id)
  diag <- model_diagnostics(fit, id)

  native <- aud$I2_percent[aud$implementation == "mada native / exact reconstruction"]
  psi_only <- aud$I2_percent[aud$implementation == "Psi only; mada denominators retained"]
  denom_only <- aud$I2_percent[aud$implementation == "mada SE^2; corrected denominators only"]
  corrected <- aud$I2_percent[aud$implementation == "full Zhou-Dendukuri reconstruction"]

  summary_rows[[id]] <- data.frame(
    analysis_id = id,
    k = nrow(x),
    published_I2_percent = unname(published_I2[id]),
    mada_reproduced_I2_percent = native,
    Psi_only_I2_percent = psi_only,
    denominators_only_I2_percent = denom_only,
    corrected_I2_percent = corrected,
    change_pp = corrected - native,
    pooled_sensitivity = diag$pooled_sensitivity,
    pooled_specificity = diag$pooled_specificity,
    var_mu_Se = diag$var_mu_Se,
    var_mu_FPR = diag$var_mu_FPR,
    tau2_Se = diag$tau2_Se,
    tau2_FPR = diag$tau2_FPR,
    rho = diag$rho_Se_FPR,
    row.names = NULL
  )

  all_variants[[id]] <- aud
}

summary_out <- do.call(rbind, summary_rows)
variants_out <- do.call(rbind, all_variants)
rownames(summary_out) <- NULL
rownames(variants_out) <- NULL

write.csv(summary_out, "results/published_reanalyses_summary.csv", row.names = FALSE)
write.csv(variants_out, "results/published_reanalyses_all_variants.csv", row.names = FALSE)

delta <- abs(summary_out$mada_reproduced_I2_percent - summary_out$published_I2_percent)
check <- data.frame(
  analysis_id = summary_out$analysis_id,
  published_I2_percent = summary_out$published_I2_percent,
  reproduced_I2_percent = summary_out$mada_reproduced_I2_percent,
  absolute_difference_pp = delta,
  within_0.15_pp = delta <= 0.15
)
write.csv(check, "results/published_reproduction_check.csv", row.names = FALSE)

if (!all(check$within_0.15_pp)) {
  warning("At least one reproduced published I2 differs by more than 0.15 percentage points.")
}

print(summary_out, digits = 6)
cat("Published reanalyses completed.\n")
