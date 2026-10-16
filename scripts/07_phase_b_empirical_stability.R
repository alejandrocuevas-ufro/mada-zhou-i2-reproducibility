# Phase B — empirical stability of near-boundary Reitsma random-effects correlation
# Reproduction target: mada 0.5.12. Run from repository root.
#
# Examples:
#   PHASE_B_BOOTSTRAP=30 Rscript scripts/07_phase_b_empirical_stability.R
#   PHASE_B_BOOTSTRAP=400 Rscript scripts/07_phase_b_empirical_stability.R
#
# Analyses: REML/ML, continuity correction, leave-one-arm/family-out,
# and nonparametric arm/family block bootstrap. This is a sensitivity
# analysis, NOT a test that the true between-study rho equals +/-1.

options(stringsAsFactors = FALSE)
stopifnot(requireNamespace("mada", quietly = TRUE))
if (as.character(utils::packageVersion("mada")) != "0.5.12") {
  stop("Requires mada 0.5.12 for comparison with the archived results.")
}
library(mada)
source("scripts/zhou_i2_utils.R")

out_dir <- "phase_B_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seed <- 20261015L
set.seed(seed)
B <- suppressWarnings(as.integer(Sys.getenv("PHASE_B_BOOTSTRAP", "100")))
if (length(B) != 1 || is.na(B) || B < 0 || B > 2000) {
  stop("PHASE_B_BOOTSTRAP must be an integer from 0 to 2000.")
}
near_threshold <- 0.999
extreme_threshold <- 0.999999

all_data <- read.csv("data/published_reanalyses_2x2.csv",
                     fileEncoding = "UTF-8", check.names = FALSE)
dat <- all_data[
  all_data$article == "Pinto-Villalba 2026" &
    grepl("TCCS", all_data$group, fixed = TRUE) &
    grepl("isqu", all_data$group, ignore.case = TRUE),
]
row.names(dat) <- NULL
stopifnot(nrow(dat) == 8L, length(unique(dat$study_arm)) == 8L,
          all(rowSums(dat[, c("TP", "FP", "TN", "FN")]) == dat$N))
# "Postert family" is a CONSERVATIVE publication-family sensitivity.
# It does NOT assert overlapping participants in Postert 1998 and 1999.
dat$family <- ifelse(grepl("^Postert ", dat$study_arm),
                     "Postert family", dat$study_arm)
dat$nD <- dat$TP + dat$FN
dat$nND <- dat$FP + dat$TN
dat$sensitivity <- dat$TP / dat$nD
dat$specificity <- dat$TN / dat$nND
dat$has_zero <- rowSums(dat[, c("TP", "FP", "TN", "FN")] == 0) > 0
stopifnot(sum(dat$has_zero) == 2L, length(unique(dat$family)) == 7L)
write.csv(dat, file.path(out_dir, "phaseB_study_table.csv"), row.names = FALSE)

new_row <- function(experiment, method, correction_control, correction,
                    removed = "", bootstrap_unit = "", replicate = NA_integer_,
                    x = dat) {
  data.frame(
    experiment = experiment, method = method,
    correction_control = correction_control, correction = correction,
    removed = removed, bootstrap_unit = bootstrap_unit,
    replicate = replicate,
    n_rows = nrow(x), n_unique_arms = length(unique(x$study_arm)),
    n_unique_families = length(unique(x$family)),
    n_zero_arms = sum(x$has_zero),
    status = "not_run", warning_count = 0L, warnings = "", error = "",
    convergence_code = NA_integer_,
    tau2_se = NA_real_, tau2_fpr = NA_real_, rho_se_fpr = NA_real_,
    distance_to_boundary = NA_real_, corr_condition_number = NA_real_,
    det_between = NA_real_, expected_within_se = NA_real_,
    expected_within_fpr = NA_real_,
    i2_zd_corrected_percent = NA_real_,
    marginal_se_percent = NA_real_, marginal_fpr_percent = NA_real_,
    directional_max_percent = NA_real_, directional_min_percent = NA_real_,
    volume_total_exploratory_percent = NA_real_,
    stringsAsFactors = FALSE
  )
}

run_fit <- function(x, experiment, method = "reml",
                    correction_control = "all", correction = 0.5,
                    removed = "", bootstrap_unit = "",
                    replicate = NA_integer_) {
  ans <- new_row(experiment, method, correction_control, correction,
                 removed, bootstrap_unit, replicate, x)
  warnings_seen <- character()
  fitted <- tryCatch(
    withCallingHandlers(
      mada::reitsma(as_mada_data(x), method = method,
                    correction = correction,
                    correction.control = correction_control),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  ans$warning_count <- length(warnings_seen)
  ans$warnings <- paste(unique(warnings_seen), collapse = " | ")
  if (inherits(fitted, "error")) {
    ans$status <- "fit_error"
    ans$error <- conditionMessage(fitted)
    return(ans)
  }

  # Diagnostic convergence status, where exposed by mvmeta.fit.
  if (!is.null(fitted$opt) && is.list(fitted$opt) &&
      !is.null(fitted$opt$convergence)) {
    ans$convergence_code <- suppressWarnings(
      as.integer(fitted$opt$convergence)[1])
  }
  if (!is.na(ans$convergence_code) && ans$convergence_code != 0L) {
    ans$status <- "nonconverged"
    return(ans)
  }

  extraction <- tryCatch({
    Psi <- fitted$Psi
    if (!is.matrix(Psi) || !identical(dim(Psi), c(2L, 2L)) ||
        any(!is.finite(Psi))) stop("Invalid Psi matrix")
    v1 <- as.numeric(Psi[1, 1])
    v2 <- as.numeric(Psi[2, 2])
    if (v1 <= 0 || v2 <= 0) stop("Non-positive marginal variance")
    rho <- safe_rho_from_Psi(Psi)
    z <- zhou_components(fitted, use_Psi = TRUE,
                        correct_denominators = TRUE)
    w1 <- z$Evar_Se
    w2 <- z$Evar_Sp
    if (!all(is.finite(c(w1, w2, z$I2))) || any(c(w1, w2) <= 0)) {
      stop("Invalid expected within-study variance or I2")
    }
    # Whitened eigenvalues: retain meaningful PSD checks.
    G <- diag(c(1/sqrt(w1), 1/sqrt(w2))) %*% Psi %*%
      diag(c(1/sqrt(w1), 1/sqrt(w2)))
    G <- (G + t(G))/2
    ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-7 * max(1, max(abs(ev)))) {
      stop("Between-study covariance materially indefinite")
    }
    lam <- sort(pmax(ev, 0))
    directional <- lam/(1 + lam)
    list(v1 = v1, v2 = v2, rho = rho, w1 = w1, w2 = w2, z = z,
         directional = directional,
         volume = 1 - 1/sqrt(prod(1 + lam)))
  }, error = function(e) e)
  if (inherits(extraction, "error")) {
    ans$status <- "extraction_error"
    ans$error <- conditionMessage(extraction)
    return(ans)
  }

  v <- extraction
  ans$status <- "ok"
  ans$tau2_se <- v$v1
  ans$tau2_fpr <- v$v2
  ans$rho_se_fpr <- v$rho
  ans$distance_to_boundary <- 1 - abs(v$rho)
  ans$corr_condition_number <- if (abs(v$rho) == 1) Inf else
    (1 + abs(v$rho))/(1 - abs(v$rho))
  ans$det_between <- v$z$det_between
  ans$expected_within_se <- v$w1
  ans$expected_within_fpr <- v$w2
  ans$i2_zd_corrected_percent <- v$z$I2_percent
  ans$marginal_se_percent <- 100 * v$v1/(v$v1 + v$w1)
  ans$marginal_fpr_percent <- 100 * v$v2/(v$v2 + v$w2)
  ans$directional_max_percent <- 100 * max(v$directional)
  ans$directional_min_percent <- 100 * min(v$directional)
  ans$volume_total_exploratory_percent <- 100 * v$volume
  ans
}

# 1. Exact archived baseline fit (REML, 0.5, all), plus ML.
baseline <- do.call(rbind, lapply(c("reml", "ml"), function(m) {
  run_fit(dat, "baseline", method = m)
}))
write.csv(baseline, file.path(out_dir, "phaseB_baseline.csv"), row.names=FALSE)

# Protect reproducibility: compare the PRIMARY original fit to the archive.
arch <- read.csv("results/published_reanalyses_all_variants.csv",
                 check.names=FALSE)
arch <- arch[arch$dataset == "pinto_tccs_ischemic" &
             arch$implementation == "full Zhou-Dendukuri reconstruction", ]
stopifnot(nrow(arch) == 1, baseline$status[1] == "ok",
          abs(baseline$rho_se_fpr[1] - arch$rho) < 1e-5,
          abs(baseline$i2_zd_corrected_percent[1] -
                arch$I2_percent) < 0.02)

# 2. Zero-cell correction sensitivity: all/single; magnitudes.
sensitivity <- do.call(rbind, lapply(c("reml", "ml"), function(m) {
  do.call(rbind, lapply(c("all", "single"), function(mode) {
    do.call(rbind, lapply(c(.25, .5, 1), function(corr) {
      run_fit(dat, "continuity_correction", m, mode, corr)
    }))
  }))
}))
write.csv(sensitivity, file.path(out_dir, "phaseB_correction_sensitivity.csv"),
          row.names=FALSE)

# 3. Leave one analysis arm out, then leave one publication family out.
loo_arm <- do.call(rbind, lapply(unique(dat$study_arm), function(nm) {
  reduced <- dat[dat$study_arm != nm, ]
  do.call(rbind, lapply(c("reml", "ml"), function(m) {
    run_fit(reduced, "leave_one_arm_out", method=m, removed=nm)
  }))
}))
loo_family <- do.call(rbind, lapply(unique(dat$family), function(nm) {
  reduced <- dat[dat$family != nm, ]
  do.call(rbind, lapply(c("reml", "ml"), function(m) {
    run_fit(reduced, "leave_one_family_out", method=m, removed=nm)
  }))
}))
loo <- rbind(loo_arm, loo_family)
write.csv(loo, file.path(out_dir, "phaseB_leave_one_out.csv"),
          row.names=FALSE)

# 4. Nonparametric resampling of 8 arms vs 7 publication families.
# With replacement, re-sampling full 2x2 units preserves within-arm counts,
# but uncertainty intervals are only descriptive with k so small.
# The family scheme resamples all rows from Postert 1998/1999 together.
boot_rows <- list()
if (B > 0) {
  k <- 1L
  family_names <- unique(dat$family)
  for (rep in seq_len(B)) {
    idx_arm <- sample(seq_len(nrow(dat)), size=nrow(dat), replace=TRUE)
    family_sample <- sample(family_names, size=length(family_names),
                            replace=TRUE)
    idx_family <- unlist(lapply(family_sample, function(nm)
      which(dat$family == nm)), use.names=FALSE)
    for (unit in c("arm", "publication_family")) {
      x <- dat[if (unit == "arm") idx_arm else idx_family, ]
      row.names(x) <- NULL
      for (method in c("reml", "ml")) {
        boot_rows[[k]] <- run_fit(x, "bootstrap", method=method,
                                   bootstrap_unit=unit, replicate=rep)
        k <- k + 1L
      }
    }
    if (rep %% 25L == 0L) {
      cat("Completed bootstrap draw", rep, "of", B, "\n")
      flush.console()
    }
  }
  bootstrap <- do.call(rbind, boot_rows)
} else {
  bootstrap <- new_row("bootstrap", "reml", "all", .5)[FALSE, ]
}
write.csv(bootstrap, file.path(out_dir, "phaseB_bootstrap_replicates.csv"),
          row.names=FALSE)

q <- function(x, prob) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs=prob, na.rm=TRUE, names=FALSE))
}
bootstrap_summary <- list()
for (unit in c("arm", "publication_family")) {
  for (method in c("reml", "ml")) {
    x <- bootstrap[bootstrap$bootstrap_unit == unit &
                   bootstrap$method == method, ]
    ok <- x[x$status == "ok", ]
    n_ok <- nrow(ok)
    bootstrap_summary[[paste(unit, method)]] <- data.frame(
      unit=unit, method=method,
      n_requested=B, n_attempted=nrow(x), n_success=n_ok,
      n_failed=nrow(x)-n_ok,
      success_fraction=if(nrow(x)) n_ok/nrow(x) else NA_real_,
      boundary_abs_rho_ge_0999_n=sum(abs(ok$rho_se_fpr)>=near_threshold),
      boundary_abs_rho_ge_0999_fraction_success=
        if(n_ok) mean(abs(ok$rho_se_fpr)>=near_threshold) else NA_real_,
      boundary_abs_rho_ge_0999999_fraction_success=
        if(n_ok) mean(abs(ok$rho_se_fpr)>=extreme_threshold) else NA_real_,
      rho_p025=q(ok$rho_se_fpr,.025),
      rho_median=q(ok$rho_se_fpr,.5),
      rho_p975=q(ok$rho_se_fpr,.975),
      I2_p025=q(ok$i2_zd_corrected_percent,.025),
      I2_median=q(ok$i2_zd_corrected_percent,.5),
      I2_p975=q(ok$i2_zd_corrected_percent,.975),
      tau2_se_median=q(ok$tau2_se,.5),
      tau2_fpr_median=q(ok$tau2_fpr,.5),
      median_distinct_arms=q(x$n_unique_arms,.5),
      stringsAsFactors=FALSE
    )
  }
}
bootstrap_summary <- do.call(rbind, bootstrap_summary)
row.names(bootstrap_summary) <- NULL
write.csv(bootstrap_summary,
          file.path(out_dir, "phaseB_bootstrap_summary.csv"), row.names=FALSE)

# Additional model diagnostics plots (base R)
png(file.path(out_dir, "phaseB_loo_correlation.png"),
    width=2200, height=1400, res=240)
chosen <- loo[loo$method == "reml" &
              loo$experiment == "leave_one_arm_out", ]
yy <- ifelse(chosen$status=="ok", chosen$rho_se_fpr, NA_real_)
plot(seq_along(yy), yy, pch=19, xlim=c(.6,length(yy)+.4),
     ylim=c(-1.05,1.05), xaxt="n",
     xlab="Excluded analysis arm", ylab="Estimated correlation (Se, FPR)",
     main="Leave-one-arm-out Reitsma REML correlation")
abline(h=c(-near_threshold,near_threshold), lty=2)
axis(1, at=seq_along(yy), labels=chosen$removed, las=2, cex.axis=.7)
dev.off()

png(file.path(out_dir, "phaseB_bootstrap_correlation.png"),
    width=2000, height=1350, res=220)
chosen <- bootstrap[bootstrap$method=="reml" &
                    bootstrap$bootstrap_unit=="arm" &
                    bootstrap$status=="ok", ]
rr <- chosen$rho_se_fpr
if (length(rr)>0) {
  hist(rr, breaks=seq(-1,1,length.out=31), xlim=c(-1,1),
       main="Nonparametric arm bootstrap: REML correlation",
       xlab="Estimated random-effects correlation (Se, FPR)",
       ylab="Successful refits")
  abline(v=c(-near_threshold,near_threshold), lty=2)
} else {
  plot.new()
  text(.5,.5,"No successful bootstrap refits")
}
dev.off()

capture.output(
  { cat("Seed:",seed,"\nBootstrap requested:",B,"\n")
    cat("Data hash:",unname(tools::md5sum("data/published_reanalyses_2x2.csv")),"\n")
    print(sessionInfo())
  },
  file=file.path(out_dir,"sessionInfo.txt")
)
cat("\nPHASE B BASELINE (I2 reported in percent)\n")
print(baseline[, c("method","status","rho_se_fpr","tau2_se","tau2_fpr",
                   "i2_zd_corrected_percent","warning_count")], row.names=FALSE)
cat("\nPHASE B CORRECTION SENSITIVITY\n")
print(sensitivity[, c("method","correction_control","correction",
                      "status","rho_se_fpr","i2_zd_corrected_percent")],
      row.names=FALSE)
cat("\nPHASE B LEAVE-ONE-OUT\n")
print(loo[,c("experiment","method","removed","n_rows","status",
             "rho_se_fpr","i2_zd_corrected_percent")], row.names=FALSE)
cat("\nPHASE B BOOTSTRAP SUMMARY\n")
print(bootstrap_summary, row.names=FALSE)
cat("\nPhase B completed. All outputs in", out_dir, "\n")
