# Test the Phase B pipeline and archived-fit reproducibility.
# Run: PHASE_B_BOOTSTRAP=20 Rscript scripts/07_test_phase_b.R
source("scripts/07_phase_b_empirical_stability.R")

stopifnot(nrow(dat) == 8L,
          sum(dat$has_zero) == 2L,
          length(unique(dat$family)) == 7L,
          identical(sort(dat$study_arm),
                    sort(unique(dat$study_arm))),
          all(c("Postert 1998", "Postert 1999") %in% dat$study_arm),
          nrow(baseline) == 2L,
          baseline$status[baseline$method=="reml"] == "ok",
          nrow(sensitivity) == 12L,
          all(sensitivity$experiment == "continuity_correction"),
          length(unique(paste(sensitivity$method,
                              sensitivity$correction_control,
                              sensitivity$correction))) == 12L,
          nrow(loo_arm) == 8L * 2L,
          nrow(loo_family) == 7L * 2L,
          nrow(loo) == 30L,
          nrow(bootstrap) == B * 4L,
          nrow(bootstrap_summary) == 4L,
          all(bootstrap_summary$n_attempted == B),
          all(bootstrap_summary$n_success +
              bootstrap_summary$n_failed == B))
good <- rbind(baseline, sensitivity, loo, bootstrap)
good <- good[good$status=="ok", ]
if (nrow(good)) {
  stopifnot(all(is.finite(good$tau2_se)),
            all(is.finite(good$tau2_fpr)),
            all(good$tau2_se > 0),
            all(good$tau2_fpr > 0),
            all(abs(good$rho_se_fpr) <= 1),
            all(good$i2_zd_corrected_percent >= 0),
            all(good$i2_zd_corrected_percent <= 100),
            all(good$directional_min_percent >= 0),
            all(good$directional_max_percent <= 100))
}
stopifnot(baseline$i2_zd_corrected_percent[1] < 0.05,
          abs(baseline$rho_se_fpr[1]) > .999999,
          all(file.exists(file.path(out_dir,
            c("phaseB_study_table.csv", "phaseB_baseline.csv",
              "phaseB_correction_sensitivity.csv",
              "phaseB_leave_one_out.csv",
              "phaseB_bootstrap_replicates.csv",
              "phaseB_bootstrap_summary.csv",
              "phaseB_loo_correlation.png",
              "phaseB_bootstrap_correlation.png", "sessionInfo.txt")))))
cat("PASS: Phase B data structure, archived REML fit, ML/REML sensitivity,\n",
    "leave-one-out design and bootstrap accounting.\n")
