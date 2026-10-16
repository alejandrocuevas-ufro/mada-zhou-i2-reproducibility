# Regression and numerical tests for Phase C, run in GitHub Actions.
source("scripts/09_phase_c_known_rho_simulation.R")
near <- function(x,y,tol=1e-9)
  isTRUE(all.equal(unname(x),unname(y),tolerance=tol))
stopifnot(nrow(design)==11L*6L*2L,
          setequal(K_grid,c(6L,10L,20L,40L,80L,160L)),
          length(rho_grid)==11L,
          all(design$pop_marginal_se>0),
          all(design$pop_marginal_fpr>0),
          all(design$pop_I2_ZD>=0),
          all(design$pop_I2_ZD<=1),
          all(sums$attempted==R),
          all(sums$successful+sums$failed==R),
          all(sums$success_fraction>=0 & sums$success_fraction<=1))
for (s in sample_schemes$scheme) {
  d <- subset(design,scheme==s & K==20)
  stopifnot(near(d$pop_I2_ZD[d$rho_true==1],0),
            near(d$pop_I2_ZD[d$rho_true== -1],0),
            near(d$pop_I2_ZD[d$rho_true==0],
                 d$pop_I2_ZD[d$rho_true==0][1]))
  for(r in c(.5,.9,.99,.9999,1)) {
    stopifnot(near(d$pop_I2_ZD[d$rho_true==r],
                   d$pop_I2_ZD[d$rho_true== -r],tol=1e-10))
  }
  stopifnot(near(d$pop_marginal_se[1],d$pop_marginal_se[nrow(d)]),
            near(d$pop_marginal_fpr[1],d$pop_marginal_fpr[nrow(d)]))
}
# True between-study matrices are held identical in marginal variances;
# only rho, K and study denominators vary.
good <- results[results$status=="ok",]
stopifnot(nrow(good) > 0L, sum(sums$successful)>0L)
if(nrow(good)>0) {
  stopifnot(all(is.finite(good$rho_hat)),
            all(abs(good$rho_hat)<=1),
            all(good$I2_hat>=0 & good$I2_hat<=1),
            all(good$I2_marginal_se_hat>=0 &
                good$I2_marginal_se_hat<=1),
            all(good$I2_marginal_fpr_hat>=0 &
                good$I2_marginal_fpr_hat<=1),
            all(good$volume_total_hat>=-1e-9 &
                good$volume_total_hat<=1+1e-9))
}
stopifnot(nrow(results)==R*(nrow(design)+18L),
          all(sums$method %in% c("reml","ml")),
          nrow(sums)==nrow(design)+18L)
outputs <- file.path(outdir,c("phaseC_design_population.csv",
                              "phaseC_replicate_results.csv",
                              "phaseC_scenario_summary.csv",
                              "phaseC_nonextreme_discordance.csv",
                              "phaseC_boundary_frequency.png",
                              "phaseC_discordant_I2.png",
                              "sessionInfo.txt"))
stopifnot(all(file.exists(outputs)),all(file.info(outputs)$size>0))
cat("PASS: Phase C population formula, symmetry, variance invariance,",
    "paired targeted ML, all scenario counts and output files.\n")
