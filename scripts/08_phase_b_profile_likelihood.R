# Phase B supplemental: profile ML and restricted likelihood in rho
# for the same normal-normal approximation used by mada::reitsma().
# No extra R packages besides mada needed. Must run from repository root.
#
# IMPORTANT: Relative objective only; 95% chi-squared cutoff is NOT a
# validated reference distribution at the singular covariance boundary.
# The goal is to assess *identifiability*, not construct a formal CI.

if (!exists("dat") || !exists("baseline")) {
  stop("Run scripts/07_phase_b_empirical_stability.R before this script.")
}

dat_c <- dat
# The exact two-stage transform in mada alpha=1 and correction.control=all:
stopifnot(any(dat_c$has_zero))
P <- as.matrix(dat_c[, c("TP","FN","FP","TN")]) + 0.5
Y <- cbind(log(P[,"TP"] / P[,"FN"]), log(P[,"FP"] / P[,"TN"]))
S <- cbind(1/P[,"TP"] + 1/P[,"FN"], 1/P[,"FP"] + 1/P[,"TN"])
colnames(Y) <- c("logit_Se", "logit_FPR")
colnames(S) <- c("var_logit_Se", "var_logit_FPR")
stopifnot(all(is.finite(Y)), all(S > 0))

profile_objective <- function(logvar, rho, likelihood = "ml") {
  # v1 and v2 are between-study variances; rho is correlation.
  v <- exp(logvar)
  c12 <- rho * sqrt(v[1] * v[2])
  a <- S[,1] + v[1]
  b <- S[,2] + v[2]
  d <- a*b - c12*c12
  if (any(!is.finite(d)) || any(d <= 0)) return(1e100)
  i11 <- b/d
  i12 <- -c12/d
  i22 <- a/d
  A <- matrix(c(sum(i11), sum(i12), sum(i12), sum(i22)), 2, 2)
  rhs <- c(sum(i11*Y[,1]+i12*Y[,2]),
           sum(i12*Y[,1]+i22*Y[,2]))
  beta <- tryCatch(solve(A, rhs), error=function(e) NA_real_)
  if (any(!is.finite(beta))) return(1e100)
  e1 <- Y[,1] - beta[1]
  e2 <- Y[,2] - beta[2]
  Q <- sum(i11*e1^2 + 2*i12*e1*e2 + i22*e2^2)
  obj <- sum(log(d)) + Q
  if (likelihood == "reml") {
    # REML adds log|sum(S_i+Psi)^(-1)|, up to a rho-independent term.
    logdetA <- determinant(A, logarithm=TRUE)
    if (logdetA$sign <= 0) return(1e100)
    obj <- obj + as.numeric(logdetA$modulus)
  }
  obj
}

# Optimise over non-negative marginal variances for each *fixed* rho.
# Multiple starts limit sensitivity to local optimizer solutions.
rho_grid <- sort(unique(c(-1,-0.999999,-.9999,-.999,-.99,-.95,
                          -.9,-.8,-.7,-.6,-.5,-.25,0,.25,.5,.75,.9,.99,1)))
records <- list()
idx <- 1L
for (method in c("ml","reml")) {
  ref <- baseline[baseline$method == method &
                  baseline$status == "ok",]
  stopifnot(nrow(ref)==1)
  starts <- rbind(
    log(c(ref$tau2_se,ref$tau2_fpr)),
    log(c(.1,2.4)), log(c(.01,.2)),
    log(c(1,1)), c(-16,log(1))
  )
  for (rho in rho_grid) {
    trials <- lapply(seq_len(nrow(starts)), function(i) {
      tryCatch(
        stats::optim(starts[i,], profile_objective, rho=rho, likelihood=method,
                     method="L-BFGS-B", lower=c(-18,-18), upper=c(8,8),
                     control=list(maxit=250)),
        error=function(e) e
      )
    })
    ok <- trials[vapply(trials, function(t)
      is.list(t) && !is.null(t$value) && is.finite(t$value), logical(1))]
    if (!length(ok)) stop("All profilers failed at rho=", rho,
                          " method=", method)
    best <- ok[[which.min(vapply(ok,function(t)t$value,numeric(1)))]]
    records[[idx]] <- data.frame(
      method=method, rho_fixed=rho, neg2logLik_up_to_constant=best$value,
      tau2_se=exp(best$par[1]), tau2_fpr=exp(best$par[2]),
      at_variance_lower_bound=any(best$par < -17.95),
      optimizer_convergence=best$convergence,
      successful_starts=length(ok), stringsAsFactors=FALSE
    )
    idx <- idx + 1L
  }
}
profile <- do.call(rbind,records)
row.names(profile) <- NULL
for (method in c("ml","reml")) {
  sel <- profile$method==method
  profile$delta_deviance[sel] <- profile$neg2logLik_up_to_constant[sel] -
    min(profile$neg2logLik_up_to_constant[sel])
  # Evaluate at the *reported* fitted covariance for a strict cross-check
  # that the function and the mada model use the same likelihood geometry.
  ref <- baseline[baseline$method==method, ]
  plug <- profile_objective(log(c(ref$tau2_se,ref$tau2_fpr)),
                            rho=ref$rho_se_fpr, likelihood=method)
  stopifnot(is.finite(plug),
            min(profile$neg2logLik_up_to_constant[sel]) <= plug + 1e-4,
            plug - min(profile$neg2logLik_up_to_constant[sel]) < .2)
}
write.csv(profile,file.path(out_dir,"phaseB_profile_likelihood.csv"),
          row.names=FALSE)

# Graph of *relative deviance* (not an inferential p-value).
png(file.path(out_dir,"phaseB_profile_likelihood.png"),
    width=2000,height=1400,res=240)
plot(c(-1,1), range(c(0,profile$delta_deviance,4)),
     type="n", xlab="Fixed between-study correlation (Se / FPR)",
     ylab="Profile difference in -2 log likelihood",
     main="Profile likelihood: small-k near-boundary covariance")
for (i in seq_along(c("ml","reml"))) {
  m <- c("ml","reml")[i]
  sub <- profile[profile$method==m,]
  lines(sub$rho_fixed,sub$delta_deviance, col=i,lty=i,lwd=2)
  points(sub$rho_fixed,sub$delta_deviance,col=i,pch=i)
}
abline(h=3.841459, lty=3, col="grey50")
legend("topleft",legend=c("ML","REML",
       "Nominal chi-squared 3.84 (NOT calibrated at boundary)"),
       col=c(1,2,"grey50"),lty=c(1,2,3),lwd=c(2,2,1),bty="n",cex=.68)
dev.off()

cat("\nPHASE B PROFILED LIKELIHOOD (relative deviance)\n")
print(profile[,c("method","rho_fixed","delta_deviance",
                "tau2_se","tau2_fpr","at_variance_lower_bound")],
      row.names=FALSE,digits=7)
cat("\nLargest profile difference ML:",
    max(subset(profile,method=="ml")$delta_deviance),
    "; REML:",max(subset(profile,method=="reml")$delta_deviance),"\n")
cat("Note: profile values are descriptive; chi-squared thresholds have",
    "nonstandard behaviour at a covariance/correlation boundary.\n")
