# Phase C: known-population-rho Monte Carlo for Zhou-Dendukuri I2.
# Run from repository root (mada 0.5.12):
#   PHASE_C_REPS=100 Rscript scripts/09_phase_c_known_rho_simulation.R
# A reproducible diagnostic simulation, not a general inferential validation.
options(stringsAsFactors=FALSE)
stopifnot(requireNamespace("mada",quietly=TRUE),
          as.character(packageVersion("mada"))=="0.5.12")
library(mada)
source("scripts/zhou_i2_utils.R")
R <- suppressWarnings(as.integer(Sys.getenv("PHASE_C_REPS","100")))
if(length(R)!=1 || is.na(R) || R<1 || R>1000) stop("PHASE_C_REPS requires 1..1000")
seed <- 20261015L
set.seed(seed)
outdir <- "phase_C_outputs"
dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
mu <- c(2.0997,-1.2637)  # logit sensitivity, logit FPR
SD <- c(1.1746,.6381)
tau2 <- SD^2

rho_grid <- c(-1,-.9999,-.99,-.9,-.5,0,.5,.9,.99,.9999,1)
K_grid <- c(6L,10L,20L,40L,80L,160L)
sample_schemes <- data.frame(scheme=c("low_precision","original_precision"),
                             nD=c(50L,150L),nND=c(50L,500L))
design <- expand.grid(rho_true=rho_grid,K=K_grid,
                      scheme=sample_schemes$scheme,
                      KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
design$nD <- sample_schemes$nD[match(design$scheme,sample_schemes$scheme)]
design$nND <- sample_schemes$nND[match(design$scheme,sample_schemes$scheme)]
design$scenario_id <- seq_len(nrow(design))

# Population expected within-study variances use the specified generating
# means and marginal random-effects variances; they do NOT depend on rho.
population <- function(rho,nD,nND) {
  w1 <- (exp(mu[1]+tau2[1]/2)+exp(-mu[1]+tau2[1]/2)+2)/nD
  theta_sp <- -mu[2]
  w2 <- (exp(theta_sp+tau2[2]/2)+exp(-theta_sp+tau2[2]/2)+2)/nND
  rootT <- prod(SD)*sqrt((1-abs(rho))*(1+abs(rho)))
  zd <- rootT/(sqrt(w1*w2)+rootT)
  # Standardized generalized eigenvalue diagnostics are exploratory,
  # distinct from the original Zhou-Dendukuri statistic.
  Psi <- matrix(c(tau2[1],rho*prod(SD),
                  rho*prod(SD),tau2[2]),2,2,byrow=TRUE)
  whitened <- diag(1/sqrt(c(w1,w2))) %*% Psi %*%
    diag(1/sqrt(c(w1,w2)))
  e <- sort(pmax(eigen(whitened,symmetric=TRUE,
                       only.values=TRUE)$values,0))
  c(pop_I2_ZD=zd,pop_marginal_se=tau2[1]/(tau2[1]+w1),
    pop_marginal_fpr=tau2[2]/(tau2[2]+w2),
    pop_direction_max=max(e/(1+e)),
    pop_volume_comparator=1-1/sqrt(prod(1+e)),
    pop_w_se=w1,pop_w_fpr=w2)
}
pm <- t(vapply(seq_len(nrow(design)), function(i) {
  d <- design[i,]
  population(d$rho_true,d$nD,d$nND)
},numeric(7)))
design <- cbind(design,as.data.frame(pm))
write.csv(design,file.path(outdir,"phaseC_design_population.csv"),
          row.names=FALSE)

simulate_counts <- function(K,rho,nD,nND) {
  # Exact PSD generator also supports rho = +/-1 without relying on
  # near-singular numerical matrix factorization.
  z1 <- rnorm(K); z2 <- rnorm(K)
  eta_se <- mu[1]+SD[1]*z1
  eta_fpr <- mu[2]+SD[2] *
    (rho*z1+sqrt((1-abs(rho))*(1+abs(rho)))*z2)
  prob_se <- plogis(eta_se); prob_fpr <- plogis(eta_fpr)
  TP <- rbinom(K,size=nD,prob=prob_se)
  FP <- rbinom(K,size=nND,prob=prob_fpr)
  data.frame(TP=TP,FN=nD-TP,FP=FP,TN=nND-FP)
}

estimate <- function(counts,method) {
  n_zero_arms <- sum(rowSums(counts==0)>0)
  warnings_seen <- character()
  fit <- tryCatch(withCallingHandlers(
    mada::reitsma(counts,method=method,correction=.5,
                  correction.control="all"),
    warning=function(w) {
      warnings_seen <<- c(warnings_seen,conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error=function(e)e)
  default <- data.frame(
    method=method,status="ok",
    warning_count=length(warnings_seen),
    warnings=paste(unique(warnings_seen),collapse=" | "),
    failure_reason="",n_zero_arms=n_zero_arms,
    tau2_se_hat=NA_real_,tau2_fpr_hat=NA_real_,rho_hat=NA_real_,
    I2_hat=NA_real_,I2_marginal_se_hat=NA_real_,
    I2_marginal_fpr_hat=NA_real_,direction_max_hat=NA_real_,
    volume_total_hat=NA_real_, stringsAsFactors=FALSE
  )
  if(inherits(fit,"error")) {
    default$status <- "fit_error"
    default$failure_reason <- conditionMessage(fit)
    return(default)
  }
  v <- tryCatch({
    P <- fit$Psi
    if(!is.matrix(P)||any(!is.finite(P))||
       any(diag(P)<=0))stop("invalid or nonpositive between-study variances")
    rho <- safe_rho_from_Psi(P)
    comp <- zhou_components(fit,use_Psi=TRUE,correct_denominators=TRUE)
    W <- c(comp$Evar_Se,comp$Evar_Sp)
    if(any(!is.finite(W))||any(W<=0))stop("invalid expected within variance")
    GG <- diag(1/sqrt(W)) %*% P %*% diag(1/sqrt(W))
    eigenval <- eigen((GG+t(GG))/2,symmetric=TRUE,
                      only.values=TRUE)$values
    if (min(eigenval) < (-1e-6 * max(1, max(abs(eigenval)))))
      stop("between-study covariance is materially indefinite")
    eigenval <- pmax(eigenval,0)
    list(rho=rho,P=P,W=W,comp=comp,
         direction=max(eigenval/(1+eigenval)),
         volume=1-1/sqrt(prod(1+eigenval)))
  },error=function(e)e)
  if(inherits(v,"error")) {
    default$status <- "extraction_error"
    default$failure_reason <- conditionMessage(v)
    return(default)
  }
  default$tau2_se_hat <- v$P[1,1]
  default$tau2_fpr_hat <- v$P[2,2]
  default$rho_hat <- v$rho
  default$I2_hat <- v$comp$I2
  default$I2_marginal_se_hat <- v$P[1,1]/(v$P[1,1]+v$W[1])
  default$I2_marginal_fpr_hat <- v$P[2,2]/(v$P[2,2]+v$W[2])
  default$direction_max_hat <- v$direction
  default$volume_total_hat <- v$volume
  default
}

# Main: REML in every scenario; ML uses the *same simulated counts*
# in targeted checks (K=6,20,80 and true rho=0,+/-.99) to limit cost.
run_ml <- function(d) d$K %in% c(6,20,80) &&
                            d$rho_true %in% c(-.99,0,.99)
result <- vector("list", nrow(design)*R)
pos <- 1L
for (i in seq_len(nrow(design))) {
  d <- design[i,]
  do_ml <- run_ml(d)
  for(rep in seq_len(R)) {
    counts <- simulate_counts(d$K,d$rho_true,d$nD,d$nND)
    methods <- if(do_ml) c("reml","ml") else "reml"
    fit_rows <- lapply(methods,function(m)estimate(counts,m))
    for(x in fit_rows) {
      for(nm in c("scenario_id","rho_true","K","scheme","nD","nND",
                  "pop_I2_ZD","pop_marginal_se","pop_marginal_fpr",
                  "pop_direction_max","pop_volume_comparator"))
        x[[nm]] <- d[[nm]]
      x$replicate <- rep
      result[[pos]] <- x
      pos <- pos+1L
    }
  }
  if(i%%10==0) {
    cat("Completed scenario",i,"/",nrow(design),"\n")
    flush.console()
  }
}
results <- do.call(rbind,result[seq_len(pos-1L)])
rownames(results) <- NULL
results$rho_near_boundary <- results$status=="ok" &
  abs(results$rho_hat)>=.999
results$rho_extreme_boundary <- results$status=="ok" &
  abs(results$rho_hat)>=.999999
results$low_estimated_I2 <- results$status=="ok" &
  results$I2_hat<.05
results$discordant_low_I2 <- results$status=="ok" &
  results$I2_hat<.05 & results$pop_I2_ZD>=.25
write.csv(results,file.path(outdir,"phaseC_replicate_results.csv"),
          row.names=FALSE)

quant <- function(x,p) {
  x <- x[is.finite(x)]
  if(!length(x))return(NA_real_)
  as.numeric(stats::quantile(x,probs=p,na.rm=TRUE,names=FALSE))
}
groups <- unique(results[,c("scenario_id","method")])
sum_rows <- lapply(seq_len(nrow(groups)),function(j){
  g <- groups[j,]
  x <- results[results$scenario_id==g$scenario_id &
               results$method==g$method,]
  ok <- x[x$status=="ok",]
  n <- nrow(x); m <- nrow(ok)
  rate <- function(v) if(m)mean(v) else NA_real_
  data.frame(
    scenario_id=g$scenario_id,method=g$method,
    rho_true=x$rho_true[1],K=x$K[1],scheme=x$scheme[1],
    nD=x$nD[1],nND=x$nND[1],
    pop_I2_ZD=x$pop_I2_ZD[1],
    pop_marginal_se=x$pop_marginal_se[1],
    pop_marginal_fpr=x$pop_marginal_fpr[1],
    pop_direction_max=x$pop_direction_max[1],
    attempted=n,successful=m,failed=n-m,
    success_fraction=m/n,
    warning_fraction=mean(x$warning_count>0),
    any_zero_arm_fraction=mean(x$n_zero_arms>0),
    boundary_fraction_success=rate(ok$rho_near_boundary),
    boundary_999999_fraction_success=rate(ok$rho_extreme_boundary),
    low_I2_fraction_success=rate(ok$low_estimated_I2),
    discordant_low_I2_fraction_success=rate(ok$discordant_low_I2),
    rho_median=quant(ok$rho_hat,.5),
    rho_p025=quant(ok$rho_hat,.025),
    rho_p975=quant(ok$rho_hat,.975),
    I2_median=quant(ok$I2_hat,.5),
    I2_p025=quant(ok$I2_hat,.025),
    I2_p975=quant(ok$I2_hat,.975),
    I2_bias_mean=if(m)mean(ok$I2_hat)-x$pop_I2_ZD[1] else NA_real_,
    I2_abs_error_mean=if(m)mean(abs(ok$I2_hat-x$pop_I2_ZD[1])) else NA_real_,
    tau2_se_mean=if(m)mean(ok$tau2_se_hat) else NA_real_,
    tau2_fpr_mean=if(m)mean(ok$tau2_fpr_hat) else NA_real_
  )
})
sums <- do.call(rbind,sum_rows)
rownames(sums) <- NULL
write.csv(sums,file.path(outdir,"phaseC_scenario_summary.csv"),
          row.names=FALSE)
write.csv(subset(sums,method=="reml" & abs(rho_true)<=.9 &
                   pop_I2_ZD>=.25),
          file.path(outdir,"phaseC_nonextreme_discordance.csv"),
          row.names=FALSE)

draw_plot <- function(value,ylabel,file) {
  png(file.path(outdir,file),width=2200,height=1450,res=220)
  par(mfrow=c(1,2),mar=c(5,5,3,2))
  for(scheme in sample_schemes$scheme) {
    # Use explicit indexing: subset(sums, scheme == scheme) silently
    # compares the data column to itself and mixes both sample regimes.
    d <- sums[sums$method == "reml" & sums$scheme == scheme, , drop=FALSE]
    stopifnot(nrow(d) == length(rho_grid)*length(K_grid),
              length(unique(d$scheme)) == 1L,
              length(unique(d$nD)) == 1L,
              length(unique(d$nND)) == 1L)
    plot(range(rho_grid),c(0,1),type="n",
         xlab="True between-study correlation",
         ylab=ylabel,main=paste(scheme,"(nD/nND:",
           paste(unique(d$nD),unique(d$nND),sep="/"),")"))
    ks <- sort(unique(d$K))
    for(k in ks) {
      z <- d[d$K==k,]
      z <- z[order(z$rho_true),]
      lines(z$rho_true,z[[value]],lty=match(k,ks),lwd=1.7)
      points(z$rho_true,z[[value]],pch=match(k,ks),cex=.55)
    }
    legend("top",legend=paste("K =",ks),lty=seq_along(ks),
           pch=seq_along(ks),bty="n",cex=.65,ncol=3)
  }
  dev.off()
}
draw_plot("boundary_fraction_success","Pr(|rho_hat| >= 0.999), successful fits",
          "phaseC_boundary_frequency.png")
draw_plot("discordant_low_I2_fraction_success",
          "Pr(I2_hat < 5% and population I2 >= 25%)",
          "phaseC_discordant_I2.png")

capture.output({
  cat("Phase C Monte Carlo\nSeed:",seed,"\nReplicates per scenario:",R,
      "\nTotal scenarios:",nrow(design),"\n")
  cat("Total attempted fits:",nrow(results),
      "\nSuccess:",sum(results$status=="ok"),
      "\nFailed:",sum(results$status!="ok"),"\n")
  cat("Package mada:",as.character(packageVersion("mada")),"\n")
  print(sessionInfo())
},file=file.path(outdir,"sessionInfo.txt"))

cat("\nPHASE C DIAGNOSTIC STATUS\n")
print(table(results$status,useNA="ifany"))
cat("Failure reasons (top entries):\n")
print(utils::head(sort(table(results$failure_reason[results$status!="ok"]),
                       decreasing=TRUE),10))
cat("\nPHASE C TOTAL: scenarios=",nrow(design),
    " reps=",R," fits=",nrow(results),
    " success=",sum(results$status=="ok"),
    " failures=",sum(results$status!="ok"),"\n",sep="")
cat("\nPHASE C SELECTED REML (population and estimated I2 in percent)\n")
selected <- subset(sums,method=="reml" &
                    rho_true %in% c(-1,-.99,-.9,0,.9,.99,1) &
                    K %in% c(6,20,80))
selected$pop_I2_ZD <- 100*selected$pop_I2_ZD
selected$I2_median <- 100*selected$I2_median
selected$I2_p975 <- 100*selected$I2_p975
print(selected[,c("scheme","K","rho_true","pop_I2_ZD",
                  "successful","failed","boundary_fraction_success",
                  "discordant_low_I2_fraction_success",
                  "I2_median","I2_p975")],row.names=FALSE,digits=4)
cat("\nCompleted outputs at",outdir,"\n")
