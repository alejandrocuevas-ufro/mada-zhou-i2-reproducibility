# Determinant and |rho| -> 1 in the Zhou-Dendukuri bivariate I^2
# Exploratory analysis v0.1 -- base R, no extra packages.
# Run from the repository root: Rscript scripts/05_determinant_boundary.R
# The theoretical/plug-in diagnostics are not Monte Carlo evidence of model fitting.

out_dir <- "analisis_determinante_rho_v01_salida"
dir.create(out_dir, showWarnings = FALSE)

# All quantities are variances on logit-Se and logit-FPR scales.
# W is the diagonal expected within-study variance from Zhou-Dendukuri.
metrics <- function(v1, v2, rho, w1, w2) {
  stopifnot(all(is.finite(c(v1,v2,rho,w1,w2))),
            v1 >= 0, v2 >= 0, w1 > 0, w2 > 0, abs(rho) <= 1)
  if (v1 == 0 || v2 == 0) rho <- 0  # correlation not identified
  c12 <- rho * sqrt(v1 * v2)
  Tmat <- matrix(c(v1,c12,c12,v2), nrow=2)
  Wmat <- diag(c(w1,w2))
  dT <- v1*v2*(1-rho^2)
  dW <- w1*w2
  zd <- sqrt(dT)/(sqrt(dT)+sqrt(dW))
  tr <- (v1+v2)/(v1+v2+w1+w2)
  dTotal <- (v1+w1)*(v2+w2)-c12^2
  stopifnot(dTotal >= dW*(1-1e-12))
  # Exploratory determinant-of-total index (NOT a published ZD or JWR estimator):
  total_volume <- 1-sqrt(dW/dTotal)
  # Generalized eigenvalues of (T, W), and directional fractions.
  lambda <- sort(eigen(diag(1/sqrt(c(w1,w2))) %*% Tmat %*%
                         diag(1/sqrt(c(w1,w2))),
                       symmetric=TRUE, only.values=TRUE)$values)
  lambda <- pmax(lambda,0)
  directional <- lambda/(1+lambda)
  c(rho=rho, tau2_1=v1, tau2_2=v2,
    within_var_1=w1, within_var_2=w2,
    det_T=dT, det_W=dW,
    I2_ZD=zd, I2_trace_exploratory=tr,
    I2_total_volume_exploratory=total_volume,
    I2_direction_small=directional[1],
    I2_direction_large=directional[2],
    I2_marginal_1=v1/(v1+w1), I2_marginal_2=v2/(v2+w2))
}

rho_values <- c(0,.50,.90,.95,.99,.999,.9999,.999999,1)
controlled <- as.data.frame(t(vapply(rho_values, function(r)
  metrics(v1=1,v2=1,rho=r,w1=.2,w2=.2), numeric(14))))
write.csv(controlled, file.path(out_dir,"escenario_controlado.csv"), row.names=FALSE)

# Real application: Pinto-Villalba, TCCS ischaemic subgroup, 8 arms.
# Fallback values are from the archived reconstruction.
real <- list(v1=0.102655734536248, v2=2.44295225756851,
             rho=-0.999999999316995,
             w1=0.305214335852786, w2=1.70776819892622)
input <- "results/published_reanalyses_all_variants.csv"
if (file.exists(input)) {
  d <- read.csv(input, check.names=FALSE)
  sel <- d[d$dataset=="pinto_tccs_ischemic" &
           d$implementation=="full Zhou-Dendukuri reconstruction", ]
  if (nrow(sel)==1) {
    real <- list(v1=sel$tau2_Se, v2=sel$tau2_FPR,
                 rho=sel$rho, w1=sel$Evar_Se, w2=sel$Evar_Sp)
  }
}
real_result <- as.data.frame(as.list(metrics(real$v1,real$v2,real$rho,
                                          real$w1,real$w2)))
real_result$dataset <- "pinto_tccs_ischemic"
write.csv(real_result, file.path(out_dir,"caso_real_tccs.csv"), row.names=FALSE)

# Graph: horizontal axis is -log10(1-|rho|); equal theoretical marginals.
x <- seq(0,10,length.out=401)
r <- 1-10^(-x)
r[1] <- 0
curves <- as.data.frame(t(vapply(r, function(z)
  metrics(1,1,z,.2,.2), numeric(14))))
png(file.path(out_dir,"sensibilidad_rho.png"),width=1800,height=1100,res=210)
matplot(x, 100*as.matrix(curves[,c("I2_ZD", "I2_trace_exploratory",
                                     "I2_total_volume_exploratory", "I2_direction_large")]),
        type="l",lwd=2,lty=1:4,
        xlab=expression(-log[10](1-abs(rho))), ylab="Index value (%)",
        ylim=c(0,100),main="Determinant collapse as |rho| approaches 1")
legend("bottomleft",lty=1:4,lwd=2,bty="n",
       legend=c("Zhou-Dendukuri", "Trace ratio (exploratory)",
                "Total determinant (exploratory)", "Strongest direction (diagnostic)"))
dev.off()

cat("\nControlled scenario (percentage):\n")
report <- controlled[,c("rho","I2_ZD","I2_trace_exploratory",
                         "I2_total_volume_exploratory","I2_direction_large")]
report[,-1] <- report[,-1]*100
print(round(report,5), row.names=FALSE)
cat("\nReal case metrics (fractions unless labeled):\n")
print(real_result, digits=10, row.names=FALSE)
cat("\nFiles saved at:",normalizePath(out_dir),"\n")
