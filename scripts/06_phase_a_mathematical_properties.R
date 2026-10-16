# Phase A: mathematical properties of determinant-based bivariate I^2.
# Base R only; no model fitting or Monte Carlo. Run from repository root.
# These are descriptive *population/matrix* identities, not inferential claims.
#
# References:
# Zhou & Dendukuri (2014), doi:10.1002/sim.6115
# Jackson, White & Riley (2012), doi:10.1002/sim.5453

phase_a_metrics <- function(Tmat, Wmat, rank_tol = 1e-10) {
  stopifnot(is.matrix(Tmat), is.matrix(Wmat),
            identical(dim(Tmat), c(2L, 2L)),
            identical(dim(Wmat), c(2L, 2L)),
            all(is.finite(Tmat)), all(is.finite(Wmat)),
            max(abs(Tmat - t(Tmat))) < 1e-12,
            max(abs(Wmat - t(Wmat))) < 1e-12)
  eigW <- eigen(Wmat, symmetric = TRUE)
  if (any(eigW$values <= 0)) stop("W must be positive definite.")
  eigT <- eigen(Tmat, symmetric = TRUE)$values
  if (any(eigT < -1e-12 * max(1, max(abs(eigT))))) {
    stop("T must be positive semidefinite.")
  }
  W_inv_sqrt <- eigW$vectors %*%
    diag(1 / sqrt(eigW$values), 2) %*% t(eigW$vectors)
  G <- W_inv_sqrt %*% Tmat %*% W_inv_sqrt
  G <- (G + t(G)) / 2
  lambda_raw <- eigen(G, symmetric = TRUE)$values
  if (any(lambda_raw < -1e-10 * max(1, max(abs(lambda_raw))))) {
    stop("Generalized eigenvalues cannot be appreciably negative.")
  }
  lambda <- sort(pmax(lambda_raw, 0))
  # |det(T)|/det(W) = lambda_1 * lambda_2.
  # The generalized eigenvalue representation avoids direct subtraction
  # in det(W+T) and exposes heterogeneity directions.
  root_ratio <- sqrt(prod(lambda))
  zd <- root_ratio / (1 + root_ratio)
  directional <- lambda / (1 + lambda)
  total_volume <- 1 - 1 / sqrt(prod(1 + lambda))
  tr_ratio <- sum(diag(Tmat)) / (sum(diag(Tmat)) + sum(diag(Wmat)))
  eff_rank <- sum(eigT > rank_tol * max(1, max(abs(eigT))))
  c(
    I2_ZD = zd,
    volume_total_exploratory = total_volume,
    trace_ratio_exploratory = tr_ratio,
    directional_min = directional[1],
    directional_max = directional[2],
    generalized_eigen_min = lambda[1],
    generalized_eigen_max = lambda[2],
    marginal_1 = Tmat[1, 1] / (Tmat[1, 1] + Wmat[1, 1]),
    marginal_2 = Tmat[2, 2] / (Tmat[2, 2] + Wmat[2, 2]),
    trace_T = sum(diag(Tmat)),
    determinant_T = prod(pmax(eigT, 0)),
    effective_rank_T = eff_rank
  )
}

phase_a_from_rho <- function(v1, v2, rho, w1, w2) {
  stopifnot(all(is.finite(c(v1, v2, rho, w1, w2))),
            v1 > 0, v2 > 0, w1 > 0, w2 > 0, abs(rho) <= 1)
  Tmat <- matrix(c(v1, rho * sqrt(v1 * v2),
                   rho * sqrt(v1 * v2), v2), 2, byrow = TRUE)
  Wmat <- diag(c(w1, w2))
  met <- phase_a_metrics(Tmat, Wmat)
  # Stable arithmetic for 1-rho^2 near the boundary.
  root_detT <- sqrt(v1 * v2) *
    sqrt((1 - abs(rho)) * (1 + abs(rho)))
  analytic_zd <- root_detT / (root_detT + sqrt(w1 * w2))
  c(rho = rho, rho_abs = abs(rho),
    distance_to_boundary = 1 - abs(rho),
    condition_correlation = if (abs(rho) == 1) Inf else
      (1 + abs(rho)) / (1 - abs(rho)),
    I2_ZD_analytic = analytic_zd, met)
}

outdir <- "phase_A_outputs"
dir.create(outdir, showWarnings = FALSE)

# 1. Boundary: both marginal variances stay at 1 (expected W_ii=0.2).
rho_values <- c(0, 0.50, 0.90, 0.95, 0.99,
                0.999, 0.9999, 0.999999, 1)
grid <- as.data.frame(t(vapply(
  rho_values, function(r) phase_a_from_rho(1, 1, r, .2, .2),
  numeric(length(phase_a_from_rho(1, 1, 0, .2, .2)))
)))
write.csv(grid, file.path(outdir, "phaseA_boundary_grid.csv"),
          row.names = FALSE)

# 2. Rank-one matrices can grow arbitrarily without changing I2_ZD = 0.
scale_values <- c(0.01, 0.1, 1, 10, 100)
rank_one <- do.call(rbind, lapply(scale_values, function(s) {
  Tmat <- s * matrix(c(1, -1, -1, 1), nrow = 2)
  c(scale = s, phase_a_metrics(Tmat, diag(c(.2, .2))))
}))
rank_one <- as.data.frame(rank_one)
write.csv(rank_one, file.path(outdir, "phaseA_rank_one_scale.csv"),
          row.names = FALSE)

# 3. Any simultaneous nonsingular congruence transforms preserve det ratio
# and generalized eigenvalues. Raw trace ratio is NOT affine invariant.
Tbase <- matrix(c(1, .3, .3, 2), nrow = 2)
Wbase <- diag(c(.2, .4))
A <- matrix(c(2, 1, .3, 1), nrow = 2, byrow = TRUE)
stopifnot(abs(det(A)) > 0.01)
Ttrans <- A %*% Tbase %*% t(A)
Wtrans <- A %*% Wbase %*% t(A)
transformation <- rbind(
  c(stage = 0, phase_a_metrics(Tbase, Wbase)),
  c(stage = 1, phase_a_metrics(Ttrans, Wtrans))
)
write.csv(as.data.frame(transformation),
          file.path(outdir, "phaseA_invariance.csv"), row.names = FALSE)

# 4. Empirical plug-in TCCS case from the ARCHIVED fitted model.
# "Reconstructed" row means the 0.5.12 implementation audit has been corrected.
p <- "results/published_reanalyses_all_variants.csv"
if (!file.exists(p)) stop("Expected archived audit file missing: ", p)
raw <- read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
case <- raw[raw$dataset == "pinto_tccs_ischemic" &
            raw$implementation == "full Zhou-Dendukuri reconstruction", ]
stopifnot(nrow(case) == 1)
emp <- phase_a_from_rho(
  case$tau2_Se, case$tau2_FPR, case$rho,
  case$Evar_Se, case$Evar_Sp
)
emp_output <- data.frame(dataset = "pinto_tccs_ischemic",
                         t(emp), archived_I2 = case$I2)
write.csv(emp_output, file.path(outdir, "phaseA_tccs.csv"),
          row.names = FALSE)

# 5. Sensitivity curve with a log boundary-distance parameterization.
log_distance <- seq(0, 14, length.out = 350)
rho_curve <- 1 - 10^(-log_distance)
curve <- as.data.frame(t(vapply(
  rho_curve, function(r) phase_a_from_rho(1, 1, r, .2, .2),
  numeric(ncol(grid))
)))
write.csv(curve, file.path(outdir, "phaseA_curve_data.csv"),
          row.names = FALSE)
draw_figure <- function() {
  matplot(log_distance,
          100 * as.matrix(curve[, c("I2_ZD_analytic",
                                    "volume_total_exploratory",
                                    "trace_ratio_exploratory",
                                    "directional_max")]),
          type = "l", lwd = 2, lty = 1:4, col = 1:4,
          ylim = c(0, 100),
          xlab = expression(-log[10](1 - abs(rho))),
          ylab = "Descriptive matrix index (%)",
          main = "Phase A: determinant collapse at the correlation boundary")
  legend("bottomleft",
         legend = c("Zhou-Dendukuri",
                    "det(W+T) diagnostic (exploratory)",
                    "trace ratio (coordinate-dependent)",
                    "largest generalized direction"),
         col = 1:4, lty = 1:4, lwd = 2, bty = "n", cex = .77)
}
png(file.path(outdir, "phaseA_boundary_curves.png"), width = 2200,
    height = 1450, res = 240)
draw_figure()
dev.off()
pdf(file.path(outdir, "phaseA_boundary_curves.pdf"), width = 9, height = 6)
draw_figure()
dev.off()

capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
cat("Phase A mathematical tables and figures saved to", outdir, "\n")
print(round(grid[, c("rho", "I2_ZD_analytic",
                     "volume_total_exploratory", "marginal_1")], 6),
      row.names = FALSE)
cat("TCCS I2_ZD (archive and eigen diagnostic):",
    case$I2, emp["I2_ZD"], "\n")
