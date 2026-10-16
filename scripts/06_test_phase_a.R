# Phase A regression tests. Does not require mada or any other package.
source("scripts/06_phase_a_mathematical_properties.R")
near <- function(a, b, tol = 1e-8) {
  all(is.finite(c(a, b))) &&
    all(abs(a - b) <= tol * pmax(1, abs(a), abs(b)))
}
W <- diag(c(.2, .2))
m0 <- phase_a_metrics(diag(c(1,1)), W)
plus <- phase_a_metrics(matrix(c(1,1,1,1), nrow=2), W)
minus <- phase_a_metrics(matrix(c(1,-1,-1,1), nrow=2), W)
stopifnot(
  near(m0["I2_ZD"], 5/6),
  near(plus["I2_ZD"], 0),
  near(minus["I2_ZD"], 0),
  near(plus["marginal_1"], 5/6),
  near(plus["marginal_2"], 5/6),
  near(plus["directional_min"], 0),
  near(plus["directional_max"], 10/11),
  near(plus["volume_total_exploratory"], 1 - 1/sqrt(11)),
  plus["effective_rank_T"] == 1,
  minus["effective_rank_T"] == 1
)
# Symmetry for +/- rho and fixed variances and W.
for (r in c(.5, .9, .99, .99999)) {
  a <- phase_a_from_rho(1,1,r,.2,.2)
  b <- phase_a_from_rho(1,1,-r,.2,.2)
  stopifnot(near(a["I2_ZD_analytic"], b["I2_ZD_analytic"]),
            near(a["I2_ZD"], b["I2_ZD"]),
            near(a["trace_T"], b["trace_T"]))
}
# Monotone decline for increasing rho while marginal variances stay fixed.
stopifnot(all(diff(grid$I2_ZD_analytic) <= 1e-10),
          all(abs(grid$marginal_1 - 5/6) < 1e-10),
          all(abs(grid$marginal_2 - 5/6) < 1e-10))

# Rank-one counterexample: between-study variance can become arbitrarily
# large while ZD is identically zero.
stopifnot(
  all(abs(rank_one$I2_ZD) < 1e-8),
  all(diff(rank_one$trace_T) > 0),
  all(rank_one$effective_rank_T == 1),
  rank_one$trace_T[nrow(rank_one)] >= 9000 *
    rank_one$trace_T[1]
)

# Simultaneous invertible affine reparameterization preserves determinant
# ratios and generalized eigenvalues, but NOT the naive trace ratio.
inv <- as.data.frame(transformation)
stopifnot(
  near(inv$I2_ZD[1], inv$I2_ZD[2]),
  near(inv$volume_total_exploratory[1],
       inv$volume_total_exploratory[2]),
  near(inv$generalized_eigen_min[1],
       inv$generalized_eigen_min[2]),
  near(inv$generalized_eigen_max[1],
       inv$generalized_eigen_max[2]),
  abs(inv$trace_ratio_exploratory[1] -
      inv$trace_ratio_exploratory[2]) > .001
)

# Algebraic near-boundary asymptotic rate at fixed T diagonals, W.
for (power in c(5, 7, 9)) {
  delta <- 10^(-power)
  z <- phase_a_from_rho(1,1,1-delta,.2,.2)
  asymptotic <- sqrt(2 * delta) / .2
  stopifnot(abs(z["I2_ZD_analytic"] / asymptotic - 1) < .04)
}
stopifnot(
  near(emp["I2_ZD_analytic"], case$I2, tol=1e-8),
  near(emp["I2_ZD"], case$I2, tol=1e-7),
  abs(case$rho) > .9999999,
  emp["marginal_2"] > .5
)
files <- file.path(outdir, c("phaseA_boundary_grid.csv",
                             "phaseA_rank_one_scale.csv",
                             "phaseA_invariance.csv",
                             "phaseA_tccs.csv",
                             "phaseA_curve_data.csv",
                             "phaseA_boundary_curves.png",
                             "phaseA_boundary_curves.pdf",
                             "sessionInfo.txt"))
stopifnot(all(file.exists(files)), all(file.info(files)$size > 0))
cat("PASS: Phase A algebraic boundary, rank-one, directional, invariance,")
cat(" asymptotic, TCCS archive and generated artifacts.\n")
