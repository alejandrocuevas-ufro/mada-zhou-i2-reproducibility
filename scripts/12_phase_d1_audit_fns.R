# Shared functions for the Phase D1 audit (sourced by the audit script and by
# the validation tests).  Kept in one place so the objective that produces the
# Hessian diagnostics is provably the one mada maximises.

# --- mada's two-stage transform, reproduced exactly -------------------------
# correction.control = "all" adds the correction to EVERY table, but only when
# at least one zero cell exists anywhere in the dataset; "single" adds it only
# to the tables that contain a zero.  Both branches are reproduced so the
# objective below is the same likelihood mada maximises.
mada_transform <- function(counts, correction = 0.5,
                           correction.control = c("all", "single")) {
  correction.control <- match.arg(correction.control)
  m <- as.matrix(counts[, c("TP", "FN", "FP", "TN")])
  if (correction.control == "all") {
    if (any(m == 0)) m <- m + correction
  } else {
    add <- apply(m == 0, 1, any) * correction
    m <- m + add
  }
  Y <- cbind(log(m[, "TP"] / m[, "FN"]), log(m[, "FP"] / m[, "TN"]))
  S <- cbind(1/m[, "TP"] + 1/m[, "FN"], 1/m[, "FP"] + 1/m[, "TN"])
  list(Y = Y, S = S)
}

# -2 log (restricted) likelihood up to an additive constant, parameterised by
# the unconstrained Cholesky factor L = (L11, L21, L22).
neg2ll_chol <- function(par, Y, S, reml = TRUE) {
  L11 <- par[1]; L21 <- par[2]; L22 <- par[3]
  v1 <- L11^2; c12 <- L11 * L21; v2 <- L21^2 + L22^2
  a <- S[, 1] + v1; b <- S[, 2] + v2
  d <- a * b - c12^2
  if (any(!is.finite(d)) || any(d <= 0)) return(1e100)
  i11 <- b / d; i12 <- -c12 / d; i22 <- a / d
  A <- matrix(c(sum(i11), sum(i12), sum(i12), sum(i22)), 2, 2)
  rhs <- c(sum(i11 * Y[, 1] + i12 * Y[, 2]), sum(i12 * Y[, 1] + i22 * Y[, 2]))
  beta <- tryCatch(solve(A, rhs), error = function(e) NA_real_)
  if (any(!is.finite(beta))) return(1e100)
  e1 <- Y[, 1] - beta[1]; e2 <- Y[, 2] - beta[2]
  obj <- sum(log(d)) + sum(i11 * e1^2 + 2 * i12 * e1 * e2 + i22 * e2^2)
  if (reml) {
    ld <- determinant(A, logarithm = TRUE)
    if (ld$sign <= 0) return(1e100)
    obj <- obj + as.numeric(ld$modulus)
  }
  obj
}

# Profile over rho: reoptimise the two marginal SDs at each fixed rho.
# Parameterised as (log s1, log s2) with L11 = s1, L21 = rho*s2, L22 =
# sqrt(1-rho^2)*s2 so that rho is held exactly.
profile_over_rho <- function(Y, S, reml = TRUE,
                             grid = c(-0.99, -0.9, -0.5, 0, 0.5, 0.9, 0.99),
                             start = c(0, 0)) {
  vapply(grid, function(rho) {
    f <- function(lp) {
      s1 <- exp(lp[1]); s2 <- exp(lp[2])
      neg2ll_chol(c(s1, rho * s2, sqrt(max(0, 1 - rho^2)) * s2), Y, S, reml)
    }
    o <- tryCatch(stats::optim(start, f, method = "L-BFGS-B",
                               lower = c(-12, -12), upper = c(5, 5),
                               control = list(maxit = 200)),
                  error = function(e) NULL)
    if (is.null(o)) NA_real_ else o$value
  }, numeric(1))
}


# --- The singular solution is always a stationary point ---------------------
# Structural fact about mvmeta's parameterisation: the objective depends on
# L22 only through v2 = L21^2 + L22^2, so it is an EVEN function of L22 and
# its first derivative in L22 vanishes identically at L22 = 0.  |rho| = 1 is
# therefore a critical point of the likelihood for EVERY dataset, not an
# artefact of particular data.  Whether the optimiser is attracted to it
# depends on the sign of the second derivative there: positive curvature of
# -2logL makes the singular solution a local MINIMUM, i.e. attracting.
singular_point_diagnostics <- function(Y, S, reml = TRUE,
                                       unconstrained_value = NA_real_) {
  f0 <- function(p) neg2ll_chol(c(exp(p[1]), p[2], 0), Y, S, reml)
  best <- NULL
  for (st in list(c(0, 0), c(0, 0.5), c(0, -0.5), c(-1, 0.1), c(1, -0.1))) {
    o <- tryCatch(stats::optim(st, f0, method = "Nelder-Mead",
                               control = list(maxit = 2000, reltol = 1e-11)),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && (is.null(best) || o$value < best$value))
      best <- o
  }
  if (is.null(best))
    return(list(curvature = NA_real_, is_local_min = NA,
                value = NA_real_, excess = NA_real_))
  L11s <- exp(best$par[1]); L21s <- best$par[2]
  h <- 1e-4
  f <- function(x) neg2ll_chol(c(L11s, L21s, x), Y, S, reml)
  curv <- (f(h) - 2 * f(0) + f(-h)) / h^2
  list(curvature = curv,
       is_local_min = isTRUE(curv > 0),
       value = best$value,
       excess = if (is.finite(unconstrained_value))
         best$value - unconstrained_value else NA_real_)
}
