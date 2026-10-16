# Dependency-free verification of the determinant/rho boundary diagnostic.
# Run from the repository root: Rscript scripts/05_test_determinant_boundary.R
source("scripts/05_determinant_boundary.R")

near <- function(a, b, tol = 1e-9) {
  isTRUE(all.equal(unname(a), unname(b), tolerance = tol))
}

# Strict boundary: positive marginal heterogeneity, zero determinant.
z0 <- metrics(1, 1, 0, .2, .2)
zp <- metrics(1, 1, 1, .2, .2)
zn <- metrics(1, 1, -1, .2, .2)
stopifnot(near(z0["I2_ZD"], 5/6),
          near(zp["I2_ZD"], 0), near(zn["I2_ZD"], 0),
          near(zp["I2_marginal_1"], 5/6),
          near(zp["I2_marginal_2"], 5/6),
          near(zp["I2_trace_exploratory"], 5/6),
          zp["I2_total_volume_exploratory"] > 0,
          zp["I2_direction_large"] > 0,
          near(zp["I2_direction_small"], 0))

# Positive/negative rho symmetry and monotone approach to singularity.
for (r in c(0.5, 0.9, 0.99, 0.999999)) {
  a <- metrics(1, 1, r, .2, .2)
  b <- metrics(1, 1, -r, .2, .2)
  stopifnot(near(a["I2_ZD"], b["I2_ZD"]),
            near(a["det_T"], b["det_T"]))
}
stopifnot(metrics(1, 1, .9, .2, .2)["I2_ZD"] >
            metrics(1, 1, .99, .2, .2)["I2_ZD"],
          metrics(1, 1, .99, .2, .2)["I2_ZD"] >
            metrics(1, 1, .999999, .2, .2)["I2_ZD"])

# Compare the diagnostic to the previously archived empirical reconstruction.
source_path <- "results/published_reanalyses_all_variants.csv"
stopifnot(file.exists(source_path))
d <- read.csv(source_path, check.names = FALSE)
selected <- subset(d, dataset == "pinto_tccs_ischemic" &
                     implementation == "full Zhou-Dendukuri reconstruction")
stopifnot(nrow(selected) == 1L,
          near(real_result$I2_ZD, selected$I2, tol = 1e-7),
          near(real_result$det_T, selected$det_between, tol = 1e-7),
          nrow(controlled) == 9L)

paths <- file.path("analisis_determinante_rho_v01_salida",
                   c("escenario_controlado.csv", "caso_real_tccs.csv",
                     "sensibilidad_rho.png"))
stopifnot(all(file.exists(paths)), all(file.info(paths)$size > 0))
cat("PASS: determinant boundary, rho symmetry, archived TCCS case and output files.\n")
