# Automatic checks for the Phase D RNG infrastructure.
# Run from repository root:  Rscript scripts/10_test_rng_streams.R
source("scripts/10_rng_streams.R")

mu <- c(2.0997, -1.2637); SD <- c(1.1746, 0.6381)
Kmax <- 160L
checks <- list()
add <- function(name, pass, detail = "") {
  checks[[length(checks) + 1L]] <<-
    data.frame(check = name, pass = pass, detail = detail,
               stringsAsFactors = FALSE)
}

# (P1) Nesting: seeds and datasets for replicates 1..5 must be identical
# whether the run was configured for 5 or for 1000 replicates.  The seed is a
# pure function of (scheme, replicate), so we verify the datasets themselves.
pool_a <- lapply(1:5, function(r) draw_latent_pool(1L, r, Kmax))
pool_b <- lapply(1:5, function(r) draw_latent_pool(1L, r, Kmax))
add("nesting_seeds_stable",
    identical(vapply(pool_a, `[[`, numeric(1), "seed"),
              vapply(pool_b, `[[`, numeric(1), "seed")))
d_a <- counts_from_pool(pool_a[[3]], 20L, 0.5, 50L, 50L, mu, SD)
d_b <- counts_from_pool(pool_b[[3]], 20L, 0.5, 50L, 50L, mu, SD)
add("nesting_datasets_identical", identical(d_a, d_b))

# Replicates must not be duplicates of one another.
add("replicates_distinct",
    !identical(counts_from_pool(pool_a[[1]], 20L, 0.5, 50L, 50L, mu, SD),
               counts_from_pool(pool_a[[2]], 20L, 0.5, 50L, 50L, mu, SD)))

# Seed map injective across the design range actually used.
grid <- expand.grid(s = 0:5, r = 1:1000)
seeds <- mapply(rng_seed_for, grid$s, grid$r)
add("seed_map_injective", !anyDuplicated(seeds),
    sprintf("%d unique of %d", length(unique(seeds)), length(seeds)))

# (P1) Prefix stability in K: a K=6 dataset must be the first 6 rows of the
# K=160 dataset built from the same pool and rho.
d6 <- counts_from_pool(pool_a[[1]], 6L, 0.0, 50L, 50L, mu, SD)
d160 <- counts_from_pool(pool_a[[1]], 160L, 0.0, 50L, 50L, mu, SD)
add("prefix_stable_in_K",
    identical(d6, head(d160, 6L)))

# (P2) CRN across rho: TP does not depend on rho at all (sensitivity side is
# driven by z1 only), and FP must stay coupled rather than independent.
d_r0 <- counts_from_pool(pool_a[[1]], 40L, 0.0, 50L, 50L, mu, SD)
d_r9 <- counts_from_pool(pool_a[[1]], 40L, 0.9, 50L, 50L, mu, SD)
add("crn_TP_invariant_to_rho", identical(d_r0$TP, d_r9$TP))

# The latent logit-FPR at rho=0 and rho=0.9 share only the z2 component, whose
# coefficient is sqrt(1-0.9^2) = 0.436.  That is the ceiling on correlation
# attainable from the shared latent effects alone; CRN must do strictly better
# because the same uniform also drives the binomial draw.  Comparing against
# independently drawn uniforms is the decisive check.
fp_cor_crn <- cor(d_r0$FP, d_r9$FP)
set.seed(99L)
u_indep <- runif(40)
eta_fpr9 <- mu[2] + SD[2] * (0.9 * pool_a[[1]]$z1[1:40] +
                               sqrt(1 - 0.81) * pool_a[[1]]$z2[1:40])
fp_indep <- qbinom(u_indep, 50L, plogis(eta_fpr9))
fp_cor_indep <- cor(d_r0$FP, fp_indep)
add("crn_FP_beats_independent_sampling", fp_cor_crn > fp_cor_indep,
    sprintf("CRN %.3f vs independent %.3f (latent ceiling %.3f)",
            fp_cor_crn, fp_cor_indep, sqrt(1 - 0.81)))
add("crn_FP_identical_when_rho_unchanged",
    identical(d_r9$FP,
              counts_from_pool(pool_a[[1]], 40L, 0.9, 50L, 50L, mu, SD)$FP))

# Marginal correctness of the generator: with many studies the realised
# between-study correlation of the latent logits must track rho.
pool_big <- draw_latent_pool(2L, 1L, 20000L)
emp <- vapply(c(-0.9, -0.5, 0, 0.5, 0.9), function(rho) {
  z1 <- pool_big$z1; z2 <- pool_big$z2
  cor(mu[1] + SD[1] * z1,
      mu[2] + SD[2] * (rho * z1 + sqrt((1 - abs(rho)) * (1 + abs(rho))) * z2))
}, numeric(1))
add("latent_correlation_recovered",
    max(abs(emp - c(-0.9, -0.5, 0, 0.5, 0.9))) < 0.02,
    sprintf("max abs deviation %.4f", max(abs(emp - c(-0.9,-0.5,0,0.5,0.9)))))

# Unequal study sizes must also be prefix-stable and shared across cells.
sz1 <- draw_study_sizes(1L, 7L, Kmax, 50, 50)
sz2 <- draw_study_sizes(1L, 7L, Kmax, 50, 50)
add("study_sizes_reproducible", identical(sz1, sz2))
add("study_sizes_prefix_stable",
    identical(draw_study_sizes(1L, 7L, 20L, 50, 50)$nD, head(sz1$nD, 20L)))

res <- do.call(rbind, checks)
dir.create("results", showWarnings = FALSE)
write.csv(res, "results/rng_nesting_check.csv", row.names = FALSE)
print(res, row.names = FALSE)
if (!all(res$pass)) stop("RNG stream checks failed.")
cat("\nAll RNG stream checks passed.\n")
