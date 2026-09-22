# Run from the repository root:
#   Rscript run_all.R

options(stringsAsFactors = FALSE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

scripts <- c(
  "scripts/00_environment.R",
  "scripts/01_validate_public_examples.R",
  "scripts/02_validate_metadta.R",
  "scripts/03_reanalyse_published.R",
  "scripts/04_simulation.R"
)

for (s in scripts) {
  cat("\n==== Running", s, "====\n")
  source(s, echo = FALSE)
}

cat("\nAll analyses completed. See results/.\n")
