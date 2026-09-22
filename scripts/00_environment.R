# Record the computational environment and the exact summary.reitsma source.

if (!requireNamespace("mada", quietly = TRUE)) {
  stop("Package 'mada' is required. See README.md.")
}

required_mada <- package_version("0.5.12")
installed_mada <- packageVersion("mada")
if (installed_mada != required_mada) {
  stop(sprintf(
    "This reproduction targets mada 0.5.12; installed version is %s. See README.md for installation instructions.",
    installed_mada
  ))
}

library(mada)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

sink("results/sessionInfo.txt")
cat("Reproducibility record\n")
cat("======================\n")
cat("mada version:", as.character(packageVersion("mada")), "\n")
if (requireNamespace("MASS", quietly = TRUE)) {
  cat("MASS version:", as.character(packageVersion("MASS")), "\n")
}
cat("\n")
print(sessionInfo())
sink()

f <- getS3method("summary", "reitsma")
writeLines(
  c(
    sprintf("# Captured from mada %s", packageVersion("mada")),
    sprintf("# R %s", getRversion()),
    deparse(f, width.cutoff = 120)
  ),
  "results/summary.reitsma_mada_0.5.12.R"
)

cat("Environment and summary.reitsma source recorded.\n")
