# Release checklist before manuscript submission

- [ ] Run the complete workflow with R and `mada` 0.5.12: `Rscript run_all.R`.
- [ ] Confirm every row in `results/published_reproduction_check.csv` is `TRUE`.
- [ ] Compare `results/published_reanalyses_summary.csv` with the manuscript table.
- [ ] Replace any preliminary simulation values in the manuscript with `results/simulation_summary.csv` from the final R run.
- [ ] Use `results/Figure1_simulation.pdf` or a journal-compliant export derived from the same data as the submitted figure.
- [ ] Inspect `results/sessionInfo.txt` and retain it in the archived release.
- [ ] Retain `results/summary.reitsma_mada_0.5.12.R` as the auditable snapshot of the function examined.
- [ ] Check study-level transcription against the source articles once more before release.
- [ ] Add the final list of manuscript authors and an appropriate code licence.
- [ ] Create a tagged repository release matching the submitted manuscript version.
- [ ] Archive that release in a DOI-granting repository (for example Zenodo).
- [ ] Replace `[DOI DEL REPOSITORIO]` in the Data and Code Availability statement with the permanent DOI.
