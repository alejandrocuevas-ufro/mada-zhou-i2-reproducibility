# Release checklist before manuscript submission

- [x] Run the complete workflow with R and `mada` 0.5.12: `Rscript run_all.R`.
- [x] Confirm every row in `results/published_reproduction_check.csv` is `TRUE`.
- [x] Compare `results/published_reanalyses_summary.csv` with the manuscript table.
- [x] Replace preliminary simulation values in the manuscript with the final validated results.
- [ ] Use `results/Figure1_simulation.pdf` or a journal-compliant export derived from the same data as the submitted figure.
- [x] Inspect `results/sessionInfo.txt` and retain it in the reproducibility record.
- [ ] Retain `results/summary.reitsma_mada_0.5.12.R` as the auditable snapshot of the function examined in the final archived manuscript release.
- [ ] Check study-level transcription against the source articles once more before manuscript submission.
- [ ] Add the final list of manuscript authors and an appropriate code licence.
- [x] Create tagged GitHub release `v1.0.0`.
- [x] Archive release `v1.0.0` in Zenodo.
- [x] Permanent DOI assigned: `10.5281/zenodo.22887677`.

Archived release: https://doi.org/10.5281/zenodo.22887677
