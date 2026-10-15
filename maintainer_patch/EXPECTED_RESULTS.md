# Expected validated results

These values come from the validated `mada` 0.5.12 reproducibility workflow archived at https://doi.org/10.5281/zenodo.22887677.

## Public examples

| Dataset | Native mada 0.5.12 | Psi only | Correct denominators only | Full reconstruction |
|---|---:|---:|---:|---:|
| AuditC6 | 52.523408% | 82.268860% | 49.594700% | 80.493610% |
| AuditC | 40.383519% | 86.046709% | 36.729971% | 84.089019% |

## Published reconstructions

| Analysis | Published | Reproduced native | Full reconstruction |
|---|---:|---:|---:|
| Dragota 2026, SII | 71.8% | 71.766781% | 90.459714% |
| Manetas-Stavrakakis 2023, PPG | 12.5% | 12.528192% | 58.091407% |
| Manetas-Stavrakakis 2023, ECG | 9.2% | 9.227995% | 53.188037% |
| Pinto-Villalba 2026, overall | 14.4% | 14.408619% | 50.568045% |
| Pinto-Villalba 2026, haemorrhagic | 13.1% | 13.118547% | 26.491638% |

## Near-singular example

Pinto-Villalba TCCS ischaemic subgroup:

- `tau2_Se = 0.1026557345`
- `tau2_FPR = 2.4429522576`
- `rho = -0.999999999316995`
- native mada I² = `0.001084535%`
- reconstructed I² = `0.002563587%`

This is a property of the determinant-based statistic and is separate from the implementation discrepancy.

## Monte Carlo experiment

Population I² = `90.557181%`.

| K | Native mada mean | Reconstructed mean |
|---:|---:|---:|
| 6 | 52.254563% | 72.243067% |
| 10 | 50.769600% | 83.302093% |
| 20 | 39.025990% | 88.114837% |
| 40 | 25.973780% | 89.761471% |
| 80 | 15.322911% | 90.142164% |
| 160 | 8.423133% | 90.298877% |
