# Phase 2, Step 7 - Plate width as epistemic uncertainty

## 1. Earlier width run and its reuse
Search of `results/` and `fe/` found it:
- **`results/old_with_width/`**: `fe_samples.csv` (256 FE runs, columns width, thickness, hole_d, E, nu, lever_arm, peak_vm, stiffness), `sobol_indices.csv`, `sensitivity_summary.csv`, `sensitivity.png`. Width uniform 16.365 +- 1.0 mm (15.369-17.359 sampled); the other five inputs from the same distributions as now; scrambled Sobol (seed 1); hole-edge size 0.55 mm, global 3.0 mm. Its Sobol result is the "about 82 %" you remember: width S1 = 0.817, ST = 0.818 for the peak stress.
- The code of that run is saved in `fe/_backup_before_width_fix/` (`uncertainty.py`, `run_samples.py`); `results/stage3_report.md` records why width was later held fixed.

**Reuse check.** I re-ran 8 of the 256 samples (including the smallest and largest width) with the current code (`verify_old_runs.csv`): peak differs by up to 3.23 MPa (0.37 % of peak; RMS 1.27 MPa), stiffness by up to 0.017 N/mm. The peak differences are of the size of the FE scatter between meshes seen before (GP noise about 2-3 MPa), so I reused all 256 runs and ran no new design. New FE runs in this step: **13** (8 verification, 5 nominal-width runs), about 1.1-1.3 s of solver time each (9.7 s and 5.9 s in total). The solver time of the 256 reused runs was not recorded.

## 2. Six-input GP (`cv_scores.csv`, `cv_predictions.csv`)
Training always included the 416 fixed-width runs (width column = 16.365 mm); 5-fold CV held out only width-varying points (256 points in total).
| | value |
|---|---|
| Q2 | 0.9989 |
| RMSE | 1.79 MPa (stage F, fixed width: 2.89) |
| bias | +0.02 MPa |
| max error | 6.65 MPa (0.80 % of peak) |
| RMSE for \|w - 16.365\| > 0.75 mm / inner points | 1.10 / 1.97 MPa |

Acceptance (ASSUMED): Q2 >= 0.95 and RMSE <= 1 % of the mean peak (8.6 MPa). The fit passes with margin, so I went on. Kernel and length scales (final GP on all 672 runs): `gp6_length_scales.csv`, `meta.json`.

## 3-4. Cases, Pf, FORM, runout, Sobol
W0: width fixed; W1: uniform 16.365 +- 0.5 mm and W2: +- 1.0 mm (both **ASSUMED** epistemic ranges, treated here as uniform scatter). Stage D scatter for everything else, Kf = Kt, 400,000 CRN samples (seed 1; the first 8 columns are the Stage D stream, the width has its own common stream).

**W0 reproduces the baseline:** Pf 0.9800 / 0.3364 / 0.9932 (required 0.9800 / 0.3365 / 0.9932; the difference 0.0001 is inside the MC SE of 0.0007). The Stage D 5-input GP gives 0.9800 / 0.3365 / 0.9932 (`W0_vs_5input_GP.csv`; RMS stress difference between the two GPs 0.22 MPa).

| case | criterion | Pf (MC SE) | FORM Pf (beta) | runout load, nominal (N, vs 1200) | load at Pf = 0.5 (N) |
|---|---|---|---|---|---|
| W0 | Goodman | 0.9800 (0.0002) | 0.9783 (-2.019) | 1024.3 (-14.6 %) | 1022.3 |
| | Gerber | 0.3364 (0.0007) | 0.3312 (0.437) | 1250.3 (+4.2 %) | 1247.8 |
| | SWT | 0.9932 (0.0001) | 0.9928 (-2.448) | 940.1 (-21.7 %) | 938.8 |
| W1 | Goodman | 0.9721 (0.0003) | 0.9646 (-1.807) | same | 1021.5 |
| | Gerber | 0.3471 (0.0008) | 0.3446 (0.400) | same | 1246.7 |
| | SWT | 0.9909 (0.0002) | 0.9886 (-2.276) | same | 937.9 |
| W2 | Goodman | 0.9487 (0.0003) | 0.9313 (-1.486) | same | 1019.4 |
| | Gerber | 0.3745 (0.0008) | 0.3719 (0.327) | same | 1243.9 |
| | SWT | 0.9834 (0.0002) | 0.9769 (-1.994) | same | 936.0 |

The nominal runout load is identical in all cases because the width distributions are symmetric about 16.365 mm (the nominal plate is the same); the "load at Pf = 0.5" column shows the shift. FORM converged in all nine cases (negative beta: the nominal plate already fails). Pf change from W0 (paired SE <= 0.0007): W1 Goodman -0.0079, Gerber +0.0107, SWT -0.0023; W2 Goodman -0.0313, Gerber +0.0382, SWT -0.0098 (`pf_differences_vs_W0.csv`).

Peak-stress sd at 1200 N: 24.7 (W0), 35.7 (W1), 57.6 MPa (W2).

**Sobol indices** (total order, base 8192, 9 variables; full tables `sobol_by_case.csv`):
| | W0 | W1 | W2 |
|---|---|---|---|
| width, peak stress | 0 | 0.522 | 0.817 |
| force, peak stress | 0.484 | 0.232 | 0.090 |
| thickness, peak stress | 0.415 | 0.198 | 0.075 |
| width, g Goodman / Gerber / SWT | 0 | 0.156 / 0.091 / 0.108 | 0.431 / 0.291 / 0.330 |
| endurance ratio, g Goodman | 0.698 | 0.589 | 0.397 |
| endurance ratio, g Gerber | 0.802 | 0.729 | 0.569 |
| endurance ratio, g SWT | (see CSV) | (see CSV) | 0.548 |

The W2 peak-stress index (0.8175) equals the earlier run's 0.818, as it should (same range).

## 5. Nominal peak vs width (other inputs nominal, 1200 N; `nominal_peak_vs_width.csv`)
| width (mm) | GP (MPa) | FE (MPa) | FE vs published 857.8 |
|---|---|---|---|
| 15.4 | 953.8 | 955.8 | +11.4 % |
| 16.0 | 890.3 | 890.1 | +3.8 % |
| 16.365 | 856.3 | 858.5 | +0.1 % |
| 16.8 | 819.4 | 819.7 | -4.4 % |
| 17.4 | 773.6 | 773.6 | -9.8 % |

The GP follows the FE runs to within 2.2 MPa. 17.4 mm lies 0.04 mm beyond the largest width in the training data; the GP and FE agree there (0.04 MPa). The slope near the calibrated width is -90.8 MPa/mm (GP): the peak moves 0.5 % away from 857.8 MPa per 0.05 mm, +-45 MPa (+-5 %) over +-0.5 mm and +-91 MPa (+-11 %) over +-1.0 mm.

## 6. Do the conclusions survive?
- **(a) Ordering of Pf: yes.** SWT > Goodman > Gerber in W0, W1 and W2. The Goodman-Gerber gap narrows from 0.644 to 0.625 (W1) and 0.574 (W2).
- **(b) Bracketing of 1200 N: yes.** Nominal runout loads 1024 N (Goodman) and 1250 N (Gerber) are unchanged by construction; the loads at Pf = 0.5 (1019-1022 N and 1244-1248 N) stay on either side of 1200 N in every case.
- **(c) Dominance of the endurance ratio: partly.** It stays first in W1 for all three criteria (width second for Goodman and SWT, third for Gerber). In W2 it stays first for Gerber (0.569 against width 0.291) and SWT (0.548 against 0.330), **but not for Goodman**, where width ranks first (0.431 against 0.397).
- **Width shares.** Of the peak-stress variance, width explains 52.2 % in W1 and 81.7 % in W2. Of the variance of g: W1 15.6 % (Goodman), 9.1 % (Gerber), 10.8 % (SWT); W2 43.1 %, 29.0 %, 33.0 %.

**Verdict.** The qualitative conclusions (criterion ordering and the bracketing of the 1200 N runout) survive a width error of +-0.5 and +-1.0 mm. The quantitative Pf values move by 0.01-0.04, much less than the Goodman-Gerber gap. The dominance of the endurance ratio survives +-0.5 mm and, for Gerber and SWT, +-1.0 mm. For Goodman at +-1.0 mm (+-6 % of the width) the width becomes the largest source of variance in g; a reviewer asking "what if the calibrated width is wrong" gets: modest changes in Pf, same conclusions, except that the sensitivity ranking under Goodman then depends on the width assumption.

## ASSUMED and limitations
- ASSUMED: the W1 (+-0.5 mm) and W2 (+-1.0 mm) ranges, and the choice of a uniform distribution; the Stage D scatter (UTS CV 3 %, ratio SD 0.05, force CV 2 %, FE input ranges), Kf = Kt; the acceptance thresholds (CV: Q2 >= 0.95, RMSE <= 1 % of peak; reuse: 1 % of peak).
- Width is epistemic: a wrong width shifts every sample in the same direction, so a uniform mixture averages over it. A Pf conditional on one wrong width (for example 15.4 mm, +11 % in peak stress) was not computed; the mixture Pf values here hide that dependence.
- The width runs come from an older code state (differences up to 3.2 MPa, 0.37 %) and cover +-1.0 mm only; solver times of those runs are unknown.
- The width acts here only through the net section and notch factor of a rectangular section; the real section is curved (Stage C) and the stiffness gap of Stage C is not addressed.
- The 6-input GP is not an extrapolation model outside 15.37-17.36 mm.

## Files
`verify_old_runs.csv`, `cv_scores.csv`, `cv_predictions.csv`, `gp6_length_scales.csv`, `nominal_peak_vs_width.csv`, `pf_form_runout_by_case.csv`, `pf_differences_vs_W0.csv`, `W0_vs_5input_GP.csv`, `sobol_by_case.csv`, `sobol_top3_g.csv`, `meta.json`, `run.log`, `verify.log`, `fe_decks/` (CalculiX decks and results of the 13 new runs), `fig_peak_vs_width.png` and `fig_pf_by_width_case.png` (600 dpi).
