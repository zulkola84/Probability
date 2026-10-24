# Phase 2, Step 8 - Strength-scatter crossover

## What was run
`fe/strength_crossover.py` (new; `stage_d.physical`, `fatigue`, `sensitivity` imported unchanged). Log: `run.log`. Baseline kept: GP on all 416 FE runs, UTS mean 900 MPa, endurance-ratio mean 0.5, Kf = Kt = 1.66, width fixed, 400,000 common-random-number samples (the Stage D stream, seed 1), Saltelli base 8192, 8 inputs. Check: the cell (ratio SD 0.05, UTS CV 3 %, rho 0) gives Pf 0.9800 / 0.3365 / 0.9932 (Goodman / Gerber / SWT) as in Stage D.

| group | cases | values |
|---|---|---|
| main grid | 35 (rho 0) | ratio SD 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08 (CV 1-16 %, **ASSUMED**) x UTS CV 0.5 % (**ASSUMED**), 1.1 % (Okazaki 2012, doi:10.3390/ma5122981), 3 % (baseline, **ASSUMED**), 5 % (**ASSUMED**), 8.5 % (Suzuki 2018, doi:10.4012/dmj.2017-125, component-level proxy) |
| correlated | 7 | the 7 ratio SDs, UTS CV 3 %, rho 0.7 (**ASSUMED**, from the brief) |
| extension | 20 | ratio SD 0.0025, 0.001, 0.0003, 0.0001 (**ASSUMED**) x the 5 UTS CVs, added because most UTS CVs showed no crossover inside the grid |
| scan | 48 | UTS CV 0.5-5 % x ratio SD 0.0001-0.05 (**ASSUMED**), to locate the UTS CV at which the crossover occurs |

Everything is in `crossover_grid_results.csv` (one row per case and criterion, 330 rows): Pf with MC SE, FORM Pf and beta, total-order Sobol index and its 95 % half-width for all eight inputs, the strength-group index (ST of endurance ratio + ST of UTS), the largest geometric/load index (thickness, force, lever arm, hole diameter) and its name. For the correlated cases the "uts" index is the shared strength driver and "endurance_ratio" the independent residual, as in Table 11. FORM did not converge in 56 of the 330 rows, all with Pf >= 0.999997, where FORM tells nothing.

## Crossover table (rho = 0; `crossover_table.csv`)
Endurance-ratio SD (and CV = SD/0.5) at which the largest geometric/load index first exceeds the index named, scanning down from large SD, linearly interpolated in log SD between grid points ("none" = no crossover down to SD 0.0001):
| criterion | UTS CV | vs strength group | vs endurance ratio alone |
|---|---|---|---|
| Goodman | 0.5 % | 0.0166 (CV 3.3 %) | 0.0170 (3.4 %) |
| | 1.1 % (Okazaki) | 0.0151 (3.0 %) | 0.0170 (3.4 %) |
| | 3 % | none | 0.0169 (3.4 %) |
| | 5 % | none | 0.0168 (3.4 %) |
| | 8.5 % (Suzuki) | none | 0.0167 (3.3 %) |
| Gerber | 0.5 % | 0.0112 (2.2 %) | 0.0116 (2.3 %) |
| | 1.1 % | 0.0094 (1.9 %) | 0.0116 (2.3 %) |
| | 3 % | none | 0.0115 (2.3 %) |
| | 5 % | none | 0.0114 (2.3 %) |
| | 8.5 % | none | 0.0113 (2.3 %) |
| SWT | 0.5 % | 0.0123 (2.5 %) | 0.0126 (2.5 %) |
| | 1.1 % | 0.0114 (2.3 %) | 0.0126 (2.5 %) |
| | 3 % | none | 0.0125 (2.5 %) |
| | 5 % | none | 0.0124 (2.5 %) |
| | 8.5 % | none | 0.0123 (2.5 %) |

- The largest geometric/load index is the **applied force** in every one of the 105 main-grid rows (its scatter, CV 2 %, is ASSUMED), so these crossovers scale with that assumption.
- At UTS CV >= 3 % the strength group never falls below the geometric/load index, even at ratio SD 0.0001: at UTS CV 3 % it is still 0.441 / 0.540 / 0.399 (Goodman / Gerber / SWT) against 0.270 / 0.222 / 0.290. This is the extension; the UTS term alone keeps it up.
- UTS CV at which the strength group equals the geometric/load maximum (`uts_cv_crossover_by_ratio_sd.csv`, linear interpolation in CV): at ratio SD 0.0001, Goodman 2.35 %, Gerber 1.92 %, SWT 2.56 %; at ratio SD 0.005 (CV 1 %) 2.25 / 1.74 / 2.35 %; at ratio SD 0.01 (CV 2 %) 1.90 / 1.01 / 1.59 %. At ratio SD >= 0.02 there is no crossover for any UTS CV up to 5 %.
- Crossover against the endurance-ratio index alone is nearly independent of the UTS CV (spread of 0.0003-0.0004 in SD).
- Sobol indices carry sampling error of up to 0.02 (95 % half-width, main grid); the interpolated SDs are correspondingly uncertain.

## Does "fatigue-strength terms explain most of the limit-state variance" hold?
(Strength-group index >= 0.5; `claim_check.csv`, `cases_where_claim_fails.csv`.)

**Not for every case.** It holds in 102 of the 126 main and correlated rows (criterion x case) and fails in 24: Goodman 10, Gerber 6, SWT 8. It fails only where both scatters are small. In ratio SD (CV) at which the strength group falls below 0.5, interpolated:
| UTS CV | Goodman | Gerber | SWT |
|---|---|---|---|
| 0.5 % | 0.0244 (4.9 %) | 0.0166 (3.3 %) | 0.0182 (3.6 %) |
| 1.1 % | 0.0233 (4.7 %) | 0.0155 (3.1 %) | 0.0175 (3.5 %) |
| 3 % | 0.0110 (2.2 %) | none (0.540 at SD 0.0001) | 0.0105 (2.1 %) |
| 5 %, 8.5 % | none | none | none |

- In the range used in the paper (ratio SD 0.03-0.08, UTS CV 3 %) the strength group is 0.70-0.92 (Goodman), 0.81-0.96 (Gerber) and 0.77-0.95 (SWT), against a geometric/load maximum of 0.02-0.15.
- With correlation 0.7 (UTS CV 3 %) the claim holds in all 21 rows: strength group 0.51-0.97 (the shared driver alone 0.49-0.68), geometric/load maximum 0.015-0.235; no crossover.
- The claim, in the stronger sense (strength group larger than any geometric/load input), fails only for UTS CV <= about 2.5 % together with ratio SD <= about 0.017 (Goodman).

## Pf below the current minimum values (`crossover_grid_results.csv`)
Pf (MC) at UTS CV 3 %, rho 0:
| ratio SD | Goodman | Gerber | SWT |
|---|---|---|---|
| 0.05 | 0.9800 | 0.3365 | 0.9932 |
| 0.03 | 0.9969 | 0.2704 | 0.9998 |
| 0.01 | 0.9998 | 0.1846 | 1.0000 |
| 0.0001 | 0.9999 | 0.1664 | 1.0000 |

At UTS CV 0.5 % and ratio SD 0.0001: 1.0000 / 0.0823 / 1.0000; at UTS CV 8.5 %: 0.9646 / 0.3381 / 0.9970. As the strength scatter shrinks, Goodman and SWT go to 1 (the nominal plate already fails under them) and Gerber falls to the value set by the stress and load scatter alone (0.08-0.17 for UTS CV 0.5-3 %). **The Goodman > Gerber order of Pf, and SWT >= Goodman, hold in all 110 distinct cases** (`pf_order_check.csv`).

## Verdict
- The claim that fatigue-strength terms explain most of the variance of g holds across the range of scatter used in the paper and for every UTS CV >= 5 %; it fails only when the UTS CV is <= 3 % and the ratio SD falls to about 0.01-0.02 (Goodman, SWT) or UTS CV <= 1.1 % and ratio SD <= 0.016 (Gerber).
- Geometric and load inputs (force first) overtake the strength group only if the UTS CV is below about 2-2.5 % *and* the endurance-ratio SD below about 0.009-0.017; against the endurance ratio alone, they overtake at SD about 0.011-0.017 (CV 2.3-3.4 %) irrespective of the UTS CV.
- Lowering the scatter does not remove the Goodman-Gerber bracket: it widens it.

## ASSUMED and limitations
- ASSUMED: every ratio SD, the UTS CVs 0.5, 3 and 5 % (and the scan values), the correlation 0.7, the extension SDs, and the other scatter (force CV 2 %, thickness SD 0.05 mm, hole +-0.1 mm, lever arm +-0.5 mm, E 5 %), with Kf = Kt. The 1.1 % and 8.5 % UTS CVs come from the sources named above (n = 3 specimens for 1.1 %; 8.5 % is a component yield-force CV, a proxy for strength).
- The crossover depends on the force and thickness scatter; if they were larger the crossover would move to larger strength scatter.
- Sobol sampling error up to 0.02; crossovers are linear interpolations; correlated-case indices refer to independent drivers.

## Files
`crossover_grid_results.csv`, `crossover_table.csv`, `uts_cv_crossover_by_ratio_sd.csv`, `claim_check.csv`, `cases_where_claim_fails.csv`, `pf_order_check.csv`, `run.log`, `fig_sobol_vs_ratio_sd.png` (600 dpi, UTS CV 3 %, 0.03-0.08 shaded) and `fig_strength_group_heatmap.png` (600 dpi, main grid, red line where the strength group equals the largest geometric/load index).
