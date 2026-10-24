# Stage D - Material scatter and model-form uncertainty

## What was run
- `fe/fatigue.py` extended (previous version: `fe/_backup/fatigue_prestageD.py`): criterion `swt`, a `kf` argument to `g_fun`, `uts_dist` / `ratio_dist` factories and a `CASES` dict. Defaults are unchanged, so the original case ("assumed": UTS CV 3 %, ratio SD 0.05, Kf = Kt = 1.66) is still what `fatigue.py` computes.
- `fe/stage_d.py`: GP surrogate refitted on the 416-run set from Stage A (`results/stageA/fe_samples_al.csv`; Q2 0.9676, CV RMSE 3.20 MPa). 400,000 MC samples, common random numbers across cases (the stress surrogate does not depend on UTS or the ratio); Sobol with 2^13 base samples. Files: `A_literature_cases.csv`, `B_criteria_kf.csv`, `C_robustness.csv`, `C_sobol_all.csv`, `meta.json`.
- Check: the assumed case gives P_f 0.9800 (Goodman) and 0.3365 (Gerber), against 0.9800 / 0.3382 at the end of Stage A (different sample set, MC SE 0.0007).
- Stage C failed its gate, so the rectangular-section FE model is still the basis.

## (1) Published scatter for ISO 5832-2 grade 4 titanium
DOIs below were read from the fetched pages. The numbers come from page text extracted by a fetch tool; I did not open the PDFs or table images.

| Source | What it gives | Use here |
|---|---|---|
| Okazaki, Materials 5(12):2981, 2012, doi:10.3390/ma5122981 | C.P. grade 4 Ti (ISO 5832-2). Annealed UTS 700 +- 8 MPa, 20 % cold-worked 870 +- 2 MPa (n = 3, so CV 1.1 % and 0.2 %). Tension-tension R = 0.1, 10 Hz, 1e7 cycles: fatigue strength / UTS = 0.69 annealed, 0.50 cold-worked. Notched (cold-worked, V-notch Kt 3.2): Kf 1.55, q = 0.25. No fatigue SD reported. | UTS CV 1.1 % ("lit_material"); q = 0.25 Kf option; 0.50-0.69 shows the ratio depends on condition |
| Suzuki et al., Dent Mater J 37(2):308-316, 2018, doi:10.4012/dmj.2017-125 | Grade 4 machined implant specimens, static yield FORCE 1088 +- 93 N (CV 8.5 %); roughened 1118 +- 96 N (8.6 %). The +- is not defined in the abstract, n not seen. | UTS CV 8.5 % ("lit_component"), a component-level proxy only |
| Bologna, Audenino, Terzini, Ann Biomed Eng 52(2):239-249, 2023, doi:10.1007/s10439-023-03363-2 | UTS is taken from material certificates, not tested; the fatigue limit is fixed at 0.5 UTS; no CV or interval given. Two pure Ti plates list 668 and 991 MPa (different plates, not a scatter estimate). Plate 08 (Ti6Al4V) needed about 30 % of the expected load, which the paper attributes to a possible certificate/material mismatch. | Context; confirms certificates carry unquantified error |

Considered and excluded: Wang et al., Materials 11(11):2261, 2018, doi:10.3390/ma11112261 (CP Ti **grade 2**, R = -1, no SD).

**Gap, stated plainly:** I found no published fatigue-strength standard deviation for grade 4 at any R, and no lot-to-lot UTS CV. The endurance-ratio SD (0.05) therefore remains **ASSUMED** in all cases. ASM Handbook, MMPDS and ASTM E739-type datasets were not accessible and not searched in depth. The UTS CVs from n = 3 specimens are very uncertain.

## (A) Literature vs assumed scatter (Kf = Kt, ratio SD 0.05 ASSUMED, no correlation)
| case | UTS CV | P_f Goodman | P_f Gerber | P_f SWT |
|---|---|---|---|---|
| assumed (original) | 3 % | 0.980 | 0.337 | 0.993 |
| lit_material | 1.1 % | 0.987 | 0.326 | 0.995 |
| lit_component | 8.5 % | 0.927 | 0.389 | 0.973 |

## (B) Model form: criterion x notch factor (original scatter)
Nominal stress is always peak / 1.66; Kf only changes the endurance limit.
| criterion | Kf = Kt (1.660) | Peterson (1.641) | Okazaki q = 0.25 (1.165) |
|---|---|---|---|
| Goodman | 0.980 | 0.974 | 0.146 |
| Gerber | 0.337 | 0.302 | 0.0008 |
| SWT | 0.993 | 0.991 | 0.168 |

- SWT, with the same endurance limit and R = 0.1, is the most conservative of the three (nominal g = -75 MPa vs -46 Goodman, +12 Gerber).
- Peterson Kf is only 1.1 % below Kt because the hole radius (3.9 mm) is large against the Peterson constant (a = 0.114 mm). That constant is a **steel** correlation applied to titanium, **ASSUMED**.
- The q = 0.25 option comes from one cold-worked grade 4 V-notch dataset; transfer to a 7.8 mm hole in annealed material is doubtful. It is shown because it is the only measured grade 4 notch sensitivity found, and it moves P_f by 0.5 to 0.8 (an upper bound on how much Kf choice matters), not as a recommended value.

## (C) Robustness grid (`C_robustness.csv`)
P_f ranges over the 18 combinations: Goodman 0.862-0.999, Gerber 0.255-0.411, SWT 0.907-1.000.

Goodman P_f (rows: ratio SD, correlation; columns: UTS CV):
| ratio SD | rho | 2 % | 3 % | 5 % |
|---|---|---|---|---|
| 0.03 | 0 | 0.999 | 0.997 | 0.988 |
| 0.03 | 0.7 | 0.994 | 0.987 | 0.964 |
| 0.05 | 0 | 0.984 | 0.980 | 0.966 |
| 0.05 | 0.7 | 0.968 | 0.954 | 0.924 |
| 0.08 | 0 | 0.929 | 0.925 | 0.913 |
| 0.08 | 0.7 | 0.905 | 0.890 | 0.862 |

Gerber P_f:
| ratio SD | rho | 2 % | 3 % | 5 % |
|---|---|---|---|---|
| 0.03 | 0 | 0.255 | 0.270 | 0.306 |
| 0.03 | 0.7 | 0.285 | 0.306 | 0.340 |
| 0.05 | 0 | 0.330 | 0.337 | 0.354 |
| 0.05 | 0.7 | 0.347 | 0.357 | 0.378 |
| 0.08 | 0 | 0.387 | 0.390 | 0.399 |
| 0.08 | 0.7 | 0.394 | 0.400 | 0.411 |

Sobol ranking (total index of g; full table in `C_sobol_all.csv`):
- With no correlation, the endurance ratio ranks first in 25 of 27 cases (criteria x grid); UTS ranks first in the other two (ratio SD 0.03, UTS CV 5 %, Goodman and Gerber), and UTS is second in 19 of the 25 ratio-first cases; force is second (Goodman and SWT) when UTS CV is 2 %.
- With rho = 0.7, the UTS driver ranks first in all 27 cases and the endurance-ratio residual second. This follows from the construction: the correlated part of the ratio is carried by the UTS driver, so the ranking depends on the correlation assumption. Indices then refer to independent drivers, not the marginal variables.
- Thickness and force are never first; the FE inputs matter little relative to the material terms. Second and third places among small indices are within Sobol sampling error (N = 8192 base, not quantified here) and should not be read as a firm order.

## Conclusions and limitations
- The qualitative picture holds: P_f is high under Goodman and moderate under Gerber, the Goodman-Gerber bracket (0.86-1.0 vs 0.26-0.41) is much wider than any change from the scatter assumptions, and the material terms (endurance ratio, UTS), not the FE inputs, drive the result.
- The largest sensitivities are the choice of criterion and the notch factor, not the distribution widths. Reported intervals therefore understate total uncertainty if only scatter is varied.
- Everything about the endurance-ratio SD, UTS-ratio correlation 0.7 (set by the brief), the Peterson constant and the Kf transfer is ASSUMED or weakly sourced.
- Definition mismatch not resolved: Okazaki's ratio is a tension-tension (R = 0.1) fatigue strength divided by UTS, while the limit state here applies 0.5 UTS as an R = -1 endurance limit following Bologna 2023; I did not verify whether Okazaki's strength is a maximum or amplitude stress.
- The FE stiffness gap from Stage C (model 38-50 % too stiff) is not propagated here.
