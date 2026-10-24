# Phase 2, Step 3 - Closed-form baseline vs Gaussian process

## What was run
`fe/closed_form_baseline.py` (new script; nothing existing edited, `fatigue.py` and `stage_f.py` imported unchanged). Run log: `run.log`. All outputs are in this folder; nothing was written to `results/`.

- Closed form: sigma_pk = Kt_net x sigma_nom, sigma_nom = 6M / ((W - d) t^2), M = F h / 2, W = 16.365 mm (calibrated, never sampled). E and nu do not enter.
- **Step 1 check passed:** at nominal inputs (t 5.12, d 7.8, h 38, F 1200) sigma_nom = **609.28 MPa**, equal to `nominal_rect_net_MPa` in `results/stress_field/summary.json`.
- Kt versions:
  - **(a)** fixed at 1.4116 (`Kt_rect_net`, `results/stress_field/summary.json`); gives 860.1 MPa at nominal, against the FE value 859.1.
  - **(b)** Howland/Heywood finite-width fit, Kt_net = 3.00 - 3.13(d/W) + 3.66(d/W)^2 - 1.53(d/W)^3 (Peterson's Stress Concentration Factors, Chart 4.1). **ASSUMED:** the coefficients are quoted from memory and I could not confirm them. Two searches did not find a bending formula for this geometry (hole axis through the thickness, bending stress varying along the hole axis). This formula is for **in-plane tension**, which is a different loading case. At d/W = 0.477 it gives 2.174.
- Hold-out: the same 84-point split as Stage F (`train_test_split`, seed 0, stratified DOE/AL). The GP for this part is trained on the other 332 runs and reproduces Stage F (RMSE 2.890, Q2 0.9735). The closed forms are not fitted to the hold-out data.
- Monte Carlo: Stage D "assumed" scatter, Kf = Kt, 400,000 common-random-number samples (seed 1), so the GP column reproduces Stage D (0.9800 / 0.3365 / 0.9932). Here the GP is the one trained on all 416 runs. FORM: HL-RF as in `fatigue.py`. Sobol: total order of g, base sample 8192, 8 variables.

## Peak stress on the 84 hold-out points (`holdout_scores.csv`)
| model | RMSE (MPa) | bias (MPa) | max error (MPa / % of peak) | Q2 |
|---|---|---|---|---|
| GP | 2.89 | -0.38 | 7.22 / 0.86 % | 0.9735 |
| closed form (a), Kt = 1.4116 | 4.65 | +2.73 | 12.49 / 1.46 % | 0.9315 |
| closed form (b), Howland | 467.5 | +467.4 | 497.4 / 56.1 % | -692 |

By subset, closed form (a) has RMSE 4.22 (DOE points) and 5.27 (active-learning points).

Where closed form (a) falls short:
- It is biased high by 2.7 MPa, which is 0.3 % of the peak. I did not test the cause. One possibility is the mesh of the single run that Kt comes from (0.3 mm hole size, against 0.55 mm for the DOE runs).
- Its residual correlates with hole diameter (Pearson r = +0.62; `holdout_residual_correlation.csv`). A fixed Kt assumes the notch factor does not change with d. The local elasticities of peak stress with respect to d at nominal inputs (`nominal_elasticities.csv`) are 0.91 (closed form) against 0.60 (GP). The elasticities for thickness (-2.00 against -1.91) and lever arm (1.00 against 1.08) are close.
- The GP shows that E and nu have no effect on the stress (elasticities 0.001 and 0.017), consistent with ignoring them.

Closed form (b) fails completely: it overpredicts by 56 % because it uses the wrong loading case, as described above. Its Q2 is negative.

## P_f, FORM (400,000 samples; `pf_by_method.csv`, `pf_differences.csv`, `form_by_method.csv`)
| method | Goodman P_f | Gerber P_f | SWT P_f |
|---|---|---|---|
| GP (416 runs) | 0.9800 (SE 0.0002) | 0.3365 (0.0007) | 0.9932 (0.0001) |
| closed form (a) | 0.9820 (0.0002) | 0.3521 (0.0008) | 0.9939 (0.0001) |
| closed form (b) | 1.0000 | 1.0000 | 1.0000 |

Difference from the GP and its size relative to the MC error:
| | Goodman | Gerber | SWT |
|---|---|---|---|
| (a) minus GP | +0.0019 | +0.0156 | +0.0006 |
| SE of the difference, independent / paired | 0.0003 / 0.0001 | 0.0011 / 0.0002 | 0.0002 / <0.0001 |
| \|diff\| / independent SE | 6.3 | 14.7 | 3.5 |
| within 1.96 SE? | no | no | no |
| samples classified differently | 897 | 6,659 | 286 |

FORM (beta is negative where the nominal plate already fails, the convention of `fatigue.py`):
| | Goodman P_f (beta) | Gerber P_f (beta) | SWT P_f (beta) |
|---|---|---|---|
| GP | 0.9782 (-2.017) | 0.3309 (0.437) | 0.9928 (-2.446) |
| (a) | 0.9806 (-2.066) | 0.3490 (0.388) | 0.9935 (-2.486) |
| (b) | 1.0000 (-6.44, not converged) | 1.0000 (-5.20) | 1.0000 (-4.95, not converged) |

Diagnostic, not a model (`diagnostic_cf_a_offset_removed.csv`): subtracting closed form (a)'s mean offset to the GP over the MC pool (3.39 MPa) gives 0.9794 / 0.3367 / 0.9931. That is within MC error for Gerber and SWT, and 0.0006 low for Goodman (about 2 SE, flagged "not within"). So most of the P_f difference is the +3.4 MPa mean offset, not a difference in how stress depends on the inputs. The offset was taken from the GP, so this shows the cause of the difference, not an independent result.

## Total-order Sobol indices of g (`sobol_by_method.csv`, `sobol_differences.csv`, `sobol_top3.csv`)
Closed form (a) against the GP, Goodman (Gerber and SWT are similar; the differences are in the CSV):
| variable | GP | (a) | difference | combined 95 % half-width |
|---|---|---|---|---|
| endurance ratio | 0.696 | 0.687 | -0.010 | 0.025 |
| UTS | 0.135 | 0.134 | -0.001 | 0.007 |
| force | 0.082 | 0.081 | 0.000 | 0.004 |
| thickness | 0.070 | 0.078 | +0.008 | 0.004 |
| lever arm | 0.014 | 0.012 | -0.002 | 0.001 |
| hole diameter | 0.004 | 0.009 | +0.005 | 0.0004 |
| E, nu | 0.000 | 0.000 | 0 | 0 |

The ranking is identical (endurance ratio > UTS > force) for all three criteria. The large indices agree within the Sobol error. Thickness and hole diameter are slightly too high in the closed form, beyond the sampling error but small (at most 0.008 for thickness and 0.005 for hole diameter, in the Goodman and Gerber/SWT rows). Closed form (b) changes the picture: for Goodman and SWT, force and thickness replace UTS in the top three, and the endurance-ratio index falls (0.41 against 0.70 for Goodman), probably because its stress is far too high (an artefact of the wrong Kt; I did not examine the mechanism).

## Verdict
1. **Closed form (a) and the GP do not agree within the Monte Carlo error.** Goodman, Gerber and SWT differ by 6, 15 and 3.5 combined standard errors. In practice the differences are small: +0.002, +0.016 and +0.001 in P_f, against a Goodman-Gerber spread of about 0.64 from the criterion choice. They have the same FORM beta to within 0.05, and the same Sobol ranking.
2. The disagreement is mostly a 3.4 MPa (0.4 %) mean offset (see the diagnostic). The closed form has no hole-diameter dependence of Kt, which matters little here because d varies by only 1.3 %.
3. **Kt version:** (a) is far closer. Version (b) is wrong by 56 % (P_f = 1, negative Q2), because the formula is for in-plane tension. This says nothing about a correct bending formula, which I could not find.
4. **Is the GP needed?** For P_f and the ranking at this accuracy, the one-run closed form (a) reproduces the GP to about 1.6 points of P_f at worst, but it is biased and its error (4.65 MPa RMSE) is 60 % above the GP's (2.89). If the target is agreement within MC error, the GP is needed (or at least a closed form with a bias/d-dependence correction).
5. **Limitations:** (a) takes its Kt from the same FE model, so it is self-consistent with the FE data and is not an independent check of the physics. The hold-out set has 84 points, 32 of them active-learning points near the failure boundary. Sobol indices have sampling error (shown in the CSV). All scatter assumptions (UTS CV 3 %, ratio SD 0.05, force CV 2 %, FE input ranges) are ASSUMED, as in Stage D. The stage C stiffness gap is not propagated.

## Files
`holdout_scores.csv`, `holdout_predictions.csv`, `holdout_residual_correlation.csv`, `nominal_elasticities.csv`, `pf_by_method.csv`, `pf_differences.csv`, `diagnostic_cf_a_offset_removed.csv`, `form_by_method.csv`, `sobol_by_method.csv`, `sobol_differences.csv`, `sobol_top3.csv`, `meta.json`, `run.log`, `comparison.png` (600 dpi; left: predicted vs FE for the 84 points, right: P_f by method).
