# Stage A - Active learning of the peak-stress surrogate

## What was run
- `fe/active_learning.py` (previous version backed up to `fe/_backup/active_learning_prestageA.py`).
- Start: `results/fe_samples_initial.csv` (original 256-run Sobol DOE; identical to the first 256 rows of the enlarged `results/fe_samples.csv` left by an earlier, non-conforming run).
- Per iteration, per criterion (Goodman, Gerber): 150,000-point MC pool of all 8 random inputs, learning function U = |g| / sigma_g (sigma_g from the GP predictive sd propagated through the limit state), 8 lowest-U points chosen greedily with a minimum separation of 0.35 input-sd (a space-filling choice, so not the 8 strictly lowest U), FE run via `run_samples.run_one`, surrogate-vs-FE error recorded before refit, refit, P_f recomputed (200,000 MC samples, fixed seed).
- 10 iterations x 16 FE runs = 160 new runs; final training set 416 runs.
- Stop rule (brief): min U >= 2 for both criteria, or 10 iterations. **min U never reached 2** (best value 0.003), so the run ended at the 10-iteration cap.
- No existing result file was modified. Outputs are only in this folder: `al_log.csv`, `al_new_points.csv`, `fe_samples_al.csv`, `pf_vs_fe_runs.png`.

## Key numbers (from `al_log.csv`, `al_new_points.csv`)
| | 256 runs (iter 0) | 416 runs (iter 10) |
|---|---|---|
| P_f Goodman, 1200 N | 0.9803 | 0.9800 |
| P_f Gerber, 1200 N | 0.3391 | 0.3382 |
| surrogate Q2 (5-fold CV) | 0.9711 | 0.9676 |
| surrogate CV RMSE (MPa) | 3.02 | 3.20 |

- P_f changed by 0.0003 (Goodman) and 0.0009 (Gerber) over 160 added runs, so the earlier P_f values are not an artefact of surrogate error near the boundary.
- Surrogate-vs-FE error at the 160 new boundary points (predicted before refit): max |error| 7.28 MPa = **0.87 % of the FE peak stress at that point** (Goodman points 0.87 %, Gerber points 0.84 %); per-iteration RMSE 2.5-3.9 MPa. FE peaks at these points span 810-894 MPa.

## Gate: surrogate error at boundary points < 1 % of peak stress
**PASS** (max 0.87 %). The margin is small.

## Limitations
- Q2 and CV RMSE did not improve with added points (RMSE ~3 MPa throughout). The error looks like a floor rather than a sampling-density problem. The earlier script's comment attributes it to ~3 MPa of mesh noise; I did not test that here.
- Because of that floor, U cannot reach 2 near g = 0 and the min-U stop rule is not satisfiable in practice: with a continuous MC pool some point lies almost exactly on g = 0 (note the pool is redrawn each iteration, so min U is itself noisy). Logged diagnostics are the mean misclassification probability mean(Phi(-U)) (about 0.002 Goodman, 0.011 Gerber, flat) and the P_f trace.
- The refit uses the same GP hyperparameter bounds as `sensitivity.py`; ConvergenceWarnings about a length scale hitting the 1000 upper bound appeared (an input the stress barely depends on). Not changed.
- The 160 new points were chosen under the surrogate's own P_f, so their error is not an unbiased estimate of surrogate error over the whole input space; the 5-fold CV figures above are the global measure.
- Input spreads (FE inputs, UTS CV 3 %, endurance-ratio SD 0.05, force CV 2 %) remain ASSUMED, as in `fatigue.py`.
- Stale files from the earlier run remain in `results/` (`active_learning_log.csv`, enlarged `fe_samples.csv`, `fatigue_summary_before_al.json`); they were left untouched. `fatigue.py` and `sensitivity.py` still read `results/fe_samples.csv` (288 rows), not the 416-run set in this folder.
