# Stage F - Surrogate benchmark

## What was run
`fe/stage_f.py` (about 40 s), project `.venv` (scikit-learn, OpenTURNS 1.27; nothing installed). The `BadmintonPose` environment you activated lacks SALib, gmsh and OpenTURNS, so it was not used.
- Data: 416 FE runs (`results/stageA/fe_samples_al.csv`: 256 Sobol DOE + 160 active-learning boundary points), 5 inputs, output peak von Mises at 1200 N.
- Split: random 20 % hold-out, stratified by source, seed 0: 332 train, 84 test (52 DOE, 32 active-learning). All tuning used 5-fold CV on the training part only.
- Models: **GP** (kernel of `sensitivity.py`); **neural network** (scikit-learn MLP, tanh, L-BFGS; grid over hidden sizes (16), (32,32), (64,64) and L2 alpha 1e-4, 1e-2, 1 -> chosen (64,64), alpha 1); **gradient boosting** (scikit-learn; grid over trees, depth, learning rate, subsample -> 600 trees, depth 2, rate 0.05, subsample 0.7); **polynomial chaos** (OpenTURNS, Hermite/Legendre basis on the input distributions, least squares or LARS, degree 2-6 -> dense least squares, degree 2).
- Files: `test_scores.csv`, `test_predictions.csv`, `test_rmse_vs_gp_bootstrap.csv`, `pf_by_model.csv`, `gp_pf_uncertainty.csv`, `model_selection.json`, `pce_cv.csv`, `benchmark.png`.

## Hold-out scores (84 test points)
| model | R2 | RMSE (MPa) | max abs error (MPa, % of local peak) | training CV RMSE (MPa) |
|---|---|---|---|---|
| GP | 0.9735 | 2.89 | 7.22 (0.86 %) | n/a (see below) |
| Neural network | 0.9730 | 2.92 | 7.17 (0.85 %) | 3.32 |
| Gradient boosting | 0.9575 | 3.66 | 9.32 (1.12 %) | 3.85 |
| Polynomial chaos | 0.9725 | 2.94 | 7.13 (0.85 %) | 3.28 |

By subset the picture is the same: DOE points RMSE 2.3-2.4 MPa (GB 3.26), active-learning points 3.6 (GB 4.23) for the three leading models.

Paired bootstrap of the test RMSE difference against the GP (2000 resamples, 95 % interval): neural network +0.03 (-0.07 to +0.13), PCE +0.06 (-0.10 to +0.20), gradient boosting +0.77 (+0.29 to +1.25). **GP, neural network and PCE are statistically indistinguishable on this test set; gradient boosting is worse.**

## Failure probability with each surrogate (400,000 common MC inputs, assumed scatter, Kf = Kt)
| model | P_f Goodman | P_f Gerber | P_f SWT |
|---|---|---|---|
| GP | 0.9800 | 0.3362 | 0.9932 |
| Neural network | 0.9799 | 0.3359 | 0.9932 |
| Gradient boosting | 0.9801 | 0.3357 | 0.9933 |
| Polynomial chaos | 0.9799 | 0.3361 | 0.9932 |

All four agree to within 0.0005, which is inside the MC standard error (0.0002 Goodman, 0.0007 Gerber, 0.0001 SWT). The surrogate type does not change the reliability result, even though gradient boosting has a visibly larger test error: errors of this size (about 3 MPa) are small against the distance of most MC samples from g = 0, and their mean on the test set is small (-0.4 to -0.7 MPa, all four models slightly under-predicting; GP -0.38), so they largely average out in P_f. A -0.4 MPa shift in stress is a much smaller effect than the Stage D/E assumptions; I did not quantify its effect on P_f.

## GP predictive variance in the P_f confidence interval (`gp_pf_uncertainty.csv`)
Epistemic uncertainty: 300 joint posterior draws of the latent GP (white-noise kernel excluded; the latent mean was checked equal to `predict`) on a 5,000-point subset of the MC pool; P_f recomputed per draw; its standard deviation is added in quadrature to the MC error. 95 % interval = plug-in +- 1.96 x sqrt(MC^2 + epistemic^2).

| GP | criterion | P_f | MC SE | epistemic SD | 95 % interval |
|---|---|---|---|---|---|
| trained on 332 | Goodman | 0.9800 | 0.0002 | 0.0004 | 0.9792-0.9808 |
| trained on 332 | Gerber | 0.3362 | 0.0007 | 0.0011 | 0.3336-0.3389 |
| trained on 332 | SWT | 0.9932 | 0.0001 | 0.0001 | 0.9929-0.9935 |
| trained on all 416 | Goodman | 0.9800 | 0.0002 | 0.0003 | 0.9793-0.9808 |
| trained on all 416 | Gerber | 0.3365 | 0.0007 | 0.0010 | 0.3341-0.3389 |
| trained on all 416 | SWT | 0.9932 | 0.0001 | 0.0001 | 0.9930-0.9935 |

The latent GP has a mean standard deviation of 0.47-0.54 MPa, against 3.2-3.3 MPa for the full predictive standard deviation, so almost all the predictive spread is the fitted noise term, which averages out over the MC pool. A per-point probabilistic version (using the total predictive SD) gives the same P_f (0.9799 / 0.3363-0.3366 / 0.9932).

## Limitations
- The 84-point test set limits what can be concluded: only the gradient-boosting deficit is resolved; ranking among GP, neural network and PCE is not. Single split, no repeated cross-validation of the comparison.
- Errors of every model sit near 3 MPa. The GP's fitted noise level corresponds to about 3 MPa and the best PCE is only quadratic, which is consistent with irreducible FE scatter (each sample is remeshed) rather than a surrogate limitation; I did not test this directly.
- Neural-network tuning ended at the grid edge on alpha (1.0) and size ((64,64)); GB and PCE had small grids. Tuning effort was modest and equal-ish, not exhaustive. TensorFlow/PyTorch networks and other kernels were not tried.
- The epistemic interval assumes the fitted GP hyperparameters and stationary kernel are right; it does not cover model-form error, FE-model error (Stage C), assumed input distributions, or the criterion/Kf choice, all of which move P_f by far more (Stages D, E). The width is also estimated on a 5,000-point subset.
- Active-learning points are not a random sample of the input space, so subset scores apply to the boundary region only.
