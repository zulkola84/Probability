# Stage 3 report: Monte Carlo sensitivity and fatigue reliability

Plate: Bologna et al. 2023 (Frontiers in Medical Engineering), grade 4 titanium diaphyseal femoral plate, ASTM F382 four-point bending, R = 0.1, experimental runout 1,200 N.
Units: mm, N, MPa. Scripts are in `fe/`; all numbers below come from files in `results/`.

## 1. Summary

- A CalculiX four-point bending model of the plate converges on mesh (peak von Mises stress changes by under 0.4% from a 0.8 mm to a 0.3 mm hole-edge size) and reproduces the paper's FE peak stress of 857.8 MPa after one fitted parameter, the plate width (16.365 mm).
- Propagating assumed manufacturing and loading scatter through a Gaussian-process surrogate (cross-validated Q² = 0.97 for stress, 1.00 for stiffness) gives a peak stress of 857 ± 25 MPa and a stiffness of 472 ± 28 N/mm.
- Sobol analysis: peak stress is driven by the applied-force accuracy and the thickness tolerance; stiffness by Young's modulus and thickness.
- Fatigue reliability at 1,200 N: the predicted failure probability is 0.98 (Goodman) and 0.34 (Gerber). The nominal-plate runout load predicted by the two criteria is 1,024 N and 1,250 N, which bracket the experimental 1,200 N. The endurance ratio, an assumed input, drives this result (Sobol index 0.70–0.80), not the FE model.
- Active learning near the failure boundary (32 extra FE runs) moved the failure probabilities by less than 0.001. The surrogate is not the limiting factor.
- Two things do not agree with the paper and are reported openly: the model is about 24% too stiff against the digitised Figure 6 (471 vs 363–394 N/mm), and a section calibrated to the paper's moment of inertia gives 8% too little stress and about 30% too much stiffness. The independent implanted-plate validation could not be done.

## 2. Model and inputs

| Item | Value | Source |
|---|---|---|
| Thickness t | 5.12 mm | paper |
| Hole diameter d | 7.8 mm | paper |
| Inner (loading) roller spacing | 80 mm | paper |
| Lever arm, inner roller to support | 38 mm (support span 156 mm) | paper |
| Young's modulus | 105 GPa | paper (Table 1) |
| Poisson's ratio | 0.33 (paper's E and G imply 0.3125) | assumed |
| Tensile strength | 900 MPa | paper (certificate) |
| Plate width | 16.365 mm | fitted so the FE von Mises peak at 1,200 N is 857.8 MPa |
| Plate length | 180 mm | assumed |
| Rollers | line loads, linear analysis | idealisation |

The model is a quarter plate (symmetry at mid-length and mid-width) with 10-node tetrahedra, refined toward the hole edge. The spans were settled by a check against the paper: M·c/I with I = 105.72 mm⁴ and a 38 mm lever arm gives 552 MPa nominal, and times the paper's Kt of 1.66 that is 916 MPa against the paper's 898 MPa analytical peak. The reverse reading (80 mm lever arm) gives about 1,160 MPa nominal and does not fit. The 2022 conference paper's runout moment of 22,800 N·mm equals 600 N × 38 mm, which confirms it.

Mesh convergence (`results/convergence.csv`, width 16.365 mm):

| Hole-edge size (mm) | Nodes | Peak von Mises (MPa) | Stiffness (N/mm) |
|---|---|---|---|
| 1.2 | 3,316 | 855.8 | 471.6 |
| 0.8 | 4,813 | 855.4 | 471.4 |
| 0.55 | 7,895 | 858.5 | 471.4 |
| 0.4 | 13,404 | 857.7 | 471.4 |
| 0.3 | 23,801 | 859.1 | 471.4 |

Mesh noise is about ±1.5 MPa. All later sampling used the 0.55 mm size.

## 3. Digitised experiment (Figure 6)

Digitised from the published image with calibrated axes (`inputs/fig6_*.csv`). Only three red strands are separable: a thick upper strand (probably two overlapping curves) and two thin lower ones.

| Curve | Stiffness (N/mm), linear fit 1–3 mm |
|---|---|
| Experimental, upper strand | about 394 |
| Experimental, middle strand | about 363 |
| Experimental, lower strand | about 363 |
| Paper FE (blue dashed) | about 372 |

Pixel digitising is good to a few N/mm at best.

## 4. Monte Carlo sensitivity (`sensitivity.py`)

Uncertain inputs (`uncertainty.py`; **all spreads are assumptions**, the paper gives no tolerances): thickness N(5.12, 0.05) mm; hole diameter ±0.1 mm uniform; E with 5% sd; ν uniform 0.30–0.36; lever arm ±0.5 mm uniform; applied force with 2% sd. Width is held fixed because it is a fitted constant, not a measured tolerance; giving it a spread made the ranking reflect my own choice of range.

The FE model was run at a 256-point scrambled Sobol-sequence design, a Gaussian-process surrogate was fitted (5-fold cross-validated Q² = 0.971 and RMSE 3.0 MPa for stress; Q² = 1.000 for stiffness), and 8,192 base samples were propagated with SALib.

| Output | Mean | SD | 5th–95th percentile |
|---|---|---|---|
| Peak von Mises (MPa) | 857.2 | 24.6 | 817–898 |
| Stiffness (N/mm) | 471.6 | 28.4 | 426–519 |

The probability that the peak stress exceeds the 900 MPa tensile strength is 4.3%.

First-order Sobol indices (total indices are equal within the confidence intervals, so effects are additive):

| Input | Peak stress | Stiffness |
|---|---|---|
| Force | 0.49 | 0 |
| Thickness | 0.41 | 0.23 |
| Lever arm | 0.08 | 0.08 |
| Hole diameter | 0.02 | 0 |
| E | 0 | 0.69 |
| ν | 0 | 0 |

Stress does not depend on E because the load is imposed. The small hole-diameter index is within the noise of the surrogate. Plot: `results/sensitivity.png`.

## 5. Fatigue reliability (`fatigue.py`)

Limit state, following the paper: nominal maximum stress = peak / 1.66; at R = 0.1 the amplitude and mean are 0.45 and 0.55 times that; the fatigue limit is 0.5 × UTS / 1.66; Goodman gives the lower bound and Gerber the upper bound of the allowable stress. Failure (no runout at 10⁶ cycles) is g < 0. Extra assumed scatter: tensile strength with a 3% coefficient of variation, and the endurance ratio N(0.5, 0.05).

| | Goodman | Gerber |
|---|---|---|
| Failure probability at 1,200 N (Monte Carlo, 400,000 samples) | 0.980 (95% CI 0.980–0.981) | 0.338 (0.337–0.339) |
| FORM check | 0.979 | 0.333 |
| Predicted runout load, nominal plate | 1,024 N (−14.7%) | 1,250 N (+4.1%) |
| Dominant input (Sobol ST of g) | endurance ratio 0.70 | endurance ratio 0.80 |

The two criteria bracket the experimental runout of 1,200 N, which is a consistency check on the limit state. The high Goodman probability reflects that this criterion is conservative by construction; it is not a prediction that the plate fails. Plot: `results/fatigue.png`.

## 6. Active learning near the failure boundary (`active_learning.py`)

Each round, candidate points were ranked by U = |g| / σ_g (AK-MCS learning function) for both criteria, a space-filling batch was run in the FE model, compared with the surrogate's prediction before retraining, and added to the training set. After 2 rounds (16 FE runs each, 288 training points in total) the failure probabilities changed by less than 0.001, which was the stop rule.

| Criterion | Failure probability before | After | Change |
|---|---|---|---|
| Goodman | 0.9803 | 0.9801 | −0.0002 |
| Gerber | 0.3382 | 0.3376 | −0.0006 |

These changes are inside the Monte Carlo uncertainty of about ±0.0005. At the new boundary points the surrogate matched FE with an RMSE of 2.5–4.0 MPa and no systematic bias; that is the mesh-noise level, so more samples cannot improve it. The criterion min U ≥ 2 cannot be used: a continuous sample pool always contains a point on g = 0, so the expected misclassification rate is logged instead (about 0.2% Goodman, 1.1% Gerber). Log: `results/active_learning_log.csv`; the original training set is kept as `results/fe_samples_initial.csv`.

## 7. Realistic cross-section (`section.py`)

The paper's I = 105.72 mm⁴ is read as the net section through the hole: its text describes the critical section at a hole as two separate polylines, and Table 1's plate area (51.21 mm²) is a net value. Both a rectangle and an ellipse were calibrated to that net I (`results/section_comparison.csv`):

| Section | Width (mm) | Gross I (mm⁴) | FE peak von Mises (MPa) | FE stiffness (N/mm) | Beam theory, gross I (N/mm) |
|---|---|---|---|---|---|
| Rectangle | 17.25 | 193 | 786 | 500 | 533 |
| Ellipse | 28.8 | 190 | 779 | 503 | 524 |
| Paper / experiment | – | – | 857.8 FE; 898.1 analytical | 363–394 | – |

- **The stiffness gap is open.** Both shapes are about 30% too stiff. Beam theory with the gross I reproduces the FE stiffness, so the hole and the mesh are not the cause. Any section with the right net I has too much gross I. Plausible causes, none tested here: compliance of the load frame and rollers in the measured displacement, and the plate's longitudinal curvature.
- **The stress gap** is about 8% (stress concentration against the net-I nominal stress of about 1.42, against about 1.55 in the paper's FE). This is consistent with the paper's note that curvature raises the stress concentration.
- The ellipse needs an unrealistic width (28.8 mm against the paper's 16.7 mm for the similar Plate 11) and is only a check, not a model to use.
- The fitted 16.365 mm rectangle matches the stress by construction and is 24% too stiff.

## 8. Not done, and limitations

1. **Implanted two-hole plate validation (paper: 811.7 MPa FE, 768.6 MPa analytical, I = 165.89 mm⁴).** Not done. The stress depends on hip and muscle forces for a 100 kg patient at 45% of the gait cycle, the femur geometry and the hole pitch, none of which are published. The paper's supplementary file only gives that model's mesh convergence (757 MPa at 1.0 mm to 815 MPa at 0.1 mm). With the loads and distances this check could be built; without them any match would be invented.
2. **All input distributions are assumptions.** The paper gives no tolerances. The failure probabilities depend most on the assumed endurance-ratio scatter and the tensile-strength scatter.
3. **Plate geometry.** The real section is curved and stepped; length, hole count and hole pitch are unknown. Width is a fitted stand-in.
4. **Idealisations.** Line loads in place of contact rollers, linear analysis (the experimental curves show slight curvature), a flat plate with no longitudinal curvature.
5. **Digitisation.** Three of the four experimental curves are separable; values are good to a few N/mm.
6. **Stress measure.** The model reports von Mises stress, as the 2023 paper does. The 2022 conference abstract reports 1,025 MPa for a similar quantity, which the 2023 paper's values do not reproduce, and I have not resolved why.

## 9. Files

- Scripts: `fe/params.py`, `make_mesh.py`, `run_ccx.py`, `convergence.py`, `uncertainty.py`, `run_samples.py`, `sensitivity.py`, `fatigue.py`, `active_learning.py`, `section.py`, `stiffness_check.py`.
- Results: `results/convergence.csv`, `fe_samples.csv` (288 points after active learning), `fe_samples_initial.csv` (256 points), `sensitivity_summary.csv`, `sobol_indices.csv`, `sensitivity.png`, `fatigue_summary.json`, `fatigue_summary_before_al.json`, `fatigue_sobol_{goodman,gerber}.csv`, `fatigue.png`, `active_learning_log.csv`, `section_comparison.csv`.
- Inputs: `inputs/bologna_inputs.md`, `inputs/fig6_*.csv`, `inputs/annals2023_table2_cross_sections.png`.

Note: `sensitivity_summary.csv`, `sobol_indices.csv` and `sensitivity.png` come from the 256-point run before active learning. The fatigue files come from the final 288-point set.
