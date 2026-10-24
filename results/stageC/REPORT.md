# Stage C - Curved (circular-arc) section and rigid rollers

## What was run
- `fe/arc_section.py`: concentric arc band, radial thickness 5.12 mm, outer radius R_out = 25 mm (**ASSUMED**), half-angle solved so the NET section at the hole (|y| < 3.9 mm removed) has I = 105.72 mm^4 about its own centroid. This is the same net-section reading of the paper's I used in `section.py`. Result: alpha = 21.9 deg, bounding-box width 18.66 mm. Width is a consequence of I, not tuned to stress or stiffness.
- `fe/stage_c.py`, `fe/stage_c_run.py`: quarter model, 10-node tets, hole refined to 0.55 mm; two rigid cylindrical rollers (radius 5 mm, **ASSUMED**), axis along y, frictionless node-to-surface penalty contact (K = 1e6 MPa/mm, **ASSUMED** numerical parameter), displacement control on the load roller. Total force = 4 x load-roller reaction. Orientation **ASSUMED**: concave side up (supports touch the apex of the convex back; load rollers touch the two inner edges). Hole drilled along z (**ASSUMED**). E = 105 GPa, nu = 0.33 as before.
- 13 FE runs (about 65 s each, 4 in parallel). Data: `stageC_runs.csv`, log `run.log`, beam-theory inverse `required_I.txt`.
- Check: load-roller and support reactions agree to all printed digits in every run.

## Key numbers
Force-displacement (base case), linear fit over 1-3 mm as for the experiment:

| | K (N/mm) |
|---|---|
| FE, arc section + rollers | **541.8** (secant 554.7 / 548.4 / 546.1 at 1 / 2 / 3 mm) |
| Experiment, upper / middle / lower strand | 393.6 / 361.5 / 363.5 |
| Paper FE (blue) | about 372 |
| Previous rectangular model | 471 |

FE is +37.7 % (upper strand) to +49.9 % (middle strand) above experiment.

Peak stress at the 1200 N load (run at 2.191 mm displacement, FE force 1200.5 N): max von Mises on the hole surface **942.1 MPa**, i.e. **+9.8 %** vs the paper's FE 857.8 MPa (+4.9 % vs the analytical 898.1). No gate was set on stress.

Mesh convergence (hole-edge size, run at 2 mm, peak scaled to 1200 N): 1.2 mm 935.5, 0.8 mm 938.2, 0.55 mm 941.1, 0.4 mm 941.7 MPa (last step +0.06 %); stiffness varies by under 1 % between meshes.

Sensitivity (each at 2 mm; K in N/mm, peak scaled to 1200 N):
| change | K | peak (MPa) |
|---|---|---|
| base | 548.4 | 941.1 |
| contact stiffness 1e5 / 1e7 | 547.3 / 549.9 | 941.4 / 938.8 |
| convex side up | 553.2 | 943.3 |
| R_out = 15 mm | 595.3 | 1065.2 |
| R_out = 60 mm | 527.5 | 829.0 |
| I = 105.72 read as GROSS section | 237.4 | 5722.7 (degenerate) |

## Gate: stiffness within 10 % of experiment
**FAIL** (+38 % to +50 %). As the gate failed, `run_samples.py`, `sensitivity.py` and `fatigue.py` were NOT rerun into `results/stageC/`, per the brief.

## Findings and limitations
- Roller contact compliance is not the missing flexibility: contact stiffness changes K by under 0.3 %, orientation by under 1 %, and the FE K is within 4 % of the beam-theory value for the gross section (arc_section.py prints 560 N/mm at R = 25 mm).
- Beam theory (`required_I.txt`) says the experiment needs an effective I of 131-143 mm^4, between the net-section I (105.72) and the gross-section I of the arc (about 203). No arc radius tried gives that: K stays between 527 and 595 N/mm for R_out 60 to 15 mm.
- The alternative reading (I = 105.72 as the gross section) gives a 10.5 mm wide plate with a 1.3 mm ligament each side of the hole, a degenerate stress result, and K = 237 N/mm (34-40 % below experiment). It also fails the gate and is not physically consistent with the paper's 16.7 mm bounding box.
- The peak stress is sensitive to the assumed curvature (829-1065 MPa for R_out 60-15 mm); +9.8 % at R = 25 mm is therefore not a validated agreement, only the value for one assumed radius. A stepped, non-circular outline like the paper's Table 2 section was not modelled.
- Untested hypotheses for the stiffness gap: the experimental displacement includes fixture/machine compliance and plate seating that this model excludes; the paper's E (calibrated against static tests) or its real section differs from this idealisation. Neither was checked here.
- Contact is a penalty method with node-to-surface slaves on a 0.6 mm local mesh; edge contact on the sharp inner corners is a point-like contact. Single runs per case, no repeats; peak stress read from CalculiX extrapolated nodal values.
- Peak stress at 1200 N is from the run whose force is 1200.5 N (scaled value 941.7 MPa, unscaled 942.1).
