# Phase 2, Step 4 - Elastic-plastic check at the hole edge

## What was run
`fe/plasticity_check.py` (new; `make_mesh`, `stress_field`, `fatigue`, `sensitivity` imported unchanged). Log: `run.log`; CalculiX decks, `.frd`, `.dat` and `.sta` files are in this folder, with names that include `ASSUMED` for the plastic cases.

- **Model:** calibrated quarter model (W = 16.365 mm), hole-edge size 0.3 mm, global 3.0 mm, 23,801 nodes, C3D10, boundary conditions as `make_mesh.write_inp`. Elastic nodal peak at the critical node: sxx = 860.08 MPa at 1200.00 N, equal to `stress_field.py` (860.08).
- **Force scaling:** displacement control on the LOAD line; total force = 4 x |RF_z of the LOAD set| (the quarter model holds a quarter of the total force over both rollers, as in `run_ccx.py`).
- **Material:** E = 105,000 MPa, nu = 0.33, *PLASTIC with isotropic hardening, bilinear from (yield, 0) to (900 MPa, 0.10). **Yield strengths 650, 750, 850 MPa are ASSUMED**; the end point (900 MPa, 0.10) is ASSUMED too. An elastic reference with the same mesh was run.
- **Loading:** 3 cycles, each a load leg to the displacement giving 1200 N and an unload leg to the displacement giving 120 N, 6 increments per leg. Geometry linear, except the NLGEOM check.
- **Critical node:** node 4 at (0, 3.9, 0) mm, the hole edge on the tension face (z = 0), the same node and face as `stress_field.py`.
- **Plastic zone:** outermost node on the net-section path (tension face) with PEEQ > 1e-5 (**ASSUMED** threshold). Resolution is the nodal spacing, about 0.145 mm near the hole.

**Displacement search** (`displacement_search.csv`): at most 4 solves per case, tolerance 1 % of each target. All four cases converged in 3 solves, with F at maximum within 0.01 N of 1200 N and F at minimum within 0.1 N of 120 N (elastic 120.00; 650: 119.90; 750: 119.97; 850: 120.00). The first two solves per case were not wasted by the plasticity: my first correction wrongly shifted the unloading displacement with the peak displacement (the 120 N displacement does not depend on the peak when unloading is elastic), so the second solve missed 120 N (F = 72 N), and the third fixed it.

**Solver times** (final solves, 2 threads each, 4 cases in parallel): elastic 25 s (6 increments); 650, 750 and 850 MPa about 249 s each (36 increments); NLGEOM about 220 s (36 increments). The plastic solves needed no increment cutbacks. **Convergence warnings:** none for the plastic and NLGEOM cases; the elastic deck gives the CalculiX note "PEEQ-output requested" (it has no plasticity, so no PEEQ block was written). My first background run crashed after the solves, in post-processing (missing PEEQ for the elastic case); I fixed this and re-ran the analysis, reusing the identical already-solved decks. The wall times for those reused solves therefore come from file time stamps, not from a timer, and include the effect of running four jobs in parallel.

## NLGEOM check (750 MPa, ASSUMED; `nlgeom_check_yield750_ASSUMED.csv`, `nlgeom_amplitude_mean_yield750_ASSUMED.csv`)
Same displacements as the linear-geometry solve:
| | linear geometry | NLGEOM | difference |
|---|---|---|---|
| force at maximum (N) | 1200.00 | 1188.86 | -0.93 % |
| sxx at maximum (MPa) | 753.53 | 753.66 | +0.02 % (same displacement) |
| sxx at maximum, scaled to the same force | 753.53 | 760.72 | +0.95 % |
| amplitude, cycle 3 (MPa) | 387.05 | 384.25 | -0.72 % |
| mean, cycle 3 (MPa) | 366.48 | 369.41 | +0.80 % |

Geometric nonlinearity changes the local stress by about 1 %, much smaller than the effect of yielding (below), so it is not the issue. The remaining results use linear geometry.

## Critical node, cycles 1 and 3 (`critical_node_table.csv`)
| | elastic | 650 (ASSUMED) | 750 (ASSUMED) | 850 (ASSUMED) |
|---|---|---|---|---|
| sxx at max load, cycle 1 / 3 (MPa) | 860.08 / 860.08 | 655.65 / 655.65 | 753.53 / 753.53 | 860.08 / 860.08 |
| sxx at min load, cycle 1 / 3 (MPa) | 86.01 / 86.01 | -118.50 / -118.50 | -20.57 / -20.57 | 86.01 / 86.01 |
| von Mises at max, cycle 1 / 3 | 858.96 | 656.51 | 752.49 / 752.49 | 858.96 |
| von Mises at min, cycle 1 / 3 | 85.90 | 116.67 | 20.85 / 20.85 | 85.90 |
| amplitude, cycle 3 (MPa) | 387.04 | 387.07 | 387.05 | 387.04 |
| mean, cycle 3 (MPa) | 473.04 | 268.58 | 366.48 | 473.04 |
| local R, cycle 3 | 0.100 | -0.181 | -0.027 | 0.100 |
| peak PEEQ after cycle 1 / 3 | 0 | 0.00238 / 0.00238 | 0.00119 / 0.00119 | 0 / 0 |
| plastic zone along the net section (mm) | none | 0.58 | 0.145 | none |

- **Shakedown:** yes. Stresses and PEEQ in cycle 3 equal those of cycle 1 (to the precision of the table), so plasticity occurs only in the first loading and later cycles are elastic. The local stress range (about 774 MPa) stays below twice any assumed yield, so there is no reverse yielding and the choice of isotropic hardening does not matter here.
- **What yielding does:** it lowers the local maximum (by 23.8 % at yield 650, 12.4 % at 750) and, because the unloading is elastic, shifts the whole cycle down by the same amount: the **amplitude is unchanged** (+0.009 %, +0.003 %), and the **mean drops by 43.2 % (650) and 22.5 % (750)**. At 850 nothing changes.
- **850 MPa:** the nodal von Mises (859 MPa) exceeds the assumed yield, yet PEEQ is zero. Nodal values are extrapolated from the integration points and the integration-point stress is lower, so 850 MPa is the onset: whether yielding starts depends on this mesh. I did not check the integration-point values.
- The plastic zone size is limited by the node spacing (the next node beyond the zone is at 0.29 mm for 750 and 0.72 mm for 650).

## Fatigue consequence (exploratory, **ASSUMED** treatment; `pf_exploratory_ASSUMED.csv`)
f = (cycle-3 sxx at 1200 N, elastic-plastic) / (elastic) = 0.762 (650), 0.876 (750), 1.000 (850). Baseline limit state (GP on 416 runs, 400,000 samples, seed 1, Kf = Kt = 1.66) gives 0.9800 / 0.3365 / 0.9932.

- **V_f:** the peak stress of every sample multiplied by f (this scales amplitude and mean together).
- **V_loc:** amplitude and mean multiplied by the cycle-3 local ratios (amplitude x1.0001 / 1.0000, mean x0.568 / 0.775). Equivalent to using the cycle-3 local values divided by Kt in place of the elastic split, applied as ratios to each sample's elastic amplitude and mean.

| yield (ASSUMED) | criterion | baseline | V_f | V_loc |
|---|---|---|---|---|
| 650 | Goodman | 0.9800 | 0.1062 | 0.6548 |
| | Gerber | 0.3365 | 0.0018 | 0.1451 |
| | SWT | 0.9932 | 0.4071 | 0.8545 |
| 750 | Goodman | 0.9800 | 0.6313 | 0.8765 |
| | Gerber | 0.3365 | 0.0424 | 0.2102 |
| | SWT | 0.9932 | 0.8620 | 0.9616 |
| 850 | all | unchanged | unchanged | unchanged |

V_f overstates the benefit because it also lowers the amplitude, which the FE shows does not change. V_loc follows the FE result (same amplitude, lower mean) and is the better of the two. **This is not an elastic-plastic fatigue model:** it keeps the elastic S-N and Goodman/Gerber/SWT rules, applies one nominal-case ratio to every Monte Carlo sample (a plate with a higher peak would yield more, a lower one less), and ignores plastic strain damage and any change in the endurance limit.

## Verdict, per assumed yield strength
- **650 MPa:** the elastic peak is **conservative**: local peak -23.8 %, amplitude unchanged, mean -43 %. With V_loc Pf falls (Goodman 0.980 to 0.655, Gerber 0.337 to 0.145, SWT 0.993 to 0.855).
- **750 MPa:** conservative: peak -12.4 %, amplitude unchanged, mean -22.5 %. V_loc: Goodman 0.877, Gerber 0.210, SWT 0.962.
- **850 MPa:** no yielding at this mesh, no change. The elastic analysis is exact in this case.
- In no case did yielding raise the amplitude, so the elastic peak is not unconservative for the fatigue limit state in these cases; the open question is only how large the conservatism is, and that depends on the unpublished yield strength. Residual compressive stress at the hole edge (sxx = -118 MPa at 120 N for 650) also lowers the local R below 0.1, which the limit state, built on R = 0.1, does not represent.

## ASSUMED / limitations
- ASSUMED: yield strengths 650/750/850 MPa, the (900 MPa, 0.10) end of the bilinear curve, the PEEQ threshold 1e-5, 6 increments per leg, and the one-ratio-for-all-samples treatment in the Pf rerun.
- Nodal values are CalculiX extrapolated values; the 0.3 mm mesh was not re-converged for the plastic cases. Plastic zone sizes are resolved only to the node spacing (about 0.15 mm).
- Isotropic hardening and rate-independent plasticity only; Bauschinger and cyclic softening were not tested (no reverse yielding occurred in these runs).
- Stage C's stiffness gap and the rectangle-section idealisation still apply: the absolute stress is that of the calibrated model.

## Files
`displacement_search.csv`, `critical_node_table.csv`, `plastic_zone.csv`, `net_section_path.csv`, `crit_node_history_*.csv`, `nlgeom_check_yield750_ASSUMED.csv`, `nlgeom_amplitude_mean_yield750_ASSUMED.csv`, `pf_exploratory_ASSUMED.csv`, `meta.json`, `run.log`, `fig_sxx_net_section_ASSUMED_yields.png` (600 dpi; the elastic curve lies under the 850 MPa curve because they are identical), `fig_stress_strain_loops_ASSUMED_yields.png` (600 dpi), and the CalculiX files.
