# Phase 2, Step 6 - Goodman life vs limit-state gap

## What was run
`fe/goodman_gap.py` (new; `stage_e`, `fatigue`, `sensitivity` imported unchanged). Log: `run.log`. Samples as in Stage E part (1): GP on all 416 runs, seed 3, 400,000 samples, same life-scatter draw. Variant V1 reproduces the Stage E file (0.9365 and 0.3420, matching `results/stageE/life_at_1200N.csv`).

Variants of the life model:
- **V1** fixed line (Stage E): N = (sa / (A x scale))^(1/b), scale = ratio x UTS / 450.
- **V2** per-sample line through (1e3, anchor_i) and (1e6, allowable_i). **ASSUMED anchor:** anchor_i = 0.9 x 0.45 x UTS_i, an *amplitude*; at UTS 900 it is 364.5 MPa, the Table 2 value. Your instruction said 0.9 x UTS_i. I read that as the Table 2 reading, because 0.9 x UTS as an amplitude (810 MPa) is not Table 2's anchor. This choice does not change classification at 1e6 cycles, only the slope and the life values.
- **V3 (diagnostic, not a proposal)** fixed line as V1, with A recalibrated once so that A x 1e6^b equals the limit-state allowable at the nominal inputs. It separates the endpoint offset from the sample-to-sample dependence.

## 1. Constants and formulas
- Table 2 constants are **fixed**: Goodman A = 670.53, b = -0.088; Gerber A = 549.36, b = -0.059. Material scatter only rescales A per sample (A_i = A x ratio_i x UTS_i / 450); b stays fixed. They are not recomputed from the sample's stress.
- Nominal amplitude: **identical** in both. sa = (1 - R)/2 x peak / 1.66 = 0.45 x peak / 1.66.
- Resistance at 1e6 cycles differs:
  - Life: A x scale x 1e6^b, which is proportional to ratio x UTS, with no Kf and no dependence on the sample's mean stress.
  - Limit state: ratio x UTS / Kf x (1 - sm/UTS) (Goodman) or x (1 - (sm/UTS)^2) (Gerber), with sm = 0.55 x peak / 1.66 and Kf = 1.66.
- So the mean-stress correction is recomputed per sample in the limit state, but is fixed in the life line at whatever stress Table 2 used for its endpoint.

## 4. Line-through-mean check (`line_through_mean_check.csv`)
| | Goodman | Gerber |
|---|---|---|
| A x 1e6^b (MPa) | 198.80 | 243.14 |
| Table 2 value | 198.14 | 241.84 |
| difference | +0.33 % | +0.54 % |
| Limit-state allowable at nominal GP peak 856.2 MPa | 185.64 | 244.15 |
| at paper FE peak 857.8 | 185.48 | 244.05 |
| at paper analytical peak 898.1 | 181.46 | 241.45 |
| limit state vs fixed endpoint at nominal | **-6.6 %** | +0.4 % |

The fixed constants reproduce the Table 2 endpoints to 0.3-0.5 %, not exactly (A and b are rounded). Life at the paper peaks (arithmetic mean of the two criteria): 5.39e5 at the analytical peak against the published 5.75e5, and 1.15e6 at the FE peak against 1e6. Per criterion at the nominal GP peak: Goodman 1.72e5, Gerber 2.20e6. The "181.5 MPa" quoted in the Stage E report is the limit-state allowable at the analytical peak 898.1 (181.46 here), not at the FE peak.

## 2. Pf and classification (`pf_limit_state_vs_life.csv`, `classification_disagreement.csv`)
| | Goodman | Gerber |
|---|---|---|
| (a) limit state P(g<0) | 0.9799 | 0.3361 |
| (b) V1 fixed line P(N<1e6) | 0.9365 | 0.3420 |
| V2 per-sample line | 0.9799 | 0.3361 |
| V3 recalibrated fixed line | 0.9877 | 0.3285 |

Disagreement of V1 with g < 0 (of 400,000):
| | total | g<0 but N>=1e6 | g>=0 but N<1e6 |
|---|---|---|---|
| Goodman | 17,370 | 17,370 | 0 |
| Gerber | 9,203 | 3,424 | 5,779 |

**V2 agrees with the limit state exactly in classification** (0 disagreements for both criteria). This holds by construction: the line passes through (1e6, allowable_i) with negative slope (mean b_i = -0.099 Goodman, -0.059 Gerber; no non-positive allowable), so N < 1e6 exactly when sa > allowable_i.

V2 life, with material scatter (`table12_with_variants.csv`): Goodman median 9.9e4, 5-95 % 2.9e4 to 5.7e5; Gerber median 2.29e6, 5-95 % 1.6e5 to 4.2e8.

## 5. Cause
Under Goodman the fixed line's 1e6-cycle endpoint is higher than the limit-state allowable. At the nominal inputs the Table 2 line gives 198.8 MPa, the limit state 185.6 MPa (-6.6 %), and across the 400,000 samples the mean ratio of fixed endpoint to allowable is 1.072 (Gerber: 0.997). A resistance about 7 % too high at 1e6 cycles removes failures, which gives 0.9365 against 0.9799; all 17,370 disagreements go in that one direction (limit state fails, life survives). Recalibrating only the endpoint (V3) moves the Goodman life Pf to 0.9877, so the endpoint offset accounts for about +0.051 of the 0.043 gap and overshoots by 0.008. The remaining 0.008 comes from the structure of the fixed line (it scales with ratio x UTS and has no dependence on the sample's own mean stress), which the limit state does not share. Under Gerber the endpoints agree to 0.4 % (243.1 against 244.2), so V1 is close to the limit state (0.342 against 0.336). **The Gerber agreement is partly a cancellation:** V1 and the limit state still classify 9,203 samples differently (3,424 one way, 5,779 the other), and V3 is 0.008 below the limit state, so the net Pf match is better than the sample-level match. The table 2 Goodman endpoint (198.14 MPa) is therefore inconsistent with the paper's own Goodman rule (sD x (1 - sm/UTS)) as implemented in `fatigue.py`, while its Gerber endpoint (241.84) is consistent with it. I did not determine why; the paper's Goodman fit may use a different stress or definition. **V2 is the internally consistent variant**, since its life is defined from the same allowable as g.

## 6. Table 12 with variants
`table12_with_variants.csv` (criterion, variant, scatter included, 5th, median, 95th percentile, P(N<1e6), MC SE). The original `results/stageE/life_at_1200N.csv` is untouched.

## Limitations
- The V2 anchor is an ASSUMED reading (above). V2 life values depend on it, though the V2 Pf does not.
- V3 is a diagnostic: the endpoint comes from the limit state itself, so it is not independent of it.
- Life scatter (log10 SD 0.2) is ASSUMED as in Stage E; the UTS CV, ratio SD and force CV are ASSUMED as in Stage D.
- Table 2's equation and stress definitions remain unverified (see the Stage E report); this analysis takes the HTML values as given.
- Gerber V2 95th percentile (4.2e8) is far beyond any tested life and has no physical meaning.

## Files
`line_through_mean_check.csv`, `pf_limit_state_vs_life.csv`, `classification_disagreement.csv`, `table12_with_variants.csv`, `goodman_log10N_hist.png` (600 dpi), `run.log`.
