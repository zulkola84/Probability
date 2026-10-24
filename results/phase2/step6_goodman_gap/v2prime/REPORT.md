# Phase 2, Step 6b - Load-independent per-sample S-N line (V2')

## What was run
`fe/goodman_gap_v2prime.py` (new; existing modules unchanged). Log: `run.log`. Same 400,000 samples as step 6 (GP on all 416 runs, seed 3).

**Formulas** (R = 0.1 load line, sm/sa = 0.55/0.45, sD = ratio x UTS / Kf, Kf = 1.66):
- Goodman: sa* = (ratio x UTS/Kf) / (1 + (0.55/0.45) x ratio/Kf)
- Gerber: sa* = (-1 + sqrt(1 + 4 k sD)) / (2 k), k = sD (0.55/0.45)^2 / UTS^2
- SWT: sa* = sqrt(0.45) x ratio x UTS/Kf (SWT has no Table 2 slope, so only classification is checked)

All three closed forms agree with a numerical root of `ft.g_fun` to better than 3e-11 MPa over 300 random (UTS, ratio) pairs (`nominal_sa_star.csv`).

**V2'**: N = 1e6 x (sa / sa*_i)^(1/b), b fixed at the Table 2 slope (-0.088 Goodman, -0.059 Gerber).

## 2. Nominal check
sa* at UTS 900, ratio 0.5: **Goodman 198.141 MPa, Gerber 241.844 MPa**, against the Table 2 values 198.14 and 241.84. So Table 2's 1e6-cycle stresses are exactly the load-line runout amplitudes at Kf = 1.66, with no tuning.

**Correction to step 6.** There I compared the fixed line's endpoint (198.8) with the limit-state allowable *at the applied 1200 N stress* (185.6) and concluded that the Table 2 Goodman endpoint was inconsistent with the Goodman rule. That was wrong: the two quantities are different things (the allowable depends on the sample's applied mean stress, the Table 2 endpoint does not), and the "6.6 % offset" was an artefact of that comparison. Table 2 is consistent with the rule. The 0.33 % difference between 198.80 and 198.14 comes from rounding of A and b. Step 6's V3 diagnostic therefore tested the wrong thing, and its attribution of the gap to an endpoint offset should be disregarded. Step 6's REPORT.md was left unchanged.

## 3. V2' results (`pf_summary.csv`, `classification_disagreement.csv`, `table12_v2prime.csv`)
| | Goodman | Gerber |
|---|---|---|
| Pf, limit state P(g<0) | 0.9799 | 0.3361 |
| Pf, V2' P(N<1e6) | 0.9799 | 0.3361 |
| Pf, V1 (fixed line) | 0.9365 | 0.3420 |
| V2' disagreements with g<0 | 0 | 0 |
| V2' median life, with material scatter | 1.62e5 | 1.94e6 |
| 5th - 95th percentile | 2.98e4 - 7.11e5 | 1.25e5 - 2.05e7 |

SWT: the criterion `sa > sa*` classifies identically to g < 0 (0 disagreements). Zero disagreements is expected because g decreases monotonically in sa along the load line, so g < 0 exactly when sa > sa*. Table 12 format with the life-scatter rows (log10 SD 0.2, ASSUMED) is in `table12_v2prime.csv`, with the V1 rows beside it.

## 4. Test of the explanation (`endpoint_ratio.csv`, `disagreement_by_ratio.csv`)
V1 endpoint / sa* at UTS 900:
| ratio | 0.4 | 0.5 | 0.6 | 0.65 | mean over samples |
|---|---|---|---|---|---|
| Goodman | 0.949 | 1.003 | 1.057 | 1.084 | 1.003 |
| Gerber | 0.969 | 1.005 | 1.047 | 1.069 | 1.006 |

- V1 is linear in ratio and in UTS. For Goodman, V1/sa* = 0.7334 + 0.5400 x ratio (fit residual 2e-16), because sa* is concave in ratio. The ratio is exactly the same for every UTS (spread 2e-16), so the UTS scaling is **not** the cause; the nonlinearity in the endurance ratio is.
- The sample mean of V1/sa* is 1.003, so the average endpoint is right. The error is a tilt: V1 is too low below the nominal ratio and too high above it.
- **Goodman:** the V1 overstatement matters because the survivors are high-ratio samples. Of the 17,370 disagreements (all of the type g < 0 but N >= 1e6), **16,105 (92.7 %) have ratio > 0.55**, which is 15.9 % of the samples; 17,369 (99.99 %) have ratio > 0.50, and 1 has ratio <= 0.50. Mean V1/sa* is 1.044 in the ratio > 0.55 region.
- **Gerber:** the tilt is milder, and the disagreements are mostly in the other half: 6,588 of 9,203 (71.6 %) have ratio <= 0.50 (V1 too low, so N < 1e6 where g >= 0) and 2,615 (28.4 %) have ratio > 0.50 (only 40 above 0.55). The two errors nearly cancel in Pf (0.3420 against 0.3361), as noted in step 6.

## Verdict
1. **Yes: the Goodman gap is caused by the linear-in-ratio scaling of the fixed line.** Table 2's line passes through the right point at ratio 0.5, but the true runout amplitude is concave in the ratio, so the line overstates the runout amplitude by 5.7 % at ratio 0.6 and 8.4 % at 0.65. Under Goodman most samples fail, and the ones that survive are those with high ratio, which is where the overstatement sits (93 % of the disagreements). Replacing the fixed endpoint by sa*_i removes the gap entirely (Pf 0.9799, 0 disagreements). Your explanation is supported by these numbers.
2. Under Gerber the same nonlinearity exists (smaller: 4.7 % at 0.6) but it acts in both directions around the median, so Pf agrees by cancellation, not because the line is right.
3. **Is V2' the right life model?** It is internally consistent: a true S-N line that depends only on the sample's material (not on the applied load), reduces exactly to Table 2 at the nominal inputs, and gives the limit-state classification for all three criteria. It is the better choice for the paper's probabilistic life than the fixed line. Open points: (a) the slope b is held at the Table 2 value for every sample, which is an ASSUMED simplification (in step 6's V2 the slope varied per sample through the 1e3 anchor, giving the same Pf but different lives); (b) the life scatter is ASSUMED, as before; (c) the Table 2 equation and stress definitions are still unconfirmed against the PDF; (d) the life tails are very wide because of the shallow slopes, so the Gerber 95th percentile (2.0e7) should not be read as a physical life.

## Files
`nominal_sa_star.csv`, `pf_summary.csv`, `classification_disagreement.csv`, `table12_v2prime.csv`, `endpoint_ratio.csv`, `disagreement_by_ratio.csv`, `endpoints_vs_ratio.png` (600 dpi), `run.log`.
