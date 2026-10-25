# Phase 2, Step 5 - Tensile strength 900 vs 991 MPa

## What was run
`fe/uts_cases.py` (new; `fatigue`, `goodman_gap_v2prime`, `stage_e`, `sensitivity` imported unchanged). Log: `run.log`. GP on all 416 runs, Stage D "assumed" scatter for all other inputs, Kf = Kt = 1.66, 400,000 common-random-number samples (seed 1), Sobol base sample 8192 with 8 variables.

| case | UTS mean | UTS CV |
|---|---|---|
| A (control) | 900 | 3 % |
| B | 991 | 3 % |
| C | 945.5 (midpoint) | **5 % ASSUMED** |

Control A reproduces Stage D: Pf 0.9800 / 0.3365 / 0.9932 (Goodman / Gerber / SWT).

## 3. Which UTS reproduces Table 2 (`table2_vs_uts.csv`)
Load-line runout amplitude sa* (ratio 0.5, Kf = 1.66):
| | Goodman | Gerber | SWT |
|---|---|---|---|
| Table 2 value | 198.14 | 241.84 | n/a |
| sa* at UTS 900 | **198.141** | **241.844** | 181.85 |
| sa* at UTS 991 | **218.18** | **266.30** | 200.24 |
| UTS implied by Table 2 | 899.996 | 899.987 | |
| ratio needed at UTS 991 to give Table 2 | 0.439 | 0.445 | |

UTS = 900 reproduces both Table 2 amplitudes exactly. The Table 2 anchor at 1e3 cycles, 364.50 MPa, also implies UTS = 900.0 (364.50 / (0.9 x 0.45)). At 991 the amplitudes would be 218.18 (Goodman) and 266.30 (Gerber). This shows that the paper's fatigue calculation used 900 (with ratio 0.5 and Kf 1.66). It says nothing about the true strength of the tested plate. I did not check the two later papers that list 991 MPa. The Stage D report notes that the Bologna certificates list 668 and 991 MPa for two different pure-titanium plates, which I took from that report, not from the source.

## 1-2. Results (`pf_form_runout.csv`)
| case | criterion | Pf (MC SE) | FORM Pf (beta) | runout load (N, vs 1200) |
|---|---|---|---|---|
| A 900 | Goodman | 0.9800 (0.0002) | 0.9782 (-2.017) | 1024 (-14.6 %) |
| | Gerber | 0.3365 (0.0007) | 0.3309 (0.437) | 1250 (+4.2 %) |
| | SWT | 0.9932 (0.0001) | 0.9928 (-2.446) | 940 (-21.7 %) |
| B 991 | Goodman | 0.7830 (0.0007) | 0.7760 (-0.759) | 1128 (-6.0 %) |
| | Gerber | 0.0820 (0.0004) | 0.0807 (1.400) | 1377 (+14.7 %) |
| | SWT | 0.9269 (0.0004) | 0.9244 (-1.435) | 1035 (-13.7 %) |
| C 945.5, CV 5 % ASSUMED | Goodman | 0.8942 (0.0005) | 0.8882 (-1.217) | 1076 (-10.3 %) |
| | Gerber | 0.1977 (0.0006) | 0.1929 (0.867) | 1314 (+9.5 %) |
| | SWT | 0.9644 (0.0003) | 0.9626 (-1.781) | 988 (-17.7 %) |

FORM converged in all nine cases. (Negative beta means the nominal plate already fails, the convention of `fatigue.py`.) Runout loads scale with the mean UTS because the limit state is proportional to UTS at a fixed ratio.

Total-order Sobol indices (`sobol_by_case.csv`), Goodman:
| variable | A | B | C |
|---|---|---|---|
| endurance ratio | 0.696 | 0.734 | 0.581 |
| UTS | 0.135 | 0.131 | 0.298 |
| force | 0.082 | 0.065 | 0.059 |
| thickness | 0.070 | 0.056 | 0.051 |

Gerber and SWT show the same pattern (all values in the CSV). Case C raises the UTS index because its CV is 5 %.

## 4. Life at 1200 N, V2' line (`life_v2prime_by_case.csv`)
| case | criterion | median | 5 % | 95 % |
|---|---|---|---|---|
| A | Goodman | 1.62e5 | 2.96e4 | 7.09e5 |
| A | Gerber | 1.94e6 | 1.23e5 | 2.04e7 |
| B | Goodman | 4.83e5 | 8.85e4 | 2.12e6 |
| B | Gerber | 9.92e6 | 6.29e5 | 1.04e8 |
| C | Goodman | 2.79e5 | 4.40e4 | 1.48e6 |
| C | Gerber | 4.36e6 | 2.27e5 | 5.99e7 |

P(N < 1e6) from V2' equals the limit-state Pf in every case, with 0 classification disagreements. The Gerber 95th percentiles are far beyond any tested life and should not be read as physical.

## 5. Do the conclusions change with UTS?
- **(a) Ranking of Pf across criteria: unchanged.** SWT > Goodman > Gerber in all three cases. The Goodman-Gerber spread is 0.644 (A), 0.701 (B) and 0.697 (C).
- **(b) Does 1200 N lie between the Goodman and Gerber runout loads? Yes in all cases:** 1024-1250 N (A), 1128-1377 N (B), 1076-1314 N (C). SWT's runout load is below 1200 N in all three (940, 1035, 988 N).
- **(c) Dominant input: unchanged.** The endurance ratio ranks first and UTS second in all nine case-criterion combinations (force third). The share of the endurance ratio rises with 991 (0.734 against 0.696 for Goodman).
- **What does change is the level:** at UTS 991 the Goodman Pf drops from 0.980 to 0.783, the Gerber Pf from 0.337 to 0.082, and the predicted lives rise by a factor of 3 (Goodman) to 5 (Gerber). Which criterion fits the observed 1200 N runout best flips: at 900 Gerber is closer (+4.2 % against -14.6 %), at 991 Goodman is closer (-6.0 % against +14.7 %).

## Verdict
The qualitative conclusions (criterion ordering, bracketing of 1200 N, dominance of the endurance ratio) do not depend on whether UTS is 900 or 991 MPa. The numbers do: a 10 % higher UTS cuts the failure probabilities sharply and moves the predicted runout loads up by about 10 %, so the paper should report results for both or state clearly which UTS it uses. Table 2 itself is consistent with 900 only, so using 991 means that Table 2's A and b would not apply unchanged. The question of which value applies to the tested plate remains open: this analysis cannot settle it.

## ASSUMED
- UTS CV 5 % in case C and the midpoint idea itself; CV 3 % in A and B, taken from Stage D.
- Endurance-ratio SD 0.05, force CV 2 %, the FE-input ranges, and a lognormal UTS (all from Stage D); Kf = Kt.
- The Table 2 slopes (fixed, in V2') and its equation (unverified, as before).
- Case B applies the same ratio distribution (mean 0.5) to the higher UTS, as instructed. If the real ratio were 0.44, as the Table 2 amplitudes would require at 991, the nominal runout amplitudes would equal case A's (the scatter would still differ; I did not run that case).

## Files
`table2_vs_uts.csv`, `pf_form_runout.csv`, `sobol_by_case.csv`, `sobol_top3.csv`, `life_v2prime_by_case.csv`, `pf_vs_load.csv`, `pf_vs_load_summary.csv`, `pf_vs_load_uts.png` (600 dpi; Goodman and Gerber at UTS 900 and 991, 800-1700 N, 1e5 samples per load, line at 1200 N), `run.log`.
