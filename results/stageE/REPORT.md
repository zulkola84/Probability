# Stage E - Life distribution and design chart

## What was run
`fe/stage_e.py` (one run, about 50 s). Surrogate: GP on the 416-run Stage A set (Q2 0.9676, CV RMSE 3.20 MPa). Files: `basquin_check.csv`, `life_at_1200N.csv`, `fe_thickness.csv`, `thickness_scaling_check.csv`, `design_chart.csv`, `design_chart.png`.

## (1) Probabilistic Basquin life at 1200 N

### Verification of Table 2
Bologna et al. 2023 (doi:10.3389/fmede.2023.1241312), checked against the article's HTML text, **not** the PDF: I could save the PDF but had no tool here to read its compressed pages (no pdftotext/poppler or Python PDF library; I did not install any). Confirmed from the text: A = 670.53 / 549.36 MPa and b = -0.088 / -0.059 (Goodman / Gerber); the stress columns 364.50 (both rows), 198.14 (Goodman) and 241.84 (Gerber); constants fixed by a straight line through N = 10^3 (oligocyclic limit, 90 % of the static failure value) and N = 10^6 (runout); average predicted life about 575 thousand cycles (analytical) and about one million (FE). Computing A*N^b gives 365.1 / 198.8 MPa (Goodman) and 365.5 / 243.1 MPa (Gerber) at 10^3 / 10^6 cycles, consistent with the table within rounding of A and b.
**Not confirmed:** Basquin's Eq. 14, the headings of the stress columns, and which stress the equation takes (the text says "the nominal alternating stress ... taking into account the R-ratio", with Goodman/Gerber for mean stress). Please check these against the PDF.

### Reading of the equation used (ASSUMED) and its consistency check
N = (sigma_a / A)^(1/b), sigma_a = 0.45 x peak / 1.66. Applied to the paper's own peak stresses (`basquin_check.csv`):
| peak source | Goodman N | Gerber N | arithmetic mean | geometric mean | paper |
|---|---|---|---|---|---|
| analytical 898.1 MPa | 1.00e5 | 9.78e5 | 5.39e5 (-6 %) | 3.13e5 | ~575,000 |
| FE 857.8 MPa | 1.68e5 | 2.13e6 | 1.15e6 (+15 %) | 5.99e5 | ~1e6 |

The reading reproduces both of the paper's statements within 15 % if "average of the predicted range" means the arithmetic mean of the two criteria. This is supporting evidence only; the paper's definition of the average is unknown. One unresolved mismatch: with this reading the Gerber runout allowable at 1200 N is 241.5 MPa (Table 2: 241.84) but the Goodman one is 181.5 MPa against 198.14, so the Goodman row is not reproduced (-8 %).

### Result at 1200 N (400,000 MC samples; `life_at_1200N.csv`)
Cycles to failure, 5th / median / 95th percentile, and P(N < 1e6):
| criterion | variant | 5 % | median | 95 % | P(N < 1e6) |
|---|---|---|---|---|---|
| Goodman | FE inputs + force only | 1.00e5 | 1.71e5 | 2.93e5 | 1.000 |
| Goodman | + material scatter (UTS CV 3 %, ratio SD 0.05, ASSUMED) | 1.94e4 | 1.69e5 | 1.14e6 | 0.937 |
| Goodman | + life scatter log10 SD 0.2 (ASSUMED) | 1.72e4 | 1.68e5 | 1.32e6 | 0.921 |
| Gerber | FE inputs + force only | 9.78e5 | 2.18e6 | 4.87e6 | 0.055 |
| Gerber | + material scatter | 8.48e4 | 2.14e6 | 3.69e7 | 0.342 |
| Gerber | + life scatter | 7.85e4 | 2.13e6 | 4.09e7 | 0.347 |

- Against the paper: our medians (Goodman 1.7e5, Gerber 2.1e6) straddle the paper's ~575,000 (analytical) and ~1e6 (FE); the paper's figures sit between the two criteria, as the consistency check shows.
- Material scatter widens the 5-95 % interval by one to two orders of magnitude because the Basquin exponents are shallow (1/b = -11 and -17): a 1 % stress change moves life by about 11 % and 17 %.
- P(N < 1e6) under Gerber (0.342) agrees with P_f from the limit state in Stage D (0.337). Under Goodman it is 0.937 against 0.980 because of the Goodman endpoint mismatch above.
- Life scatter about the S-N line (log10 SD 0.2) has **no source**: no published grade 4 value was found. It is shown only as a sensitivity.

## (2) Design chart: thickness for target beta (`design_chart.csv`, `design_chart.png`)
Method: limit state of Stage D (Goodman or Gerber, Kf = Kt, assumed scatter), beta = -Phi^-1(P_f), 2,000,000 MC samples with common random numbers, thickness solved by root finding. Peak stress at thickness t = GP(x) x R(t) x load/1200, with R(t) from 12 FE runs at nominal inputs (t = 4-12 mm, `fe_thickness.csv`; nominal peak 1367.9 MPa at 4 mm, 858.5 at 5.12, 479.1 at 7, 178.5 at 12). Load scatter (CV 2 %) is kept about the design load; width 16.365 mm and hole 7.8 mm are unchanged.

Required thickness (mm):
| load (N) | Goodman b2 | b3 | b4 | Gerber b2 | b3 | b4 |
|---|---|---|---|---|---|---|
| 800 | 4.96 | 5.26 | 5.61 | 4.50 | 4.80 | 5.16 |
| 900 | 5.28 | 5.60 | 5.98 | 4.79 | 5.11 | 5.50 |
| 1000 | 5.59 | 5.92 | 6.33 | 5.07 | 5.41 | 5.82 |
| 1100 | 5.88 | 6.23 | 6.66 | 5.33 | 5.69 | 6.12 |
| 1200 | 6.16 | 6.53 | 6.98 | 5.59 | 5.96 | 6.41 |
| 1300 | 6.43 | 6.82 | 7.29 | 5.83 | 6.22 | 6.69 |
| 1400 | 6.69 | 7.10 | 7.59 | 6.07 | 6.47 | 6.97 |

Checks: re-computed P_f at each solved thickness matches the target within 0.6 %; at 5.12 mm and 1200 N the chart's P_f is 0.980 (Goodman) and 0.336 (Gerber), as in Stage D. The R(t) x GP scaling was tested against fresh FE at random inputs with thickness shifted to 6.0 and 7.5 mm (24 points each): bias -0.04 / -0.16 %, RMSE 0.43 / 0.64 %, max |error| 0.78 / 1.57 % (`thickness_scaling_check.csv`).

## Limitations
- Required thicknesses depend on the assumed scatter (endurance ratio SD, UTS CV, force CV), on Kf = Kt, and on the spread between Goodman and Gerber, which is larger (0.46-0.62 mm at equal beta and load) than the change between adjacent beta levels (0.3-0.5 mm). Treat the two columns as a bracket. Stage D showed Kf choice can change P_f far more; not applied here.
- Thickness scaling assumes the same width, hole, plane section and calibrated model (Stage C's stiffness mismatch applies). Thicknesses above about 7.6 mm are outside the range of the scaling validation (up to 7.5 mm); the FE sweep itself goes to 12 mm. Thick plates may not be clinically feasible; stiffness and bone-stress effects of thickness are not considered.
- beta = 4 corresponds to P_f = 3.2e-5, about 63 failures in 2e6 samples (relative MC error about 13 %, a thickness effect of about 0.01 mm), but the tail depends on the Gaussian and lognormal input assumptions.
- A, b and the 0.45 R-factor reading of the equation are partly ASSUMED (see above); the comparison with 575,000 and 1e6 cycles is a consistency check, not a validation, since the paper's averaging is undefined.
