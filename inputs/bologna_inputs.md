# Bologna et al. 2023 inputs

Source: Bologna et al. 2023, Frontiers in Medical Engineering,
https://www.frontiersin.org/journals/medical-engineering/articles/10.3389/fmede.2023.1241312/full

Values were read via a summarizing fetch; spot-check against the PDF.

| Item | Value | Status |
|---|---|---|
| **Plate** | | |
| Length | – | assumed (not stated) |
| Width | – | assumed (not stated) |
| Thickness | 5.12 mm | paper |
| Hole diameter | 7.8 mm | paper |
| Number of holes | – | assumed (not stated; implanted model has two holes) |
| Hole spacing | – | assumed (not stated) |
| Countersink | – | assumed (not stated) |
| **Test set-up** | | |
| Loading span | 38 mm | paper |
| Other span | 80 mm | paper (see note 2) |
| Roller diameter | – | assumed (not stated; rollers modeled as rigid cylinders) |
| Load ratio R | 0.1 | paper |
| **Material** | | |
| Material | Medical grade 4 titanium | paper (ISO 5832-2) |
| Tensile strength | 900 MPa | paper (raw material certificate) |
| Young's modulus | 105 GPa | paper (Table 1; calibrated against static tests) |
| Poisson's ratio | – | assumed (not stated; typical grade 4 Ti about 0.32–0.34) |
| **Results** | | |
| Runout load | 1,200 N | paper (1 million cycles) |
| Analytical peak stress | 898.1 MPa | paper |
| FE peak stress | 857.8 MPa | paper |
| Peak stress difference | about 4.5–4.7% | calculated |
| Stress concentration factor | 1.66 analytical, 1.59 FE | paper |
| Implanted plate peak stress | 768.6 MPa analytical, 811.7 MPa FE (5.3%) | paper |

## Notes

1. The paper's 4.1% is the difference between the stress concentration factors, not the peak stresses. The peak stresses differ by about 4.7% relative to FE, or 4.5% relative to analytical.
2. Spans: 80 mm is the working length between the two loading (inner) rollers, with the hole at its centre; 38 mm is the distance from each inner roller to its outer support (lever arm), so the outer span is 156 mm. Supporting check: the paper's I = 105.72 mm^4 (four-point bending section) gives nominal stress 600 N x 38 mm x 2.56 mm / 105.72 = 552 MPa, x Kt 1.66 = 916 MPa vs the paper's 898 MPa. The 80 mm lever arm does not fit.
3. Plate cross-section is not a plain rectangle (the paper mentions curvature and I = 105.72 mm^4 for bending, 165.89 mm^4 implanted). Plate width is not stated; the FE model uses an assumed 16 mm rectangle tuned to the paper's peak stress.
4. Source of the experiments, Bologna, Terzini & Audenino 2022, is a conference paper (ESB Porto, https://esbiomech.org/conference/archive/2022porto/papers/porto_1919.pdf), not retrieved.
5. Figure 6 digitised from the published image (axes calibrated from tick marks). Only three red strands are separable (1-3 mm): the thick upper strand (probably two overlapping curves) and two thin lower strands. Linear-fit stiffness, 1.0-3.0 mm: upper about 394 N/mm, middle about 363, lower about 363; FE (blue) about 372. Files: fig6_strand_{upper,mid,lower}.csv, fig6_fe_blue.csv, fig6_red_{upper,lower}.csv.
6. Plate geometry from Bologna, Audenino & Terzini 2023, Ann. Biomed. Eng., "Bone Plates Runout Prediction Through Tensile Strength and Geometric Properties for Regulatory Mechanical Testing" (DOI 10.1007/s10439-023-03363-2, PMC10810033), Table 2 (image saved as annals2023_table2_cross_sections.png). Plate 11 (pure titanium, diaphyseal femur, sigma_UTS 991 MPa) is probably the same plate family as the Frontiers one but is not confirmed to be the same plate: t = 5.2 mm, w = 16.7 mm, d = 6.5 mm (minimum hole diameter extracted from CAD), Ixx = 109.0 mm^4. Frontiers values for comparison: t = 5.12, d = 7.8, I = 105.72, sigma_UTS 900. Width is the bounding box of the section, which is non-rectangular with a curved, stepped outline. Plate length and hole count are still not given.
