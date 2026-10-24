# Stage B - Independent validation, implanted two-hole plate

## What was run
- `fe/make_mesh.py` extended with optional `hole_spacing` (previous version: `fe/_backup/make_mesh_prestageB.py`). Holes at x = +-spacing/2 on the length axis; quarter model keeps the x=0 symmetry plane. `hole_spacing=None` is the original single hole.
- `fe/stage_b.py`: geometry identical to the calibrated intact model (width 16.365 mm NOT retuned, t 5.12, hole 7.8, 80/38 mm rollers, 1200 N) plus a second hole. Centre-to-centre spacing 20, 30, 40 mm (**ASSUMED**; not in the paper), plus the single-hole control. Each case at hole-edge mesh sizes 0.55 (DOE mesh) and 0.35 mm.
- Peak = max von Mises on the hole surface (same metric as the intact calibration). Data: `stageB_results.csv`.

## Key numbers (MPa; paper implanted FE 811.7, analytical 768.6)
| spacing (ASSUMED) | peak, 0.55 mm | peak, 0.35 mm | vs 811.7 (0.35 mm) | stiffness N/mm |
|---|---|---|---|---|
| single hole (control) | 858.5 | 858.8 | +5.8 % | 471 |
| 20 mm | 857.0 | 860.1 | +6.0 % | 440 |
| 30 mm | 855.7 | 859.5 | +5.9 % | 440 |
| 40 mm | 856.0 | 859.1 | +5.8 % | 440 |

- Control reproduces the existing calibrated value (858.45 MPa at 0.55 mm), so the mesh edit did not change the intact model.
- Mesh sensitivity of the peak: <= 0.4 % between the two hole sizes.

## Gate: FE peak within 10 % of 811.7 MPa
**PASS** for all three spacings (+5.4 to +6.0 %; allowed band 730.5-892.9 MPa).

## Limitations - read before citing this as validation
- The pass is weak. Both holes sit in the constant-moment span (between the 80 mm rollers), so the peak is essentially independent of spacing (spread ~1 MPa) and equals the intact-plate peak. The model does not represent the paper's implanted section (I = 165.89 mm^4); I is not used by the FE model. The paper's own implanted peak is 5.4 % below its intact peak (811.7 vs 857.8); this model reproduces no such reduction. The +5.8 % offset is that missing reduction.
- Compared with the paper's analytical value 768.6, the model is +11.5 to +11.9 %; no gate applies to it.
- Assumed that the paper's 811.7 MPa refers to the same 1200 N load; the inputs file does not say so.
- Stiffness drops 7 % with the second hole (471 -> 440 N/mm); 440 is closer to the intact-plate experiment (363-394 N/mm) than 471 but still 12-21 % above it. No implanted-plate stiffness reference was available.
- Spacing, and the absence of the implanted section shape, remain open; Stage C (curved section) is the place to address the section.
