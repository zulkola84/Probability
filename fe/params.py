"""Shared parameters for the Bologna et al. 2023 four-point bending model.

Every value tagged "assumed" is not given in the paper (see inputs/bologna_inputs.md).
Units: mm, N, MPa.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FE_DIR = ROOT / "fe"
RESULTS_DIR = ROOT / "results"

# Plate
THICKNESS = 5.12        # paper
HOLE_DIAMETER = 7.8     # paper
LENGTH = 180.0          # assumed (must exceed the support span of 156)
WIDTH = 16.365          # assumed; tuned so FE von Mises peak at 1200 N = 857.8 MPa (paper FE)

# Four-point bending. Paper: 80 mm working length between the loading (inner)
# rollers, 38 mm loading span from each inner roller to its support (lever arm).
# Check: M*c/I with I = 105.72 mm^4 (paper) and a 38 mm lever arm gives 552 MPa
# nominal, x Kt 1.66 = 916 MPa vs the paper's 898 MPa analytical peak.
# The reverse reading (80 mm lever arm) gives ~1160 MPa nominal and does not fit.
INNER_SPAN = 80.0       # between loading rollers (hole sits at the centre)
LEVER_ARM = 38.0        # inner roller to outer support
SUPPORT_SPAN = INNER_SPAN + 2 * LEVER_ARM
I_BENDING = 105.72      # paper, mm^4 (section is not a plain rectangle)

# Material (grade 4 titanium)
E = 105000.0            # paper, Table 1
NU = 0.33               # assumed
UTS = 900.0             # paper

# Loading
RUNOUT_FORCE = 1200.0   # paper, total force over both loading rollers
APPLIED_DISP = 1.0      # displacement imposed on the loading line (linear model)

# Paper reference values
PAPER_PEAK_ANALYTICAL = 898.1
PAPER_PEAK_FE = 857.8
