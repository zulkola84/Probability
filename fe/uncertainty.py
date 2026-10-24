"""Uncertain inputs for the Monte Carlo sensitivity study.

All distributions are assumptions (the paper gives no tolerances or scatter), except
that the nominal values come from params.py. Edit here to change the study.
FE_INPUTS enter the FE model; FORCE scales the stress analytically (the model is linear).
"""
from scipy import stats

import params as p

FE_INPUTS = {
    # width is NOT uncertain: it is a calibrated model constant (params.WIDTH), not a measured
    # tolerance, so giving it a spread would make the Sobol ranking reflect that assumption.
    "thickness": stats.norm(p.THICKNESS, 0.05),                  # mm, machining tolerance
    "hole_d": stats.uniform(p.HOLE_DIAMETER - 0.1, 0.2),         # mm, drilling tolerance
    "E": stats.norm(p.E, 0.05 * p.E),                            # MPa, 5 % sd
    "nu": stats.uniform(0.30, 0.06),                             # grade 4 Ti range
    "lever_arm": stats.uniform(p.LEVER_ARM - 0.5, 1.0),          # mm, roller placement
}
FORCE = stats.norm(p.RUNOUT_FORCE, 0.02 * p.RUNOUT_FORCE)        # N, load-cell accuracy

FE_NAMES = list(FE_INPUTS)
ALL_NAMES = FE_NAMES + ["force"]


def to_physical(u):
    """Map unit-cube samples (n, len(FE_NAMES)) to physical FE inputs."""
    import numpy as np
    u = np.clip(u, 1e-6, 1 - 1e-6)
    return np.column_stack([FE_INPUTS[k].ppf(u[:, i]) for i, k in enumerate(FE_NAMES)])
