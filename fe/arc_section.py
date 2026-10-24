"""Circular-arc (curved) plate section for Stage C, calibrated to the paper's moment of inertia.

Section: concentric arc band, radial thickness h (= paper thickness 5.12 mm), outer radius R_OUT (ASSUMED),
bounded by two radial cuts at +-alpha. The half-angle alpha is solved so that the NET section at the hole
(strip |y| < d/2 removed) has I = p.I_BENDING about its own centroidal axis (same reading as section.py).
Coordinates here: y across the width, z up, apex of the band at y = 0.
"""
import numpy as np
from scipy.optimize import brentq

import params as p

R_OUT = 25.0          # mm, ASSUMED outer radius of curvature (not in the paper)
STEP = 0.02           # mm, grid step for the numerical section integrals


def _band_mask(y, z, r_out, h, alpha):
    r = np.hypot(y, z)
    ang = np.arctan2(np.abs(y), z)           # angle from the apex direction
    return (r >= r_out - h) & (r <= r_out) & (ang <= alpha)


def inertia(r_out, h, alpha, d, net=True):
    """Second moment about the centroidal horizontal axis (bending in z)."""
    ymax = r_out * np.sin(alpha)
    y = np.arange(d / 2 if net else 0.0, ymax, STEP) + STEP / 2
    z = np.arange((r_out - h) * np.cos(alpha) * 0.99, r_out, STEP) + STEP / 2
    Y, Z = np.meshgrid(y, z)
    m = _band_mask(Y, Z, r_out, h, alpha)
    A = m.sum() * STEP**2
    zc = (Z * m).sum() * STEP**2 / A
    I = ((Z - zc) ** 2 * m).sum() * STEP**2
    return 2 * I, 2 * A, zc          # both half-widths


def solve_alpha(r_out=R_OUT, h=None, d=None, I_target=None, net=True):
    h, d, I_target = h or p.THICKNESS, d or p.HOLE_DIAMETER, I_target or p.I_BENDING
    f = lambda a: inertia(r_out, h, a, d, net)[0] - I_target
    return brentq(f, np.arcsin((d / 2 + 0.5) / r_out), 1.45, xtol=1e-5)


if __name__ == "__main__":
    from stiffness_check import delta_per_N
    h, d = p.THICKNESS, p.HOLE_DIAMETER
    for R in (15, 20, 25, 35, 60, 200):
        a = solve_alpha(R)
        Ig = inertia(R, h, a, d, net=False)[0]
        print(f"R_out={R:>4}: alpha={np.degrees(a):5.1f} deg  bbox width={2 * R * np.sin(a):5.2f}  "
              f"gross I={Ig:6.1f}  beam K(gross I)={1 / delta_per_N(Ig):.0f} N/mm")
