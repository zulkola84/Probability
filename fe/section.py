"""Non-rectangular plate sections calibrated to the paper's moment of inertia.

The paper's I = 105.72 mm^4 is read as the NET section through the hole (its text: the
critical section "at the location of a hole ... consists of two distinct polylines"; Table 1's
plate area of 51.21 mm^2 is also a net value). Each section shape is therefore calibrated so
that I_net (hole strip |y| < d/2 removed) equals p.I_BENDING, with thickness t at mid-width.

Shapes: "rect"    t x w rectangle
        "ellipse" semi-axes (w/2, t/2), a crude stand-in for the paper's curved outline
Loading and support lines on the ellipse follow the curved surface (a flat roller touches only
the apex, so this is an idealisation: nodes along the arc get the same z displacement).

Usage: python section.py   (prints the calibrated sizes, then runs FE for both shapes)
"""
import numpy as np
import gmsh
from scipy.integrate import quad
from scipy.optimize import brentq

import make_mesh
import params as p
import run_ccx

TOL = 1e-4


def _height(y, half_w, t, shape):
    if shape == "rect":
        return t if abs(y) <= half_w else 0.0
    return t * np.sqrt(max(0.0, 1 - (y / half_w) ** 2))


def inertia(width, t, d, shape, net=True):
    """Second moment about the width axis (z = thickness direction); net removes the hole strip."""
    hw = width / 2
    lo = d / 2 if net else 0.0
    val, _ = quad(lambda y: _height(y, hw, t, shape) ** 3 / 12, lo, hw)
    return 2 * val


def area(width, t, d, shape):
    hw = width / 2
    val, _ = quad(lambda y: _height(y, hw, t, shape), d / 2, hw)
    return 2 * val


def calibrate_width(shape, I_target=None, t=None, d=None):
    I_target = I_target or p.I_BENDING
    t, d = t or p.THICKNESS, d or p.HOLE_DIAMETER
    return brentq(lambda w: inertia(w, t, d, shape) - I_target, d + 0.5, 60.0)


def build_ellipse(hole_size, global_size, name="ell", out_dir=None):
    """Quarter model of an elliptical-section plate (centre z = t/2) with a hole through z."""
    out_dir = out_dir or p.FE_DIR
    Lh, t = p.LENGTH / 2, p.THICKNESS
    a = p.WIDTH / 2                      # semi-width (p.WIDTH is the full section width here)
    r = p.HOLE_DIAMETER / 2
    xs, xl = p.SUPPORT_SPAN / 2, p.INNER_SPAN / 2

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)
    occ = gmsh.model.occ

    # three conforming segments so nodes exist on the load (x=xl) and support (x=xs) lines
    disk = occ.addDisk(0, 0, t / 2, a, t / 2, zAxis=[1, 0, 0], xAxis=[0, 1, 0])
    vols, face = [], (2, disk)
    for x0, x1 in ((0, xl), (xl, xs), (xs, Lh)):
        out = occ.extrude([face], x1 - x0, 0, 0)
        vols.append(next(e for e in out if e[0] == 3))
        face = next(e for e in out if e[0] == 2 and e != face)
    cutter = occ.addBox(-1, -a - 1, -1, Lh + 2, a + 1, t + 2)         # removes y < 0
    cyl = occ.addCylinder(0, 0, -1, 0, 0, t + 2, r)
    res, _ = occ.cut(vols, [(3, cutter), (3, cyl)])
    occ.fragment(res, [])
    occ.synchronize()

    hole_s = [tg for _, tg in gmsh.model.getEntitiesInBoundingBox(
        -TOL, -TOL, -TOL, r + TOL, r + TOL, t + TOL, dim=2)]
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", hole_s)
    gmsh.model.mesh.field.setNumber(1, "Sampling", 200)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", hole_size)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", global_size)
    gmsh.model.mesh.field.setNumber(2, "DistMin", r * 0.3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 4 * r)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    for o in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
        gmsh.option.setNumber(o, 0)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(2)

    conn = []
    for _, tag in gmsh.model.getEntities(3):
        et, _, en = gmsh.model.mesh.getElements(3, tag)
        assert et[0] == 11, "expected 10-node tetrahedra"
        conn.append(en[0].reshape(-1, 10))
    conn = np.vstack(conn).astype(int)
    order = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]

    ntags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = coords.reshape(-1, 3)
    X, Y, Z = xyz.T
    on_surface = np.abs((Y / a) ** 2 + ((Z - t / 2) / (t / 2)) ** 2 - 1) < 1e-4
    ids = lambda m: [int(i) for i in ntags[m]]
    sets = {
        "SYMX": ids(np.abs(X) < TOL),
        "SYMY": ids(np.abs(Y) < TOL),
        "SUPPORT": ids(on_surface & (np.abs(X - xs) < TOL) & (Z <= t / 2 + TOL)),
        "LOAD": ids(on_surface & (np.abs(X - xl) < TOL) & (Z >= t / 2 - TOL)),
        "HOLE": ids(np.abs(np.hypot(X, Y) - r) < TOL),
    }
    path = make_mesh.write_inp(out_dir / f"{name}.inp", ntags, coords, conn, order, sets)
    gmsh.finalize()
    return path, len(ntags)


def fe_run(shape, width, hole_size=0.4, global_size=3.0):
    p.WIDTH = width
    if shape == "rect":
        make_mesh.build(hole_size, global_size, "sec")
    else:
        build_ellipse(hole_size, global_size, "sec")
    return run_ccx.run("sec")


if __name__ == "__main__":
    import csv

    from stiffness_check import delta_per_N

    t, d = p.THICKNESS, p.HOLE_DIAMETER
    print(f"target net I = {p.I_BENDING} mm^4 (t={t}, d={d})")
    fitted = p.WIDTH
    rows = []
    for shape in ("rect", "ellipse"):
        w = calibrate_width(shape)
        gross = inertia(w, t, d, shape, net=False)
        r = fe_run(shape, w)
        nominal = p.RUNOUT_FORCE / 2 * p.LEVER_ARM * t / 2 / p.I_BENDING
        rows.append({"shape": shape, "width_mm": round(w, 2), "net_I": round(inertia(w, t, d, shape), 2),
                     "gross_I": round(gross, 1), "net_area_mm2": round(area(w, t, d, shape), 1),
                     "beam_K_gross_I": round(1 / delta_per_N(gross)), "fe_K": round(float(r["stiffness"])),
                     "fe_peak_vm": round(float(r["peak_vm"])), "kt_vs_netI_nominal": round(float(r["peak_vm"]) / nominal, 3)})
        print(rows[-1])
    p.WIDTH = fitted
    with open(p.RESULTS_DIR / "section_comparison.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("paper: FE peak 857.8 MPa, analytical 898.1 (Kt 1.66; FE/nominal 1.55); experiment K 363-394 N/mm")
