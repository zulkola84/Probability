"""Stage C: curved (circular-arc) plate section on rigid cylindrical rollers with frictionless contact.

Section: arc_section.py (net I at the hole = 105.72 mm^4, radial thickness 5.12, R_out ASSUMED).
Rollers: rigid cylinders, axis along y, radius ROLLER_R = 5 mm (ASSUMED), frictionless penalty contact
(node-to-surface, linear pressure-overclosure with stiffness CONTACT_K, ASSUMED numerical parameter).
Orientation ASSUMED: CONVEX_UP = False -> concave (bone) side up: supports touch the convex back at the apex,
load rollers touch the two inner edges of the concave face. Hole drilled along z (ASSUMED).
Quarter model (symmetry x = 0 and y = 0). Displacement control; total force = 4 x reaction of the load roller.
Usage: python stage_c.py   (see main)
"""
import re
import subprocess
import tempfile
from pathlib import Path

import gmsh
import numpy as np

import arc_section
import params as p
import run_ccx

ROLLER_R = 5.0          # mm, ASSUMED
CONTACT_K = 1.0e6       # MPa/mm, ASSUMED penalty stiffness
GAP = 0.01              # mm initial gap, closed by ADJUST
CONVEX_UP = False       # ASSUMED orientation
FACE_NODES = {1: (0, 1, 2), 2: (0, 3, 1), 3: (1, 3, 2), 4: (2, 3, 0)}   # CalculiX C3D10 faces (corner indices)
TOL = 1e-4


def geometry(r_out=None, convex_up=None, net=True):
    r_out = r_out or arc_section.R_OUT
    convex_up = CONVEX_UP if convex_up is None else convex_up
    h = p.THICKNESS
    alpha = arc_section.solve_alpha(r_out, net=net)
    return {"r_out": r_out, "r_in": r_out - h, "alpha": alpha, "s": 1 if convex_up else -1,
            "ye": r_out * np.sin(alpha), "yc": (r_out - h) * np.sin(alpha)}


def build_mesh(geo, hole_size, global_size=3.0, contact_size=0.6):
    """Mesh plate + two rollers (quarter model). Returns a dict with nodes, elements, sets and surfaces."""
    Lh = p.LENGTH / 2
    xl, xs = p.INNER_SPAN / 2, p.SUPPORT_SPAN / 2
    s, Ro, Ri, a = geo["s"], geo["r_out"], geo["r_in"], geo["alpha"]
    r = p.HOLE_DIAMETER / 2
    zc = 0.0
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("c")
    occ = gmsh.model.occ

    def disk(R):
        return occ.addDisk(0, 0, zc, R, R, zAxis=[1, 0, 0], xAxis=[0, 1, 0])
    ring, _ = occ.cut([(2, disk(Ro))], [(2, disk(Ri))])
    pts = [occ.addPoint(0, 0, zc), occ.addPoint(0, 0, zc + s * 2 * Ro),
           occ.addPoint(0, 2 * Ro * np.sin(a), zc + s * 2 * Ro * np.cos(a))]
    loop = occ.addCurveLoop([occ.addLine(pts[0], pts[1]), occ.addLine(pts[1], pts[2]), occ.addLine(pts[2], pts[0])])
    wedge = occ.addPlaneSurface([loop])
    band, _ = occ.intersect(ring, [(2, wedge)])
    out = occ.extrude(band, Lh, 0, 0)
    plate = [e for e in out if e[0] == 3]
    cyl = occ.addCylinder(0, 0, zc - Ro - 1, 0, 0, 2 * Ro + 2, r)
    plate, _ = occ.cut(plate, [(3, cyl)])

    # roller centres (z), tangent to the highest / lowest contact point plus GAP
    z_apex_top, z_apex_bot = (zc + Ro, zc + Ri * np.cos(a)) if s > 0 else (zc - Ri * np.cos(a), zc - Ro)
    zl, zs = z_apex_top + ROLLER_R + GAP, z_apex_bot - ROLLER_R - GAP
    Lr = geo["ye"] + 2.0
    roll_l = occ.addCylinder(xl, 0, zl, 0, Lr, 0, ROLLER_R)
    roll_s = occ.addCylinder(xs, 0, zs, 0, Lr, 0, ROLLER_R)
    occ.synchronize()

    plate_tags = [t for _, t in plate]
    hole_s = [t for _, t in gmsh.model.getEntitiesInBoundingBox(-TOL, -TOL, zc - Ro - 2, r + TOL, r + TOL, zc + Ro + 2, dim=2)]
    f = gmsh.model.mesh.field
    f.add("Distance", 1); f.setNumbers(1, "SurfacesList", hole_s); f.setNumber(1, "Sampling", 200)
    f.add("Threshold", 2); f.setNumber(2, "InField", 1); f.setNumber(2, "SizeMin", hole_size)
    f.setNumber(2, "SizeMax", global_size); f.setNumber(2, "DistMin", r * 0.3); f.setNumber(2, "DistMax", 4 * r)
    # refinement around the contact points: load rollers at the edge (concave up) or apex (convex up)
    edge = (geo["yc"] - 2.5, geo["ye"] + 2.5)
    apex = (-1.0, 3.0)
    load_y, supp_y = (apex, edge) if s > 0 else (edge, apex)
    ids = [2]
    for k, (x0, (y0, y1)) in enumerate(((xl, load_y), (xs, supp_y)), start=3):
        f.add("Box", k); f.setNumber(k, "VIn", contact_size); f.setNumber(k, "VOut", global_size)
        for nm, v in zip(("XMin", "XMax", "YMin", "YMax", "ZMin", "ZMax"), (x0 - 3, x0 + 3, y0, y1, -100, 100)):
            f.setNumber(k, nm, v)
        f.setNumber(k, "Thickness", 2.0)
        ids.append(k)
    f.add("Min", 9); f.setNumbers(9, "FieldsList", ids); f.setAsBackgroundMesh(9)
    for o in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
        gmsh.option.setNumber(o, 0)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")             # removes sliver tets that invert when mid-nodes go onto curved faces
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)
    gmsh.model.mesh.setOrder(2)

    ntags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = coords.reshape(-1, 3)
    index = {int(t): i for i, t in enumerate(ntags)}

    def tets(tag):
        et, _, en = gmsh.model.mesh.getElements(3, tag)
        assert et[0] == 11
        return en[0].reshape(-1, 10).astype(int)
    conn_plate = np.vstack([tets(t) for t in plate_tags])
    bodies = {"ROLL_L": tets(roll_l), "ROLL_S": tets(roll_s)}
    # roller boundary faces -> (element number, S-face) via the tet corner nodes
    faces = {}
    offset = len(conn_plate)
    surf = {}
    for name, conn, tag in (("ROLL_L", bodies["ROLL_L"], roll_l), ("ROLL_S", bodies["ROLL_S"], roll_s)):
        lookup = {}
        for e, c in enumerate(conn):
            for fi, idx in FACE_NODES.items():
                lookup[frozenset(c[list(idx)])] = (e, fi)
        bnd = gmsh.model.getBoundary([(3, tag)], oriented=False)
        lst = []
        for _, st in bnd:
            et, _, en = gmsh.model.mesh.getElements(2, abs(st))
            tri = en[0].reshape(-1, 6).astype(int)
            lst += [lookup[frozenset(t[:3])] for t in tri]
        surf[name] = (offset, lst)
        offset += len(conn)
    gmsh.finalize()
    X, Y, Z = xyz.T
    rad = np.hypot(Y, Z - zc)
    plate_nodes = np.unique(conn_plate)
    pn = np.zeros(len(ntags), bool); pn[[index[int(n)] for n in plate_nodes]] = True
    top_r, bot_r = (Ro, Ri) if s > 0 else (Ri, Ro)
    sl_top = pn & (np.abs(rad - top_r) < 1e-3) & (np.abs(X - xl) < 3.5)
    sl_bot = pn & (np.abs(rad - bot_r) < 1e-3) & (np.abs(X - xs) < 3.5)
    hole = pn & (np.abs(np.hypot(X, Y) - r) < 1e-3)
    ids = lambda m: [int(t) for t in ntags[m]]
    return {"ntags": ntags, "xyz": xyz, "plate": conn_plate, "bodies": bodies, "surf": surf,
            "SYMX": ids(pn & (np.abs(X) < TOL)), "SYMY": ids(pn & (np.abs(Y) < TOL)),
            "SL_TOP": ids(sl_top), "SL_BOT": ids(sl_bot), "HOLE": ids(hole),
            "refs": (int(ntags.max()) + 1, int(ntags.max()) + 2, int(ntags.max()) + 3, int(ntags.max()) + 4),
            "ref_xyz": [(xl, 0, zl), (xl, 0, zl), (xs, 0, zs), (xs, 0, zs)], "n_nodes": len(ntags)}


def write_deck(path, m, disp):
    order = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]
    nl, rl, ns, rs = m["refs"]
    with open(path, "w") as f:
        f.write("*NODE\n")
        for n, (x, y, z) in zip(m["ntags"], m["xyz"]):
            f.write(f"{int(n)}, {x:.8g}, {y:.8g}, {z:.8g}\n")
        for n, (x, y, z) in zip(m["refs"], m["ref_xyz"]):
            f.write(f"{n}, {x:.8g}, {y:.8g}, {z:.8g}\n")
        e0 = 1
        for name, conn in (("PLATE", m["plate"]), ("ROLL_L", m["bodies"]["ROLL_L"]), ("ROLL_S", m["bodies"]["ROLL_S"])):
            f.write(f"*ELEMENT, TYPE=C3D10, ELSET={name}\n")
            for i, c in enumerate(conn, e0):
                f.write(f"{i}, " + ", ".join(str(c[k]) for k in order) + "\n")
            e0 += len(conn)
        for k in ("SYMX", "SYMY", "SL_TOP", "SL_BOT", "HOLE"):
            f.write(f"*NSET, NSET={k}\n")
            ids = m[k]
            for i in range(0, len(ids), 8):
                f.write(", ".join(map(str, ids[i:i + 8])) + "\n")
        f.write(f"*NSET, NSET=REFL\n{nl}\n*NSET, NSET=REFS\n{ns}\n")
        for name, key in (("M_LOAD", "ROLL_L"), ("M_SUPP", "ROLL_S")):
            off, lst = m["surf"][key]
            f.write(f"*SURFACE, NAME={name}, TYPE=ELEMENT\n")
            for e, fi in lst:
                f.write(f"{off + e + 1}, S{fi}\n")
        f.write("*SURFACE, NAME=S_TOP, TYPE=NODE\nSL_TOP\n*SURFACE, NAME=S_BOT, TYPE=NODE\nSL_BOT\n")
        f.write(f"*MATERIAL, NAME=TI4\n*ELASTIC\n{p.E}, {p.NU}\n*MATERIAL, NAME=STEEL\n*ELASTIC\n210000., 0.3\n")
        f.write("*SOLID SECTION, ELSET=PLATE, MATERIAL=TI4\n*SOLID SECTION, ELSET=ROLL_L, MATERIAL=STEEL\n"
                "*SOLID SECTION, ELSET=ROLL_S, MATERIAL=STEEL\n")
        f.write(f"*RIGID BODY, ELSET=ROLL_L, REF NODE={nl}, ROT NODE={rl}\n"
                f"*RIGID BODY, ELSET=ROLL_S, REF NODE={ns}, ROT NODE={rs}\n")
        f.write("*SURFACE INTERACTION, NAME=SI\n*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR\n"
                f"{CONTACT_K:.4g}\n")
        f.write("*CONTACT PAIR, INTERACTION=SI, TYPE=NODE TO SURFACE, ADJUST=0.05\nS_TOP, M_LOAD\n"
                "*CONTACT PAIR, INTERACTION=SI, TYPE=NODE TO SURFACE, ADJUST=0.05\nS_BOT, M_SUPP\n")
        f.write(f"*BOUNDARY\nSYMX, 1, 1\nSYMY, 2, 2\n{nl}, 1, 2\n{nl}, 3, 3, {-disp}\n{rl}, 1, 3\n"
                f"{ns}, 1, 3\n{rs}, 1, 3\n")
        f.write("*STEP, INC=200\n*STATIC\n0.1, 1., 1e-5, 0.25\n*NODE FILE\nU, S\n"
                "*NODE PRINT, NSET=REFL\nRF\n*NODE PRINT, NSET=REFS\nRF\n*END STEP\n")


def read_rf(dat, nset):
    """Last-increment RF_z of the single node in nset."""
    val, on = None, False
    for line in open(dat):
        if "forces (fx,fy,fz) for set" in line:
            on = nset in line.upper().split("SET")[1].split()[0]
        elif on:
            t = line.split()
            if len(t) == 4:
                val = float(t[3])
    return val


def run(m, disp, workdir):
    deck = Path(workdir) / "c.inp"
    write_deck(deck, m, disp)
    r = subprocess.run(["ccx", "c"], cwd=workdir, capture_output=True, text=True)
    if r.returncode != 0 or "error" in r.stdout.lower().replace("*error", "ERROR")[:0]:
        raise RuntimeError(r.stdout[-1500:] + r.stderr[-500:])
    rf_l, rf_s = read_rf(Path(workdir) / "c.dat", "REFL"), read_rf(Path(workdir) / "c.dat", "REFS")
    stress = run_ccx.read_frd_stress(Path(workdir) / "c.frd")
    hole = [n for n in m["HOLE"] if n in stress]
    return {"disp": disp, "force": 4 * abs(rf_l), "force_support": 4 * abs(rf_s),
            "peak_vm": max(run_ccx.von_mises(stress[n]) for n in hole),
            "peak_p1": max(run_ccx.max_principal(stress[n]) for n in hole), "stdout_tail": r.stdout[-300:]}


if __name__ == "__main__":
    geo = geometry()
    print(geo)
    m = build_mesh(geo, 0.55)
    print("nodes", m["n_nodes"], "elements", len(m["plate"]), "slave", len(m["SL_TOP"]), len(m["SL_BOT"]))
    with tempfile.TemporaryDirectory() as tmp:
        print(run(m, 1.0, tmp))
