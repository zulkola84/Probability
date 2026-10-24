"""Build the quarter-plate mesh with Gmsh and write a CalculiX input deck.

Quarter model: symmetry planes x=0 (mid-length) and y=0 (mid-width), hole
centred on the origin, thickness along z. The rollers are idealised as lines:
the support line (bottom face, x=SUPPORT_SPAN/2) is fixed in z and the loading
line (top face, x=INNER_SPAN/2) is displaced in -z. Roller diameter is unused.

Optional hole_spacing (centre-to-centre, mm) gives the two-hole plate of Stage B: holes at x = +-spacing/2
on the length axis; the quarter model keeps the x=0 symmetry plane and holds the hole at +spacing/2.
hole_spacing=None is the original single hole at the origin.

Usage: python make_mesh.py [hole_size] [global_size] [name]
"""
import sys

import gmsh

import params as p

TOL = 1e-6


def _nodes(dim, tags):
    out = set()
    for tag in tags:
        ids, _, _ = gmsh.model.mesh.getNodes(dim, tag, includeBoundary=True)
        out.update(int(i) for i in ids)
    return sorted(out)


def _entities(dim, box):
    return [t for _, t in gmsh.model.getEntitiesInBoundingBox(*box, dim=dim)]


def write_inp(path, ntags, coords, conn, order, sets):
    """Write the CalculiX deck: nodes, 10-node tets, node sets, material, BCs and the static step."""
    with open(path, "w") as f:
        f.write("*NODE\n")
        for n, (x, y, z) in zip(ntags, coords.reshape(-1, 3)):
            f.write(f"{int(n)}, {x:.8g}, {y:.8g}, {z:.8g}\n")
        f.write("*ELEMENT, TYPE=C3D10, ELSET=PLATE\n")
        for i, c in enumerate(conn, 1):
            f.write(f"{i}, " + ", ".join(str(c[k]) for k in order) + "\n")
        for k, ids in sets.items():
            f.write(f"*NSET, NSET={k}\n")
            for i in range(0, len(ids), 8):
                f.write(", ".join(map(str, ids[i:i + 8])) + "\n")
        f.write(f"*MATERIAL, NAME=TI4\n*ELASTIC\n{p.E}, {p.NU}\n")
        f.write("*SOLID SECTION, ELSET=PLATE, MATERIAL=TI4\n")
        f.write("*BOUNDARY\nSYMX, 1, 1\nSYMY, 2, 2\nSUPPORT, 3, 3\n")
        f.write(f"LOAD, 3, 3, {-p.APPLIED_DISP}\n")
        f.write("*STEP\n*STATIC\n")
        f.write("*NODE FILE\nU, RF, S\n")
        f.write("*NODE PRINT, NSET=LOAD, TOTALS=ONLY\nRF\n")
        f.write("*END STEP\n")
    return path


def build(hole_size, global_size, name="plate", out_dir=None, hole_spacing=None):
    """Mesh the plate and write <name>.inp; returns the .inp path and node count."""
    out_dir = out_dir or p.FE_DIR
    L, W, t = p.LENGTH / 2, p.WIDTH / 2, p.THICKNESS
    r = p.HOLE_DIAMETER / 2
    cx = 0.0 if hole_spacing is None else hole_spacing / 2
    xs, xl = p.SUPPORT_SPAN / 2, p.INNER_SPAN / 2

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)
    occ = gmsh.model.occ

    box = occ.addBox(0, 0, 0, L, W, t)
    cyl = occ.addCylinder(cx, 0, -1, 0, 0, t + 2, r)
    plate, _ = occ.cut([(3, box)], [(3, cyl)])
    # lines on the faces that become the load / support node sets
    top = occ.addLine(occ.addPoint(xl, 0, t), occ.addPoint(xl, W, t))
    bot = occ.addLine(occ.addPoint(xs, 0, 0), occ.addPoint(xs, W, 0))
    occ.fragment(plate, [(1, top), (1, bot)])
    occ.synchronize()

    vol = [t_ for _, t_ in gmsh.model.getEntities(3)]
    load_c = _entities(1, (xl - TOL, -TOL, t - TOL, xl + TOL, W + TOL, t + TOL))
    supp_c = _entities(1, (xs - TOL, -TOL, -TOL, xs + TOL, W + TOL, TOL))
    sym_x = _entities(2, (-TOL, -TOL, -TOL, TOL, W + TOL, t + TOL))
    sym_y = _entities(2, (-TOL, -TOL, -TOL, L + TOL, TOL, t + TOL))
    hole_s = _entities(2, (cx - r - TOL, -TOL, -TOL, cx + r + TOL, r + TOL, t + TOL))

    # refine toward the hole edge, grade out to global_size
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
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(2)

    # elements (10-node tets); gmsh and CalculiX differ in the last two mid-edge nodes
    etypes, etags, enodes = gmsh.model.mesh.getElements(3, vol[0])
    assert etypes[0] == 11, "expected 10-node tetrahedra"
    conn = enodes[0].reshape(-1, 10).astype(int)
    order = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]

    ntags, coords, _ = gmsh.model.mesh.getNodes()
    sets = {
        "SYMX": _nodes(2, sym_x),
        "SYMY": _nodes(2, sym_y),
        "SUPPORT": _nodes(1, supp_c),
        "LOAD": _nodes(1, load_c),
        "HOLE": _nodes(2, hole_s),
    }
    n_nodes = len(ntags)

    path = write_inp(out_dir / f"{name}.inp", ntags, coords, conn, order, sets)

    gmsh.finalize()
    return path, n_nodes


if __name__ == "__main__":
    hs = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    gs = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
    nm = sys.argv[3] if len(sys.argv) > 3 else "plate"
    print(build(hs, gs, nm))
