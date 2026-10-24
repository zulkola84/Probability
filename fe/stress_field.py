"""Stress field of the nominal plate at the 1200 N runout load (manuscript stress figure).

Reruns the nominal quarter model ONCE with a fine hole mesh, keeps the CalculiX files, and writes
to results/stress_field/:
  stress_field.inp / .frd / .dat   CalculiX input and results (load-line displacement 1 mm)
  quarter.vtu                      ParaView file, stresses and displacements scaled to 1200 N
                                   (use ParaView's Reflect filter twice for the full plate)
  fig_stress_contour.png           von Mises on the tension face: (a) full plate, (b) hole zoom
  fig_stress_paths.png             sxx (a) across the net section, (b) through the thickness at the hole edge
  fig_deflection.png               FE deflection of the plate edge vs Euler-Bernoulli beam theory
  summary.json                     peak stress, location, Kt checks, stiffness, mesh size

Nothing in results/ outside results/stress_field/ is touched. Nominal inputs come from params.py.
Usage: python stress_field.py [hole_size_mm]     (default 0.3; global size 3.0 mm as in the DOE)
       python stress_field.py --post             (re-plot from the existing files, no new solve)
"""
import json
import sys
from pathlib import Path

import numpy as np

import params as p

GLOBAL_SIZE = 3.0
NAME = "stress_field"
DPI = 600


# ------------------------------------------------------------------ reading
def read_inp(path):
    """Nodes {id: (x, y, z)} and C3D10 connectivity (CalculiX order) from the deck."""
    nodes, elems, mode = {}, [], None
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("*"):
            u = s.upper()
            mode = "node" if u.startswith("*NODE") and "PRINT" not in u and "FILE" not in u \
                else "elem" if u.startswith("*ELEMENT") else None
            continue
        v = [t for t in s.split(",") if t.strip()]
        if mode == "node":
            nodes[int(v[0])] = tuple(float(t) for t in v[1:4])
        elif mode == "elem":
            elems.append([int(t) for t in v[1:11]])
    return nodes, elems


def read_nset(path, name):
    ids, on = [], False
    for line in open(path):
        if line.startswith("*"):
            on = line.strip().upper() == f"*NSET, NSET={name}"
        elif on:
            ids += [int(t) for t in line.split(",") if t.strip()]
    return ids


def read_frd_block(path, key):
    """{node: [values]} from the LAST result block whose header contains key (DISP or STRESS)."""
    out, active = {}, False
    for line in open(path):
        if line.startswith(" -4"):
            active = key in line
            if active:
                out = {}
        elif active and line.startswith(" -3"):
            active = False
        elif active and line.startswith(" -1"):
            node = int(line[3:13])
            rest = line[13:].rstrip("\n")
            out[node] = [float(rest[12 * i:12 * i + 12]) for i in range(len(rest) // 12)]
    return out


def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = s.T
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * (sxy ** 2 + syz ** 2 + szx ** 2))


# ------------------------------------------------------------------ beam theory
def beam_deflection(x, force, inertia):
    """Euler-Bernoulli deflection (positive down) of the simply supported span under two loads force/2 at +-a/2."""
    L = p.SUPPORT_SPAN
    s = np.asarray(x) + L / 2

    def one(c):
        P = force / 2
        b = L - c
        return np.where(s <= c, P * b * s * (L**2 - b**2 - s**2), P * c * (L - s) * (L**2 - c**2 - (L - s) ** 2)) / (6 * p.E * inertia * L)
    return one(p.LEVER_ARM) + one(L - p.LEVER_ARM)


# ------------------------------------------------------------------ post-processing
def postprocess(xyz, s, u, hole_mask, out, meta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.linewidth": 0.8})
    t, W, r = p.THICKNESS, p.WIDTH, p.HOLE_DIAMETER / 2
    tol = 1e-5
    vm = von_mises(s)
    x, y, z = xyz.T

    # peaks on the hole surface: fatigue is checked on the TENSION face (largest positive sxx);
    # the overall von Mises peak can sit on the compression face, since bending is nearly antisymmetric
    hid = np.flatnonzero(hole_mask)
    on0 = hid[np.abs(z[hid]) < tol]
    onT = hid[np.abs(z[hid] - t) < tol]
    z_face = 0.0 if s[on0, 0].max() >= s[onT, 0].max() else t
    fid = on0 if z_face == 0.0 else onT
    k = fid[np.argmax(s[fid, 0])]
    kv = hid[np.argmax(vm[hid])]
    peak = dict(vm=float(vm[k]), sxx=float(s[k, 0]), x=float(x[k]), y=float(y[k]), z=float(z[k]))
    overall = dict(vm=float(vm[kv]), sxx=float(s[kv, 0]), z=float(z[kv]))

    # nominal stresses at the hole section
    M = p.RUNOUT_FORCE * p.LEVER_ARM / 2
    s_rect = 6 * M / ((W - 2 * r) * t**2)
    s_paper = M * (t / 2) / p.I_BENDING

    # ---------- (1) contour on the tension face, mirrored to the full plate
    f = np.abs(z - z_face) < tol
    tri = mtri.Triangulation(x[f], y[f])
    cx_, cy_ = x[f][tri.triangles].mean(1), y[f][tri.triangles].mean(1)
    mask = np.hypot(cx_, cy_) < r * 0.999
    vmax = np.ceil(vm[f].max() / 50) * 50
    levels = np.arange(0, vmax + 1, 50)
    fig, axs = plt.subplots(2, 1, figsize=(7.2, 4.3), gridspec_kw={"height_ratios": [1, 2.2], "hspace": 0.45})
    for ax, xl in ((axs[0], (-p.LENGTH / 2, p.LENGTH / 2)), (axs[1], (-14, 14))):
        for sx in (1, -1):
            for sy in (1, -1):
                tr = mtri.Triangulation(sx * x[f], sy * y[f], triangles=tri.triangles, mask=mask)
                cs = ax.tricontourf(tr, vm[f], levels=levels, cmap="turbo")
        ax.set_xlim(*xl)
        ax.set_ylim(-W / 2, W / 2)
        ax.set_aspect("equal")
        ax.add_patch(plt.Circle((0, 0), r, fill=False, lw=0.6, ec="k"))
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
    for xv, lab in ((p.INNER_SPAN / 2, "load line"), (p.SUPPORT_SPAN / 2, "support line")):
        for sgn in (1, -1):
            axs[0].axvline(sgn * xv, color="k", lw=0.6, ls="--")
        axs[0].text(xv, W / 2 + 1.0, lab, ha="center", fontsize=7)
    axs[0].add_patch(plt.Rectangle((-14, -W / 2), 28, W, fill=False, lw=0.8, ec="k"))
    axs[0].set_title("(a) von Mises stress on the tension face at 1200 N", fontsize=9, loc="left", pad=12)
    axs[1].set_title("(b) Critical hole", fontsize=9, loc="left")
    axs[1].plot([peak["x"]], [abs(peak["y"])], marker="v", ms=7, mfc="w", mec="k")
    axs[1].annotate(f"peak {peak['vm']:.0f} MPa", (peak["x"], abs(peak["y"])), xytext=(4.5, W / 2 - 1.2),
                    fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.6))
    cb = fig.colorbar(cs, ax=axs, shrink=0.9, pad=0.02, ticks=np.arange(0, vmax + 1, 100))
    cb.set_label("von Mises stress (MPa)")
    fig.savefig(out / "fig_stress_contour.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    # ---------- (2) stress paths
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 2.9))
    net = (np.abs(x) < tol) & (np.abs(z - z_face) < tol) & (y >= r - tol)
    o = np.argsort(y[net])
    a1.plot(y[net][o] - r, s[net][o, 0], "o-", ms=2.5, lw=1, color="C0", label="FE, $\\sigma_{xx}$")
    a1.axhline(s_rect, color="k", ls="--", lw=0.8, label=f"nominal, rectangle {s_rect:.0f} MPa")
    a1.axhline(s_paper, color="0.45", ls=":", lw=1, label=f"nominal, $I$ = {p.I_BENDING} mm$^4$: {s_paper:.0f} MPa")
    a1.set_xlabel("distance from hole edge (mm)")
    a1.set_ylabel("$\\sigma_{xx}$ (MPa)")
    a1.set_title("(a) Net section, tension face", fontsize=9, loc="left")
    a1.legend(fontsize=7, frameon=False)
    edge = (np.abs(x) < tol) & (np.abs(y - r) < 1e-3)
    o = np.argsort(z[edge])
    zz = np.linspace(0, t, 50)
    sign = 1 if z_face == 0 else -1
    nom_z = sign * M * (t / 2 - zz) / ((W - 2 * r) * t**3 / 12)
    a2.plot(s[edge][o, 0], z[edge][o], "o-", ms=2.5, lw=1, color="C0", label="FE at hole edge")
    a2.plot(nom_z, zz, "k--", lw=0.8, label="nominal (rectangle, net)")
    a2.axvline(0, color="0.7", lw=0.6)
    a2.set_xlabel("$\\sigma_{xx}$ (MPa)")
    a2.set_ylabel("z (mm)")
    a2.set_title("(b) Through the thickness at the hole edge", fontsize=9, loc="left")
    a2.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fig_stress_paths.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    # ---------- (3) deflection of the outer edge, bottom face
    e = (np.abs(y - W / 2) < tol) & (np.abs(z) < tol)
    o = np.argsort(x[e])
    xe, we = x[e][o], -u[e][o, 2]
    xb = np.linspace(-p.SUPPORT_SPAN / 2, p.SUPPORT_SPAN / 2, 200)
    I_gross = W * t**3 / 12
    fig, ax = plt.subplots(figsize=(4.6, 2.8))
    ax.plot(np.r_[-xe[::-1], xe], np.r_[we[::-1], we], "o", ms=2, color="C0", label="FE, plate edge (bottom face)")
    ax.plot(xb, beam_deflection(xb, p.RUNOUT_FORCE, I_gross), "k--", lw=0.9, label="beam theory, gross section")
    for sgn in (1, -1):
        ax.axvline(sgn * p.INNER_SPAN / 2, color="0.6", lw=0.6, ls=":")
        ax.axvline(sgn * p.SUPPORT_SPAN / 2, color="0.6", lw=0.6)
    ax.invert_yaxis()
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("deflection at 1200 N (mm)")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fig_deflection.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    w_load = float(np.interp(p.INNER_SPAN / 2, xe, we))
    summary = {**meta,
               "peak_von_mises_MPa": peak["vm"], "peak_sxx_MPa": peak["sxx"],
               "peak_location_mm": [peak["x"], peak["y"], peak["z"]],
               "overall_peak_von_mises_MPa": overall["vm"], "overall_peak_sxx_MPa": overall["sxx"],
               "overall_peak_z_mm": overall["z"],
               "tension_face_z_mm": z_face,
               "nominal_rect_net_MPa": s_rect, "nominal_paper_I_MPa": s_paper,
               "Kt_rect_net": peak["sxx"] / s_rect, "Kt_vs_paper_nominal": peak["vm"] / s_paper,
               "paper_peak_FE_MPa": p.PAPER_PEAK_FE, "diff_vs_paper_FE_pct": 100 * (overall["vm"] / p.PAPER_PEAK_FE - 1),
               "deflection_at_load_line_edge_mm": w_load,
               "beam_deflection_at_load_line_mm": float(beam_deflection(p.INNER_SPAN / 2, p.RUNOUT_FORCE, I_gross))}
    json.dump(summary, open(out / "summary.json", "w"), indent=2)
    return summary


def write_vtu(path, xyz, conn, s, u, vm):
    import meshio
    meshio.write(path, meshio.Mesh(points=xyz, cells=[("tetra10", conn)],
                                   point_data={"von_mises_MPa": vm, "S_MPa": s, "U_mm": u}))


def main():
    import run_ccx
    if "--post" not in sys.argv:
        import make_mesh
    post_only = "--post" in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    hs = float(args[0]) if args else 0.3
    out = p.RESULTS_DIR / "stress_field"
    out.mkdir(parents=True, exist_ok=True)
    if post_only:      # reuse the existing .inp/.frd/.dat (no meshing, no solve)
        n_nodes = sum(1 for _ in read_inp(out / f"{NAME}.inp")[0])
        force = 4 * abs(run_ccx.read_total_rf(out / f"{NAME}.dat"))
        res = {"stiffness": force / p.APPLIED_DISP, "peak_vm": None}
    else:
        inp, n_nodes = make_mesh.build(hs, GLOBAL_SIZE, NAME, out_dir=out)
        res = run_ccx.run(NAME, out_dir=out)
        force = res["stiffness"] * p.APPLIED_DISP
    sf = p.RUNOUT_FORCE / force

    nodes, elems = read_inp(out / f"{NAME}.inp")
    S = read_frd_block(out / f"{NAME}.frd", "STRESS")
    U = read_frd_block(out / f"{NAME}.frd", "DISP")
    ids = np.array(sorted(n for n in nodes if n in S))
    index = {n: i for i, n in enumerate(ids)}
    xyz = np.array([nodes[n] for n in ids])
    s = np.array([S[n][:6] for n in ids]) * sf
    u = np.array([U[n][:3] for n in ids]) * sf
    hole = np.zeros(len(ids), bool)
    hole[[index[n] for n in read_nset(out / f"{NAME}.inp", "HOLE") if n in index]] = True

    meta = {"hole_size_mm": hs, "global_size_mm": GLOBAL_SIZE, "nodes": int(n_nodes),
            "force_at_unit_disp_N": force, "scale_to_1200N": sf, "stiffness_N_per_mm": res["stiffness"],
            "peak_vm_run_ccx_MPa": res["peak_vm"], "width_mm": p.WIDTH, "E_MPa": p.E, "nu": p.NU}
    summary = postprocess(xyz, s, u, hole, out, meta)
    try:
        conn = np.array([[index[n] for n in e] for e in elems])
        write_vtu(out / "quarter.vtu", xyz, conn, s, u, von_mises(s))
    except Exception as exc:  # meshio missing or a node without results
        print("VTU export skipped:", exc)
    for k in ("peak_von_mises_MPa", "overall_peak_von_mises_MPa", "diff_vs_paper_FE_pct", "Kt_rect_net",
              "Kt_vs_paper_nominal", "stiffness_N_per_mm"):
        print(f"{k:28s} {summary[k]:.3f}")
    print("Peak location (x, y, z) mm:", [round(v, 3) for v in summary["peak_location_mm"]])
    print("Files written to", out)


if __name__ == "__main__":
    main()
