"""Phase 2, step 4: elastic-plastic check at the hole edge of the calibrated quarter model.

Mesh: make_mesh.build (hole size 0.3 mm, global 3.0 mm, W = 16.365 mm, C3D10), boundary conditions as make_mesh.write_inp
(SYMX, SYMY, SUPPORT, displacement of the LOAD line). Material: E, nu from params; *PLASTIC, HARDENING=ISOTROPIC, bilinear from
(yield, 0) to (900 MPa, plastic strain 0.10). The three YIELD strengths 650 / 750 / 850 MPa are ASSUMED (not published).
Displacement control; total force = 4 x |RF_z of LOAD| (the quarter model holds 1/4 of the total force over both rollers,
same scaling as run_ccx.run / stress_field.py). A displacement search (at most 4 solves per case, tolerance 1 % of each
target: 1200 N and 120 N) sets the peak and trough displacements; the final solve is 3 cycles (load leg + unload leg each),
10 increments per leg. The unloading displacement for 120 N is found on the unloading branch of the same solve.
Critical node: node of the HOLE set on the tension face (z = 0, as stress_field.py) with the largest sxx in the elastic run.
Output is restricted to the node set OUTN (hole nodes on z = 0 or x = 0, and the net-section path x = 0, z = 0, y >= r).
Plastic zone: nodes with PEEQ > PEEQ_TOL (ASSUMED numerical threshold) - resolution is the nodal spacing.
All files go to OUT (outside results/).
"""
import json
import os
import re
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fatigue as ft
import make_mesh
import params as p
import sensitivity as sens
import stress_field as sf
import uncertainty as unc

OUT = Path(__file__).resolve().parents[1] / "results" / "phase2" / "step4_plasticity"
YIELDS = (650, 750, 850)            # ASSUMED
UTS_END = 900.0                     # MPa, end of the bilinear curve (the paper's UTS), plastic strain 0.10 (ASSUMED)
EPS_END = 0.10                      # ASSUMED
PEEQ_TOL = 1e-5                     # ASSUMED numerical threshold for "plastic"
NINC = 6                            # increments per leg
LOAD_MAX, LOAD_MIN = 1200.0, 120.0  # N
HOLE_SIZE, GLOBAL_SIZE = 0.3, 3.0
CCX_ENV = {**os.environ, "OMP_NUM_THREADS": "2", "CCX_NPROC_EQUATION_SOLVER": "2"}
CYCLES = 3


def label(y):
    return "elastic" if y is None else f"yield{y}_ASSUMED"


# ---------------------------------------------------------------- deck
def head_and_sets():
    path, n = make_mesh.build(HOLE_SIZE, GLOBAL_SIZE, "mesh_base", out_dir=OUT)
    txt = open(path).read()
    head = txt[:txt.index("*MATERIAL")]
    nodes, _ = sf.read_inp(path)
    hole = sf.read_nset(path, "HOLE")
    r, tol = p.HOLE_DIAMETER / 2, 1e-5
    outn = set(n_ for n_ in hole if abs(nodes[n_][2]) < tol or abs(nodes[n_][0]) < tol)
    outn |= {n_ for n_, c in nodes.items() if abs(c[0]) < tol and abs(c[2]) < tol and c[1] >= r - tol}
    return head, sorted(outn), nodes, hole, n


def write_deck(name, head, outn, yld, legs, nlgeom=False, fname=None):
    with open(OUT / (fname or f"{name}.inp"), "w") as f:
        f.write(head)
        f.write("*NSET, NSET=OUTN\n")
        for i in range(0, len(outn), 8):
            f.write(", ".join(map(str, outn[i:i + 8])) + "\n")
        f.write(f"*MATERIAL, NAME=TI4\n*ELASTIC\n{p.E}, {p.NU}\n")
        if yld is not None:
            f.write(f"*PLASTIC, HARDENING=ISOTROPIC\n{yld}., 0.\n{UTS_END}, {EPS_END}\n")
        f.write("*SOLID SECTION, ELSET=PLATE, MATERIAL=TI4\n")
        f.write("*BOUNDARY\nSYMX, 1, 1\nSYMY, 2, 2\nSUPPORT, 3, 3\n")
        for u in legs:
            dt = 1.0 / NINC
            f.write("*STEP, INC=2000" + (", NLGEOM" if nlgeom else "") + "\n*STATIC\n")
            f.write(f"{dt}, 1., 1e-6, {dt}\n")
            f.write(f"*BOUNDARY\nLOAD, 3, 3, {-u}\n")
            f.write("*NODE FILE, NSET=OUTN\nS, E, PEEQ\n")
            f.write("*NODE PRINT, NSET=LOAD, TOTALS=ONLY\nRF\n")
            f.write("*END STEP\n")


# ---------------------------------------------------------------- parsing
def parse_dat(path):
    out, t = [], None
    for line in open(path):
        m = re.search(r"total force \(fx,fy,fz\) for set LOAD and time\s+([-+\d.Ee]+)", line)
        if m:
            t = float(m.group(1)); continue
        if t is not None:
            v = line.split()
            if len(v) == 3:
                try:
                    out.append((t, 4 * abs(float(v[2])))); t = None
                except ValueError:
                    pass
    return out


def parse_frd(path):
    """{time: {"S": {node: [6]}, "E": {...}, "PE": {node: value}}}"""
    res, t, key = {}, None, None
    for line in open(path):
        if line.startswith("  100CL"):
            t = float(line.split()[2]); res.setdefault(t, {"S": {}, "E": {}, "PE": {}}); key = None
        elif line.startswith(" -4"):
            nm = line.split()[1]
            key = {"STRESS": "S", "TOSTRAIN": "E", "PE": "PE"}.get(nm)
        elif line.startswith(" -3"):
            key = None
        elif key and line.startswith(" -1"):
            rest = line[13:].rstrip("\n")
            v = [float(rest[12 * i:12 * i + 12]) for i in range(len(rest) // 12)]
            res[t][key][int(line[3:13])] = v if key != "PE" else v[0]
    return res


def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = s
    return float(np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * (sxy**2 + syz**2 + szx**2)))


# ---------------------------------------------------------------- one solve
def solve(name, head, outn, yld, u_max, u_min, nlgeom=False):
    legs = [u_max, u_min] * CYCLES
    inp = OUT / f"{name}.inp"
    write_deck(name, head, outn, yld, legs, nlgeom, fname=f"{name}.tmp_inp")
    same = inp.exists() and (OUT / f"{name}.frd").exists() and (OUT / f"{name}.stdout.txt").exists() \
        and inp.read_text() == (OUT / f"{name}.tmp_inp").read_text()
    if same:      # identical deck already solved (re-analysis after a post-processing error): reuse; wall time from file times
        (OUT / f"{name}.tmp_inp").unlink()
        wall = (OUT / f"{name}.frd").stat().st_mtime - inp.stat().st_mtime
        log = (OUT / f"{name}.stdout.txt").read_text()
        pr = type("P", (), {"returncode": 0})()
    else:
        (OUT / f"{name}.tmp_inp").replace(inp)
        t0 = time.perf_counter()
        pr = subprocess.run(["ccx", name], cwd=OUT, capture_output=True, text=True, env=CCX_ENV)
        wall = time.perf_counter() - t0
        log = pr.stdout + pr.stderr
        (OUT / f"{name}.stdout.txt").write_text(log)
    warns = sorted({l.strip() for l in log.splitlines() if re.search(r"warning|error|not converge|cutback|cut back|divergen", l, re.I)})
    if pr.returncode != 0 or "*ERROR" in log:
        raise RuntimeError(f"{name}: ccx failed\n{log[-1500:]}")
    dat = parse_dat(OUT / f"{name}.dat")
    frd = parse_frd(OUT / f"{name}.frd")
    n_inc = sum(1 for l in open(OUT / f"{name}.sta") if re.match(r"\s+\d+\s+\d+\s+\d+\s+\d+", l))
    return {"name": name, "wall_s": wall, "dat": dat, "frd": frd, "warnings": warns, "n_increments": n_inc,
            "legs": legs}


def leg_curve(res, leg_no):
    """(u, F) along leg leg_no (1-based) from the printed totals; u from the linear ramp inside the step."""
    legs = res["legs"]
    start = 0.0 if leg_no == 1 else legs[leg_no - 2]
    end = legs[leg_no - 1]
    pts = [(start + (t - (leg_no - 1)) * (end - start), F) for t, F in res["dat"] if leg_no - 1 - 1e-9 < t <= leg_no + 1e-9]
    return np.array(pts)


def force_at(res, t_end):
    return [F for t, F in res["dat"] if abs(t - t_end) < 1e-6][0]


def search(yld, head, outn, u_el, k_el):
    """Displacement search (<= 4 solves). Returns the final solve result and the search log."""
    u_max = 1.04 * u_el                         # first guess; plasticity can only lower the force at a given displacement
    u_min = u_max - (LOAD_MAX - LOAD_MIN) / k_el
    log = []
    for it in range(1, 5):
        res = solve(f"{label(yld)}_s{it}", head, outn, yld, u_max, u_min)
        fmax, fmin = force_at(res, 1.0), force_at(res, 2.0)
        ok = abs(fmax - LOAD_MAX) <= 0.01 * LOAD_MAX and abs(fmin - LOAD_MIN) <= 0.01 * LOAD_MIN
        log.append({"case": label(yld), "solve": it, "u_max": u_max, "u_min": u_min, "F_max_N": fmax, "F_min_N": fmin,
                    "converged_search": ok, "wall_s": res["wall_s"], "n_increments": res["n_increments"]})
        if ok or it == 4:
            break
        c1, c2 = leg_curve(res, 1), leg_curve(res, 2)
        if c1[:, 1].min() <= LOAD_MAX <= c1[:, 1].max():
            um = float(np.interp(LOAD_MAX, c1[:, 1], c1[:, 0]))
        else:
            um = u_max * LOAD_MAX / fmax
        c2s = c2[np.argsort(c2[:, 1])]
        umin = float(np.interp(LOAD_MIN, c2s[:, 1], c2s[:, 0])) if c2s[0, 1] <= LOAD_MIN <= c2s[-1, 1] \
            else u_min + (LOAD_MIN - fmin) / k_el
        u_min = umin + (um - u_max)             # unloading branch shifts with the peak displacement
        u_max = um
    return res, log


# ---------------------------------------------------------------- analysis of one finished case
def analyse(res, crit, path_nodes, nodes, yld):
    """Critical-node table, plastic zone and path data from the final solve."""
    frd = res["frd"]
    times = sorted(frd)
    hist = []
    for t in times:
        s, e = frd[t]["S"][crit], frd[t]["E"][crit]
        F = [f for tt, f in res["dat"] if abs(tt - t) < 1e-6]
        hist.append({"time": t, "force_N": F[0] if F else np.nan, "sxx": s[0], "vm": von_mises(s), "exx": e[0], "peeq": frd[t]["PE"].get(crit, 0.0)})
    h = pd.DataFrame(hist)
    r = p.HOLE_DIAMETER / 2
    path_sorted = sorted(path_nodes, key=lambda n: nodes[n][1])
    states = {}
    for key, t in (("c1_max", 1.0), ("c1_min", 2.0), ("c3_max", 5.0), ("c3_min", 6.0)):
        tt = min(times, key=lambda x: abs(x - t))
        frd_t = frd[tt]
        peeq_path = np.array([frd_t["PE"].get(n, 0.0) for n in path_sorted])
        y = np.array([nodes[n][1] for n in path_sorted]) - r
        zone = float(y[peeq_path > PEEQ_TOL].max()) if (peeq_path > PEEQ_TOL).any() else 0.0
        nxt = float(y[y > zone + 1e-9].min()) if (y > zone + 1e-9).any() else np.nan
        states[key] = {"time": tt, "path_y_minus_r": y, "path_sxx": np.array([frd_t["S"][n][0] for n in path_sorted]),
                       "path_peeq": peeq_path, "zone_mm": zone, "next_node_mm": nxt}
    return h, states


def row_at(h, t):
    return h.iloc[(h.time - t).abs().argmin()]


def case_worker(args):
    """Displacement search for one case (and the NLGEOM check for the 750 MPa case); analysis is done in main."""
    yld, head, outn, u_el, k_el = args
    res, log = search(yld, head, outn, u_el, k_el)
    res_nl = None
    if yld == 750:
        res_nl = solve("yield750_ASSUMED_nlgeom", head, outn, 750, res["legs"][0], res["legs"][1], nlgeom=True)
    return yld, res, log, res_nl


# ---------------------------------------------------------------- limit state with explicit sa, sm
def g_general(sa, sm, uts, ratio, crit, kf=ft.KT_B):
    sd = ratio * uts / kf
    if crit == "swt":
        return sd - np.sqrt((sa + sm) * sa)
    k = sm / uts
    return (sd * (1 - k) if crit == "goodman" else sd * (1 - k**2)) - sa


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summ = json.load(open(p.RESULTS_DIR / "stress_field" / "summary.json"))
    k_el = summ["force_at_unit_disp_N"] / p.APPLIED_DISP
    u_el = summ["scale_to_1200N"] * p.APPLIED_DISP
    print(f"elastic stiffness (stress_field) {k_el:.3f} N/mm, displacement for 1200 N {u_el:.4f} mm; quarter model force scaling x4 (RF of LOAD set)")
    head, outn, nodes, hole, n_nodes = head_and_sets()
    print(f"mesh: {n_nodes} nodes, output set {len(outn)} nodes, hole size {HOLE_SIZE} mm")
    r, tol = p.HOLE_DIAMETER / 2, 1e-5
    path_nodes = [n for n in outn if abs(nodes[n][0]) < tol and abs(nodes[n][2]) < tol and nodes[n][1] >= r - tol]

    # ---- elastic reference and the three elastic-plastic cases run in parallel (search guesses from stress_field.py)
    with Pool(4) as pool:
        out = pool.map(case_worker, [(y, head, outn, u_el, k_el) for y in (None,) + YIELDS])
    res_e = out[0][1]
    frd_e = res_e["frd"]
    t1 = min(frd_e, key=lambda x: abs(x - 1.0))
    tens = [n for n in hole if abs(nodes[n][2]) < tol and n in frd_e[t1]["S"]]
    crit = max(tens, key=lambda n: frd_e[t1]["S"][n][0])
    print(f"critical node {crit} at {nodes[crit]}; elastic sxx at the first peak {frd_e[t1]['S'][crit][0]:.2f} MPa "
          f"at F = {force_at(res_e, 1.0):.2f} N (stress_field.py, scaled to 1200 N: {summ['peak_sxx_MPa']:.2f})")
    h_e, st_e = analyse(res_e, crit, path_nodes, nodes, None)
    h_e.to_csv(OUT / "crit_node_history_elastic.csv", index=False)
    cases, logs, res_nl = {}, list(out[0][2]), None
    for y, res, log, rnl in out[1:]:
        h, st = analyse(res, crit, path_nodes, nodes, y)
        h.to_csv(OUT / f"crit_node_history_{label(y)}.csv", index=False)
        cases[y] = dict(name=res["name"], wall=res["wall_s"], warns=res["warnings"], ninc=res["n_increments"], legs=res["legs"], h=h, st=st)
        logs += log
        if rnl is not None:
            res_nl = rnl
    searchlog = pd.DataFrame(logs)
    searchlog.to_csv(OUT / "displacement_search.csv", index=False); print(searchlog.round(3).to_string(index=False))

    # ---- NLGEOM check, 750 MPa: same displacements as the final linear-geometry solve
    c750 = cases[750]
    h_nl, _ = analyse(res_nl, crit, path_nodes, nodes, 750)
    h_nl.to_csv(OUT / "crit_node_history_yield750_ASSUMED_nlgeom.csv", index=False)
    rows = []
    for lab, t in (("cycle 1 max", 1.0), ("cycle 1 min", 2.0), ("cycle 3 max", 5.0), ("cycle 3 min", 6.0)):
        a, b = row_at(c750["h"], t), row_at(h_nl, t)
        rows.append({"state": lab, "F_linear_geom_N": a.force_N, "F_nlgeom_N": b.force_N,
                     "sxx_linear_geom": a.sxx, "sxx_nlgeom": b.sxx, "sxx_diff_pct": 100 * (b.sxx / a.sxx - 1),
                     "sxx_nlgeom_scaled_to_same_force": b.sxx * a.force_N / b.force_N,
                     "diff_pct_same_force": 100 * (b.sxx * a.force_N / b.force_N / a.sxx - 1),
                     "vm_linear_geom": a.vm, "vm_nlgeom": b.vm, "peeq_linear_geom": a.peeq, "peeq_nlgeom": b.peeq})
    nl = pd.DataFrame(rows); nl.to_csv(OUT / "nlgeom_check_yield750_ASSUMED.csv", index=False)
    print("\nNLGEOM check (750 MPa ASSUMED)"); print(nl.round(4).to_string(index=False))
    a5, a6, b5, b6 = row_at(c750["h"], 5.0), row_at(c750["h"], 6.0), row_at(h_nl, 5.0), row_at(h_nl, 6.0)
    am = pd.DataFrame([{"geometry": "linear", "amplitude_c3": (a5.sxx - a6.sxx) / 2, "mean_c3": (a5.sxx + a6.sxx) / 2, "F_max_N": a5.force_N},
                       {"geometry": "NLGEOM", "amplitude_c3": (b5.sxx - b6.sxx) / 2, "mean_c3": (b5.sxx + b6.sxx) / 2, "F_max_N": b5.force_N}])
    am["amplitude_vs_linear_pct"] = 100 * (am.amplitude_c3 / am.amplitude_c3.iloc[0] - 1)
    am["mean_vs_linear_pct"] = 100 * (am.mean_c3 / am.mean_c3.iloc[0] - 1)
    am.to_csv(OUT / "nlgeom_amplitude_mean_yield750_ASSUMED.csv", index=False); print(am.round(4).to_string(index=False))

    # ---- critical-node table (item 4)
    rows = []
    for y in (None,) + YIELDS:
        h = h_e if y is None else cases[y]["h"]; st = st_e if y is None else cases[y]["st"]
        r1, r2, r5, r6 = (row_at(h, t) for t in (1.0, 2.0, 5.0, 6.0))
        amp, mean = (r5.sxx - r6.sxx) / 2, (r5.sxx + r6.sxx) / 2
        rows.append({"case": label(y), "F_c1_max": r1.force_N, "F_c1_min": r2.force_N, "F_c3_max": r5.force_N, "F_c3_min": r6.force_N,
                     "sxx_c1_max": r1.sxx, "sxx_c1_min": r2.sxx, "sxx_c3_max": r5.sxx, "sxx_c3_min": r6.sxx,
                     "vm_c1_max": r1.vm, "vm_c1_min": r2.vm, "vm_c3_max": r5.vm, "vm_c3_min": r6.vm,
                     "amplitude_c3_sxx": amp, "mean_c3_sxx": mean,
                     "peeq_after_c1": r2.peeq, "peeq_after_c3": r6.peeq, "peeq_at_c1_max": r1.peeq, "peeq_at_c3_max": r5.peeq,
                     "plastic_zone_mm_c1_max": st["c1_max"]["zone_mm"], "plastic_zone_mm_c3_max": st["c3_max"]["zone_mm"],
                     "next_nonplastic_node_mm": st["c3_max"]["next_node_mm"],
                     "solver_wall_s_final_solve": (res_e["wall_s"] if y is None else cases[y]["wall"]),
                     "n_increments": (res_e["n_increments"] if y is None else cases[y]["ninc"])})
    tab = pd.DataFrame(rows)
    # elastic reference scaled linearly to exactly 1200 / 120 N for the comparison (linear model)
    e = tab[tab.case == "elastic"].iloc[0]
    tab["sxx_c3_max_vs_elastic_pct"] = 100 * (tab.sxx_c3_max / e.sxx_c3_max - 1)
    tab["amplitude_vs_elastic_pct"] = 100 * (tab.amplitude_c3_sxx / e.amplitude_c3_sxx - 1)
    tab["mean_vs_elastic_pct"] = 100 * (tab.mean_c3_sxx / e.mean_c3_sxx - 1)
    tab["local_R_c3"] = tab.sxx_c3_min / tab.sxx_c3_max
    tab["f_cycle1"] = tab.sxx_c1_max / e.sxx_c1_max
    tab["f_cycle3"] = tab.sxx_c3_max / e.sxx_c3_max
    tab.to_csv(OUT / "critical_node_table.csv", index=False)
    print("\nCritical node"); print(tab.T.to_string())
    zrows = []
    for y in YIELDS:
        for k, v in cases[y]["st"].items():
            zrows.append({"case": label(y), "state": k, "zone_mm_along_net_section_tension_face": v["zone_mm"],
                          "next_node_beyond_zone_mm": v["next_node_mm"], "max_peeq_on_path": float(v["path_peeq"].max())})
    pd.DataFrame(zrows).to_csv(OUT / "plastic_zone.csv", index=False)
    pd.DataFrame(zrows).pipe(lambda d: print("\nplastic zone\n" + d.round(4).to_string(index=False)))
    path_rows = []
    for y in (None,) + YIELDS:
        st = st_e if y is None else cases[y]["st"]
        for k in ("c1_max", "c3_max", "c3_min"):
            for yy, sx, pe in zip(st[k]["path_y_minus_r"], st[k]["path_sxx"], st[k]["path_peeq"]):
                path_rows.append({"case": label(y), "state": k, "dist_from_hole_edge_mm": yy, "sxx": sx, "peeq": pe})
    pd.DataFrame(path_rows).to_csv(OUT / "net_section_path.csv", index=False)

    # ---- fatigue consequence (item 5)
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    gp, q2, rmse = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
    rng = np.random.default_rng(1)
    xs = ft.to_physical(rng.random((400_000, len(ft.NAMES))))
    peak = ft.peak_stress(gp, xs); uts, ratio = xs[:, -2], xs[:, -1]
    sa_el, sm_el = (1 - ft.R) / 2 * peak / ft.KT_B, (1 + ft.R) / 2 * peak / ft.KT_B
    prow = []
    for crit_name in ("goodman", "gerber", "swt"):
        base = (ft.g_fun(peak, uts, ratio, crit_name) < 0).mean()
        chk = (g_general(sa_el, sm_el, uts, ratio, crit_name) < 0).mean()
        assert abs(base - chk) < 1e-12, "g_general does not reproduce g_fun"
        prow.append({"case": "elastic (baseline)", "criterion": crit_name, "variant": "-", "Pf": base, "dPf_vs_baseline": 0.0,
                     "mc_se": np.sqrt(base * (1 - base) / len(peak))})
        for y in YIELDS:
            r = tab[tab.case == label(y)].iloc[0]
            f3 = r.f_cycle3
            v1 = (ft.g_fun(peak * f3, uts, ratio, crit_name) < 0).mean()
            fa, fm = r.amplitude_c3_sxx / e.amplitude_c3_sxx, r.mean_c3_sxx / e.mean_c3_sxx
            v2 = (g_general(fa * sa_el, fm * sm_el, uts, ratio, crit_name) < 0).mean()
            for var, v, note in (("V_f: peak x f(cycle 3)", v1, f"f={f3:.4f}"), ("V_loc: amplitude, mean from cycle 3", v2, f"amp x{fa:.4f}, mean x{fm:.4f}")):
                prow.append({"case": label(y), "criterion": crit_name, "variant": var, "Pf": v, "dPf_vs_baseline": v - base,
                             "mc_se": np.sqrt(v * (1 - v) / len(peak)), "factors": note})
    pf = pd.DataFrame(prow); pf.to_csv(OUT / "pf_exploratory_ASSUMED.csv", index=False)
    print("\nPf exploratory (ASSUMED)"); print(pf.round(4).to_string(index=False))

    # ---------------- figures ----------------
    cols = {None: "k", 650: "#4c78a8", 750: "#f58518", 850: "#54a24b"}
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for y in (None,) + YIELDS:
        st = st_e if y is None else cases[y]["st"]
        d = st["c3_max"]
        ax.plot(d["path_y_minus_r"], d["path_sxx"], "-o", ms=3, color=cols[y], label="elastic" if y is None else f"yield {y} MPa (ASSUMED)")
        if y is not None and d["zone_mm"] > 0:
            ax.axvline(d["zone_mm"], color=cols[y], ls=":", lw=1)
    ax.set_xlim(0, 3); ax.set_xlabel("distance from the hole edge along the net section, tension face (mm)")
    ax.set_ylabel("sigma_xx at maximum load, cycle 3 (MPa)")
    ax.set_title("dotted lines: outermost node with PEEQ > 1e-5", fontsize=8); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig_sxx_net_section_ASSUMED_yields.png", dpi=600); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    for a, y in zip(axs, YIELDS):
        a.plot(h_e.exx * 1e3, h_e.sxx, color="0.5", ls="--", lw=1, label="elastic")
        h = cases[y]["h"]
        a.plot(h.exx * 1e3, h.sxx, "-o", ms=2, color=cols[y], label=f"yield {y} MPa (ASSUMED)")
        a.set_xlabel("total strain exx at the critical node (1e-3)"); a.legend(fontsize=8)
    axs[0].set_ylabel("sigma_xx (MPa)")
    fig.tight_layout(); fig.savefig(OUT / "fig_stress_strain_loops_ASSUMED_yields.png", dpi=600); plt.close(fig)

    json.dump({"critical_node": crit, "critical_node_xyz": nodes[crit], "n_nodes": n_nodes, "yields_ASSUMED": list(YIELDS),
               "peeq_tol_ASSUMED": PEEQ_TOL, "uts_end_ASSUMED": UTS_END, "eps_end_ASSUMED": EPS_END, "ninc_per_leg": NINC,
               "warnings": {label(y): cases[y]["warns"] for y in YIELDS} | {"elastic": res_e["warnings"], "nlgeom": res_nl["warnings"]},
               "increments": {label(y): cases[y]["ninc"] for y in YIELDS} | {"elastic": res_e["n_increments"], "nlgeom": res_nl["n_increments"]},
               "wall_s": {label(y): cases[y]["wall"] for y in YIELDS} | {"elastic": res_e["wall_s"], "nlgeom": res_nl["wall_s"]}},
              open(OUT / "meta.json", "w"), indent=2)
    print("\nwarnings / wall times:"); print(json.load(open(OUT / "meta.json"))["warnings"], json.load(open(OUT / "meta.json"))["wall_s"])


if __name__ == "__main__":
    main()
