"""Stage C driver: F(delta) curve, stress at 1200 N, mesh convergence and sensitivity runs for stage_c.py.
Everything is written to results/stageC/. Usage: python stage_c_run.py
Sensitivity cases vary ASSUMED inputs (contact stiffness, orientation, R_out) and the reading of I (net vs gross).
"""
import tempfile
from multiprocessing import Pool

import numpy as np
import pandas as pd

import params as p
import stage_c as c

OUT = p.RESULTS_DIR / "stageC"
BASE = dict(r_out=25.0, convex_up=False, K=1.0e6, hole=0.55, net=True)


def task(t):
    cfg = {**BASE, **t["cfg"]}
    c.CONTACT_K = cfg["K"]
    geo = c.geometry(cfg["r_out"], cfg["convex_up"], cfg["net"])
    m = c.build_mesh(geo, cfg["hole"])
    with tempfile.TemporaryDirectory() as tmp:
        r = c.run(m, t["disp"], tmp)
    r.pop("stdout_tail")
    return {"case": t["case"], **cfg, "alpha_deg": np.degrees(geo["alpha"]), "bbox_width": 2 * geo["ye"],
            "nodes": m["n_nodes"], **r}


def main():
    OUT.mkdir(exist_ok=True)
    tasks = [{"case": "base", "cfg": {}, "disp": d} for d in (1.0, 2.0, 3.0)]
    tasks += [{"case": "conv", "cfg": {"hole": h}, "disp": 2.0} for h in (1.2, 0.8, 0.4)]
    tasks += [{"case": "sens_K1e5", "cfg": {"K": 1e5}, "disp": 2.0}, {"case": "sens_K1e7", "cfg": {"K": 1e7}, "disp": 2.0},
              {"case": "sens_convex_up", "cfg": {"convex_up": True}, "disp": 2.0},
              {"case": "sens_R15", "cfg": {"r_out": 15.0}, "disp": 2.0},
              {"case": "sens_R60", "cfg": {"r_out": 60.0}, "disp": 2.0},
              {"case": "sens_gross_I", "cfg": {"net": False}, "disp": 2.0}]
    with Pool(4) as pool:
        rows = pool.map(task, tasks, chunksize=1)
    base = [r for r in rows if r["case"] == "base"]
    d, F = np.array([r["disp"] for r in base]), np.array([r["force"] for r in base])
    k, b = np.polyfit(d, F, 1)
    d_star = (p.RUNOUT_FORCE - b) / k
    rows.append(task({"case": "base_at_1200N", "cfg": {}, "disp": float(d_star)}))
    df = pd.DataFrame(rows)
    df["peak_vm_scaled_to_1200N"] = df.peak_vm * p.RUNOUT_FORCE / df.force
    df["K_secant"] = df.force / df.disp
    df.to_csv(OUT / "stageC_runs.csv", index=False)
    print(f"linear fit 1-3 mm: K = {k:.1f} N/mm, intercept {b:.1f} N, delta at 1200 N = {d_star:.3f} mm")
    print(df.drop(columns=["net", "convex_up"]).round(3).to_string())


if __name__ == "__main__":
    main()
