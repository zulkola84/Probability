"""Run the FE model over a Sobol-sequence design of the uncertain inputs.

Writes results/fe_samples.csv (inputs + peak von Mises at the runout force + stiffness).
Usage: python run_samples.py [n_samples (power of 2)] [workers]
"""
import csv
import os
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.stats import qmc

import make_mesh
import params as p
import run_ccx
import uncertainty as unc

HOLE_SIZE, GLOBAL_SIZE = 0.55, 3.0


def run_one(x):
    """Run one FE sample; x is a physical input vector in FE_NAMES order."""
    v = dict(zip(unc.FE_NAMES, x))
    p.THICKNESS, p.HOLE_DIAMETER = v["thickness"], v["hole_d"]
    p.E, p.NU, p.LEVER_ARM = v["E"], v["nu"], v["lever_arm"]
    p.SUPPORT_SPAN = p.INNER_SPAN + 2 * p.LEVER_ARM
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        make_mesh.build(HOLE_SIZE, GLOBAL_SIZE, "s", out_dir=tmp)
        r = run_ccx.run("s", out_dir=tmp)
    return float(r["peak_vm"]), float(r["stiffness"])


def main():
    os.environ["OMP_NUM_THREADS"] = "1"
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 128
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    u = qmc.Sobol(d=len(unc.FE_NAMES), scramble=True, seed=1).random(n)
    x = unc.to_physical(u)
    with Pool(workers) as pool:
        out = pool.map(run_one, x, chunksize=2)
    p.RESULTS_DIR.mkdir(exist_ok=True)
    with open(p.RESULTS_DIR / "fe_samples.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(unc.FE_NAMES + ["peak_vm", "stiffness"])
        for xi, (pk, k) in zip(x, out):
            w.writerow(list(np.round(xi, 6)) + [round(pk, 3), round(k, 3)])
    print(f"{n} samples written")


if __name__ == "__main__":
    main()
