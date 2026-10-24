"""Stage B: independent validation on the implanted two-hole plate (Bologna et al. 2023).

Geometry unchanged from the calibrated intact-plate model (width 16.365 mm NOT retuned, thickness 5.12,
hole 7.8, same rollers and 1200 N); only a second hole is added. Hole spacing is not in the paper:
three centre-to-centre spacings are run, all ASSUMED. Paper implanted peak: 811.7 MPa FE, 768.6 MPa
analytical (ASSUMED to correspond to the same 1200 N load; the paper excerpt in inputs/ does not say).
Gate: FE peak (max von Mises on the hole surface) within 10 % of 811.7 MPa.
Each case is run at two hole-edge mesh sizes as a mesh check. Writes results/stageB/stageB_results.csv.
"""
import tempfile
from pathlib import Path

import pandas as pd

import make_mesh
import params as p
import run_ccx

PAPER_FE, PAPER_AN = 811.7, 768.6
SPACINGS = [None, 20.0, 30.0, 40.0]       # None = intact single-hole control; others ASSUMED
HOLE_SIZES = [0.55, 0.35]                 # 0.55 is the mesh used in the DOE


def main():
    rows = []
    for sp in SPACINGS:
        for hs in HOLE_SIZES:
            with tempfile.TemporaryDirectory() as tmp:
                tmp = Path(tmp)
                _, nn = make_mesh.build(hs, 3.0, "b", out_dir=tmp, hole_spacing=sp)
                r = run_ccx.run("b", out_dir=tmp)
            rows.append({"spacing_mm": sp if sp else 0, "assumed": sp is not None, "hole_size": hs, "nodes": nn,
                         "peak_vm": r["peak_vm"], "peak_principal": r["peak"], "stiffness": r["stiffness"],
                         "diff_vs_811.7_pct": (r["peak_vm"] / PAPER_FE - 1) * 100,
                         "diff_vs_768.6_pct": (r["peak_vm"] / PAPER_AN - 1) * 100,
                         "gate_pass": abs(r["peak_vm"] / PAPER_FE - 1) < 0.10 if sp else None})
            print(rows[-1])
    pd.DataFrame(rows).to_csv(p.RESULTS_DIR / "stageB" / "stageB_results.csv", index=False)


if __name__ == "__main__":
    main()
