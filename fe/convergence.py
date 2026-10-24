"""Mesh convergence study: refine the hole-edge element size and track K and peak stress.

Writes results/convergence.csv. Usage: python convergence.py
"""
import csv

import make_mesh
import params as p
import run_ccx

HOLE_SIZES = [1.2, 0.8, 0.55, 0.4, 0.3]
GLOBAL_SIZE = 3.0


def main():
    p.RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for hs in HOLE_SIZES:
        name = f"conv_h{hs}"
        _, n_nodes = make_mesh.build(hs, GLOBAL_SIZE, name)
        r = run_ccx.run(name)
        rows.append({"hole_size_mm": hs, "nodes": n_nodes, **r})
        print(f"h={hs:<5} nodes={n_nodes:<7} K={r['stiffness']:.1f} N/mm "
              f"peak={r['peak']:.1f} vM={r['peak_vm']:.1f} MPa Kt_net={r['kt']:.3f}")

    for prev, cur in zip(rows, rows[1:]):
        cur["peak_change_pct"] = 100 * (cur["peak"] - prev["peak"]) / prev["peak"]
    with open(p.RESULTS_DIR / "convergence.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["hole_size_mm", "nodes", "stiffness",
                                          "peak", "peak_vm", "nominal", "kt", "kt_gross", "nominal_I", "peak_change_pct"])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
