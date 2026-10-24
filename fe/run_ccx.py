"""Run CalculiX on a deck from make_mesh.py and extract stiffness and hole-edge peak stress.

The model is linear, so results are scaled from the applied displacement to the
runout force. Stresses are CalculiX nodal values (extrapolated and averaged).

Usage: python run_ccx.py <name>
"""
import re
import subprocess
import sys

import numpy as np

import params as p


def read_frd_stress(path):
    """Return {node: [sxx, syy, szz, sxy, syz, szx]} from the last STRESS block."""
    stress, active = {}, False
    for line in open(path):
        if line.startswith(" -4") and "STRESS" in line:
            stress, active = {}, True
        elif active and line.startswith(" -3"):
            active = False
        elif active and line.startswith(" -1"):
            vals = [float(v) for v in re.findall(r"-?\d\.\d{5}E[+-]\d{2}", line[13:])]
            stress[int(line[3:13])] = vals
    return stress


def read_hole_nodes(inp):
    nodes, on = [], False
    for line in open(inp):
        if line.startswith("*"):
            on = line.strip().upper() == "*NSET, NSET=HOLE"
        elif on:
            nodes += [int(v) for v in line.split(",") if v.strip()]
    return nodes


def read_total_rf(dat):
    """Sum of RF_z over the LOAD node set (quarter model)."""
    for line in open(dat):
        m = re.match(r"\s*([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s*$", line)
        if m:
            rf = float(m.group(3))
    return rf


def max_principal(s):
    sxx, syy, szz, sxy, syz, szx = s
    m = np.array([[sxx, sxy, szx], [sxy, syy, syz], [szx, syz, szz]])
    return np.linalg.eigvalsh(m)[-1]


def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = s
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
                   + 3 * (sxy**2 + syz**2 + szx**2))


def run(name, out_dir=None):
    out_dir = out_dir or p.FE_DIR
    subprocess.run(["ccx", name], cwd=out_dir, check=True, capture_output=True)

    # quarter model holds 1/4 of the total force over both loading rollers
    force = 4 * abs(read_total_rf(out_dir / f"{name}.dat"))
    stiffness = force / p.APPLIED_DISP

    stress = read_frd_stress(out_dir / f"{name}.frd")
    hole = read_hole_nodes(out_dir / f"{name}.inp")
    peak = max(max_principal(stress[n]) for n in hole if n in stress)
    peak *= p.RUNOUT_FORCE / force
    peak_vm = max(von_mises(stress[n]) for n in hole if n in stress)
    peak_vm *= p.RUNOUT_FORCE / force

    # nominal bending stress on the net section at the hole
    moment = p.RUNOUT_FORCE / 2 * p.LEVER_ARM
    nominal = 6 * moment / ((p.WIDTH - p.HOLE_DIAMETER) * p.THICKNESS**2)
    gross = 6 * moment / (p.WIDTH * p.THICKNESS**2)
    return {"stiffness": stiffness, "peak": peak, "peak_vm": peak_vm, "nominal": nominal, "kt": peak / nominal,
            "kt_gross": peak / gross,
            "nominal_I": moment * p.THICKNESS / 2 / p.I_BENDING}


if __name__ == "__main__":
    print(run(sys.argv[1]))
