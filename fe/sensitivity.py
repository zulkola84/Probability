"""Surrogate-based Monte Carlo and Sobol sensitivity of peak von Mises stress and stiffness.

Reads results/fe_samples.csv, fits a Gaussian-process surrogate per output (checked by
5-fold cross-validation), propagates the input distributions in uncertainty.py with
Monte Carlo, and computes Sobol indices with SALib. Peak stress scales linearly with the
applied force, so force enters analytically: stress = surrogate(x) * force / 1200.

Writes results/sensitivity_summary.csv, results/sobol_indices.csv and results/sensitivity.png.
Usage: python sensitivity.py [N (power of 2, default 8192)]
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SALib.analyze import sobol
from SALib.sample import sobol as sobol_sample
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.model_selection import cross_val_predict

import params as p
import uncertainty as unc

OUTPUTS = {"peak_vm": "Peak von Mises at hole edge (MPa)", "stiffness": "Stiffness (N/mm)"}
# surrogate inputs are scaled by the nominal distribution mean and sd
MU = np.array([unc.FE_INPUTS[k].mean() for k in unc.FE_NAMES])
SD = np.array([unc.FE_INPUTS[k].std() for k in unc.FE_NAMES])


def fit_surrogate(x, y):
    kernel = ConstantKernel(1.0) * RBF(np.ones(x.shape[1]), (1e-2, 1e3)) + WhiteKernel(1e-4, (1e-8, 1e-1))
    gp = GaussianProcessRegressor(kernel, normalize_y=True, n_restarts_optimizer=3, random_state=0)
    xs = (x - MU) / SD
    pred = cross_val_predict(gp, xs, y, cv=5)
    q2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
    rmse = np.sqrt(np.mean((y - pred) ** 2))
    return gp.fit(xs, y), q2, rmse


def model(gps, u):
    """Evaluate both outputs for unit-cube samples u (n, 7): 6 FE inputs + force."""
    x = unc.to_physical(u[:, :-1])
    force = unc.FORCE.ppf(np.clip(u[:, -1], 1e-6, 1 - 1e-6))
    xs = (x - MU) / SD
    return {
        "peak_vm": gps["peak_vm"].predict(xs) * force / p.RUNOUT_FORCE,
        "stiffness": gps["stiffness"].predict(xs),
    }


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8192
    df = pd.read_csv(p.RESULTS_DIR / "fe_samples.csv")
    x = df[unc.FE_NAMES].to_numpy()

    gps, rows = {}, []
    for out in OUTPUTS:
        gps[out], q2, rmse = fit_surrogate(x, df[out].to_numpy())
        print(f"{out}: surrogate Q2={q2:.4f}  CV RMSE={rmse:.3g}")
        rows.append({"output": out, "q2": q2, "cv_rmse": rmse})

    problem = {"num_vars": len(unc.ALL_NAMES), "names": unc.ALL_NAMES,
               "bounds": [[0, 1]] * len(unc.ALL_NAMES)}
    u = sobol_sample.sample(problem, n, calc_second_order=False, seed=1)
    y = model(gps, u)

    indices = []
    summary = []
    for out in OUTPUTS:
        s = sobol.analyze(problem, y[out], calc_second_order=False, seed=1)
        for i, name in enumerate(unc.ALL_NAMES):
            indices.append({"output": out, "input": name, "S1": s["S1"][i], "S1_conf": s["S1_conf"][i],
                            "ST": s["ST"][i], "ST_conf": s["ST_conf"][i]})
        v = y[out]
        summary.append({"output": out, "mean": v.mean(), "sd": v.std(), "p05": np.percentile(v, 5),
                        "p50": np.percentile(v, 50), "p95": np.percentile(v, 95)})
    idx = pd.DataFrame(indices)
    summ = pd.DataFrame(summary).merge(pd.DataFrame(rows), on="output")

    # reliability-style numbers for the stress output
    v = y["peak_vm"]
    summ.loc[summ.output == "peak_vm", "P(peak > paper FE 857.8)"] = np.mean(v > p.PAPER_PEAK_FE)
    summ.loc[summ.output == "peak_vm", "P(peak > UTS 900)"] = np.mean(v > p.UTS)

    idx.to_csv(p.RESULTS_DIR / "sobol_indices.csv", index=False)
    summ.to_csv(p.RESULTS_DIR / "sensitivity_summary.csv", index=False)
    print(summ.round(4).to_string(index=False))
    print(idx.round(3).to_string(index=False))

    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
    for j, out in enumerate(OUTPUTS):
        a = ax[0, j]
        a.hist(y[out], bins=60, color="#4c78a8")
        a.set_xlabel(OUTPUTS[out]); a.set_ylabel("count")
        if out == "peak_vm":
            a.axvline(p.PAPER_PEAK_FE, color="k", ls="--", label="paper FE 857.8")
            a.axvline(p.UTS, color="r", ls=":", label="UTS 900")
            a.legend()
        b = ax[1, j]
        d = idx[idx.output == out]
        pos = np.arange(len(d))
        b.bar(pos - 0.2, d.S1, 0.4, yerr=d.S1_conf, label="S1", color="#4c78a8")
        b.bar(pos + 0.2, d.ST, 0.4, yerr=d.ST_conf, label="ST", color="#f58518")
        b.set_xticks(pos); b.set_xticklabels(d.input, rotation=30)
        b.set_ylabel("Sobol index"); b.legend()
    fig.tight_layout()
    fig.savefig(p.RESULTS_DIR / "sensitivity.png", dpi=150)


if __name__ == "__main__":
    main()
