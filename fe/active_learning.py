"""Active learning of the peak-stress surrogate near the fatigue failure boundary g = 0.

g depends on the FE model only through the peak von Mises stress (see fatigue.py), so the FE
training points live in the five FE inputs. Each iteration, for the Goodman and the Gerber criterion:
    1. draw a Monte Carlo pool of all random inputs and evaluate g with the current surrogate;
    2. rank pool points by the AK-MCS learning function U = |g| / sigma_g (small U = the surrogate
       cannot tell which side of g = 0 the point is on), sigma_g from the surrogate's predictive sd;
    3. take a space-filling batch of the lowest-U points, run the FE model there, and compare with the
       surrogate's prediction BEFORE retraining;
    4. retrain on the enlarged set and recompute the failure probability.
Stops when min U >= U_STOP for both criteria, or after MAX_ITER iterations (Stage A brief). Caveat: a
continuous pool always holds a point near g = 0, so min U tends to stay small; the mean misclassification
probability mean(Phi(-U)) and the P_f change between iterations are logged as diagnostics, not stop rules.

Starts from results/fe_samples_initial.csv (the original 256-run DOE) and writes ONLY to results/stageA/;
no existing result file is modified.
Usage: python active_learning.py [batch_per_criterion (default 8)] [workers (default 8)]
Writes results/stageA/{al_log.csv, al_new_points.csv, fe_samples_al.csv, pf_vs_fe_runs.png}
"""
import sys
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import fatigue as ft
import params as p
import run_samples
import sensitivity as sens
import uncertainty as unc

N_POOL = 150_000          # candidate pool per criterion
N_PF = 200_000            # Monte Carlo size for the failure probability inside the loop
MAX_ITER = 10
U_STOP = 2.0
OUT_DIR = p.RESULTS_DIR / "stageA"
MIN_DIST = 0.35           # min distance between chosen points, in units of input sd
N_FE = len(unc.FE_NAMES)


def predict(gp, x_fe, return_std=False):
    xs = (x_fe - sens.MU) / sens.SD
    mean, sd = [], []
    for i in range(0, len(xs), 50_000):
        m, s = gp.predict(xs[i:i + 50_000], return_std=True)
        mean.append(m); sd.append(s)
    return (np.concatenate(mean), np.concatenate(sd)) if return_std else np.concatenate(mean)


def learning_function(gp, x, crit):
    """U = |g| / sigma_g for physical samples x (n, 8)."""
    pk, sd = predict(gp, x[:, :N_FE], return_std=True)
    scale = x[:, N_FE] / p.RUNOUT_FORCE
    g = ft.g_fun(pk * scale, x[:, -2], x[:, -1], crit)
    slope = (ft.g_fun(pk * scale + 1.0, x[:, -2], x[:, -1], crit) - g)     # dg per MPa of peak stress
    return g, np.abs(g) / np.maximum(np.abs(slope) * sd * scale, 1e-12)


def pick(x_fe, u, n, existing):
    """Greedy space-filling choice of n low-U points, away from each other and from existing samples."""
    chosen, keep = [], list(existing)
    for i in np.argsort(u):
        z = (x_fe[i] - sens.MU) / sens.SD
        if all(np.linalg.norm(z - k) > MIN_DIST for k in keep):
            chosen.append(i); keep.append(z)
        if len(chosen) == n:
            break
    return chosen


def failure_probability(gp, crit):
    pf, _ = ft.mc(gp, N_PF, crit, np.random.default_rng(7))
    return pf


def main():
    batch = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    OUT_DIR.mkdir(exist_ok=True)
    df = pd.read_csv(p.RESULTS_DIR / "fe_samples_initial.csv")
    n0 = len(df)
    pts = []
    rng = np.random.default_rng(11)
    log, history = [], []

    for it in range(MAX_ITER + 1):
        gp, q2, rmse = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
        pf = {c: failure_probability(gp, c) for c in ("goodman", "gerber")}
        pool = ft.to_physical(rng.random((N_POOL, len(ft.NAMES))))
        gu = {c: learning_function(gp, pool, c) for c in pf}
        min_u = {c: float(gu[c][1].min()) for c in pf}
        misclass = {c: float(np.mean(stats.norm.cdf(-gu[c][1]))) for c in pf}
        row = {"iter": it, "n_train": len(df), "q2": q2, "cv_rmse": rmse,
               "pf_goodman": pf["goodman"], "pf_gerber": pf["gerber"],
               "min_u_goodman": min_u["goodman"], "min_u_gerber": min_u["gerber"],
               "misclass_goodman": misclass["goodman"], "misclass_gerber": misclass["gerber"]}
        history.append(pf)
        print(f"iter {it}: n={len(df)} Q2={q2:.4f} Pf goodman={pf['goodman']:.4f} gerber={pf['gerber']:.4f} "
              f"misclass {misclass['goodman']:.4f}/{misclass['gerber']:.4f}")
        if it == MAX_ITER or all(min_u[c] >= U_STOP for c in pf):
            log.append(row)
            break

        existing = [(r - sens.MU) / sens.SD for r in df[unc.FE_NAMES].to_numpy()]
        new_idx = []
        for c in pf:
            idx = pick(pool[:, :N_FE], gu[c][1], batch, existing + [(pool[i, :N_FE] - sens.MU) / sens.SD for i in new_idx])
            new_idx += idx
        x_new = pool[new_idx, :N_FE]
        with Pool(workers) as pl:
            out = pl.map(run_samples.run_one, x_new, chunksize=1)
        fe_peak = np.array([o[0] for o in out]); fe_k = np.array([o[1] for o in out])
        err = predict(gp, x_new) - fe_peak
        fe_pk = pd.DataFrame({"iter": it, "criterion": np.repeat(list(pf), batch)[:len(x_new)],
                              "fe_peak_MPa": fe_peak, "surrogate_MPa": fe_peak + err})
        fe_pk["err_MPa"] = err; fe_pk["err_pct_of_peak"] = 100 * err / fe_peak
        pts.append(fe_pk)
        row.update({"new_points": len(x_new), "surrogate_bias_MPa": err.mean(),
                    "surrogate_rmse_MPa": np.sqrt(np.mean(err**2)), "surrogate_max_abs_MPa": np.abs(err).max(),
                    "surrogate_max_abs_pct_of_peak": np.max(np.abs(err) / fe_peak) * 100})
        print(f"   FE vs surrogate at {len(x_new)} boundary points: bias {err.mean():+.2f}, "
              f"rmse {row['surrogate_rmse_MPa']:.2f}, max |err| {np.abs(err).max():.2f} MPa")
        log.append(row)
        add = pd.DataFrame(x_new, columns=unc.FE_NAMES)
        add["peak_vm"], add["stiffness"] = fe_peak, fe_k
        df = pd.concat([df, add], ignore_index=True)

    df.to_csv(OUT_DIR / "fe_samples_al.csv", index=False)
    lg = pd.DataFrame(log)
    lg.to_csv(OUT_DIR / "al_log.csv", index=False)
    if pts:
        pd.concat(pts, ignore_index=True).to_csv(OUT_DIR / "al_new_points.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(lg.n_train, lg.pf_goodman, "-o", label="Goodman")
    ax.plot(lg.n_train, lg.pf_gerber, "-s", label="Gerber")
    ax.set_xlabel("number of FE runs"); ax.set_ylabel("P_f at 1200 N"); ax.legend()
    fig.tight_layout(); fig.savefig(OUT_DIR / "pf_vs_fe_runs.png", dpi=150)
    print(f"{len(df)} samples (initial {n0}) written to {OUT_DIR}")


if __name__ == "__main__":
    main()
