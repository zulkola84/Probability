"""Phase 2, step 7: plate width as epistemic uncertainty.

Earlier width-sampling data (found in results/old_with_width/): 256 FE runs, width uniform 16.365 +- 1.0 mm (ASSUMED range at the
time), the five other FE inputs from the Table 2 distributions (same as uncertainty.py), scrambled Sobol (seed 1), hole-edge
size 0.55 mm, global 3.0 mm, peak von Mises at 1200 N. `verify_old` re-runs a few of them with the CURRENT code to confirm the
data are reusable; if they are, no new design is run. Then:
  - 6-input GP (thickness, hole_d, E, nu, lever_arm, width) on 416 fixed-width runs (results/stageA, width column = 16.365)
    plus the 256 width-varying runs; 5-fold CV on the width-varying runs only (the held-out folds contain only
    width-varying points; training always contains the 416 fixed-width runs).
  - Cases W0 (fixed 16.365), W1 (+-0.5 mm uniform, ASSUMED epistemic range), W2 (+-1.0 mm uniform, ASSUMED).
Usage: python width_cases.py --verify | python width_cases.py
All outputs go to OUT (outside results/).
"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SALib.analyze import sobol
from SALib.sample import sobol as sobol_sample
from scipy import stats
from scipy.optimize import brentq
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.model_selection import KFold

import fatigue as ft
import make_mesh
import params as p
import run_ccx
import sensitivity as sens
import uncertainty as unc

OUT = Path(__file__).resolve().parents[1] / "results" / "phase2" / "step7_width"
OLD = p.RESULTS_DIR / "old_with_width" / "fe_samples.csv"
W0 = 16.365
HOLE_SIZE, GLOBAL_SIZE = 0.55, 3.0
FE5 = unc.FE_NAMES                       # thickness, hole_d, E, nu, lever_arm
N_MC, N_SOBOL = 400_000, 2**13
CRITS = ("goodman", "gerber", "swt")
CASES = {"W0": 0.0, "W1": 0.5, "W2": 1.0}     # half-widths (mm); W1 and W2 ASSUMED epistemic ranges
WIDTHS_TABLE = (15.4, 16.0, 16.365, 16.8, 17.4)
Q2_MIN, RMSE_MAX_FRAC = 0.95, 0.01      # ASSUMED acceptance thresholds for "poor fit"


# ---------------------------------------------------------------- FE
def run_fe(args):
    """One FE run (current code path of run_samples.run_one, plus the width); deck kept under OUT/fe_decks/<tag>."""
    tag, w, x = args
    v = dict(zip(FE5, x))
    p.WIDTH, p.THICKNESS, p.HOLE_DIAMETER = w, v["thickness"], v["hole_d"]
    p.E, p.NU, p.LEVER_ARM = v["E"], v["nu"], v["lever_arm"]
    p.SUPPORT_SPAN = p.INNER_SPAN + 2 * p.LEVER_ARM
    d = OUT / "fe_decks" / tag
    d.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    make_mesh.build(HOLE_SIZE, GLOBAL_SIZE, "s", out_dir=d)
    r = run_ccx.run("s", out_dir=d)
    return float(r["peak_vm"]), float(r["stiffness"]), time.perf_counter() - t0


def verify_old():
    old = pd.read_csv(OLD)
    idx = [int(old.width.idxmin()), int(old.width.idxmax()), 0, 1, 50, 100, 200, 255]
    jobs = [(f"verify_{i}", float(old.width[i]), [float(old.loc[i, k]) for k in FE5]) for i in idx]
    t0 = time.perf_counter()
    with Pool(4) as pool:
        res = pool.map(run_fe, jobs)
    rows = [{"old_row": i, "width": old.width[i], "peak_old": old.peak_vm[i], "peak_rerun_current_code": r[0],
             "diff_MPa": r[0] - old.peak_vm[i], "stiffness_old": old.stiffness[i], "stiffness_rerun": r[1], "solver_s": r[2]}
            for i, r in zip(idx, res)]
    v = pd.DataFrame(rows); v.to_csv(OUT / "verify_old_runs.csv", index=False)
    print(v.round(4).to_string(index=False)); print(f"max |diff| {v.diff_MPa.abs().max():.4f} MPa; wall {time.perf_counter() - t0:.1f} s for {len(jobs)} runs")
    return v


# ---------------------------------------------------------------- GP with width
def scale6(x5, w):
    return np.column_stack([(x5 - sens.MU) / sens.SD, (np.asarray(w) - W0) / 0.5])


def fit_gp6(x5, w, y):
    k = ConstantKernel(1.0) * RBF(np.ones(6), (1e-2, 1e3)) + WhiteKernel(1e-4, (1e-8, 1e-1))
    return GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=3, random_state=0).fit(scale6(x5, w), y)


def predict6(gp, x5, w):
    s = scale6(x5, w)
    return np.concatenate([gp.predict(s[i:i + 100_000]) for i in range(0, len(s), 100_000)])


# ---------------------------------------------------------------- reliability
def build_x(u, half):
    """u: (n, 9) unit-cube samples -> (x8 as in fatigue.py, width). Width is the 9th column (fixed if half = 0)."""
    x = ft.to_physical(u[:, :8])
    w = W0 + half * (2 * u[:, 8] - 1)
    return x, w


def peak_at(gp, x, w, force=None):
    f = x[:, 5] if force is None else force
    return predict6(gp, x[:, :5], w) * f / p.RUNOUT_FORCE


def form9(gp, crit, half, iters=100):
    def g_u(z):
        x, w = build_x(stats.norm.cdf(z)[None, :], half)
        return float(ft.g_fun(peak_at(gp, x, w), x[:, -2], x[:, -1], crit)[0])
    z, done = np.zeros(9), False
    for _ in range(iters):
        g0 = g_u(z)
        grad = np.array([(g_u(z + 1e-4 * e) - g0) / 1e-4 for e in np.eye(9)])
        if grad @ grad < 1e-18:
            break
        z_new = (grad @ z - g0) / (grad @ grad) * grad
        done = np.linalg.norm(z_new - z) < 1e-6
        z = z_new
        if done:
            break
    beta = np.linalg.norm(z)
    if g_u(np.zeros(9)) < 0:
        return stats.norm.cdf(beta), -beta, done
    return stats.norm.cdf(-beta), beta, done


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if "--verify" in sys.argv:
        verify_old(); return
    old = pd.read_csv(OLD)
    ver = pd.read_csv(OUT / "verify_old_runs.csv")
    assert (ver.diff_MPa.abs() / ver.peak_old).max() < 0.01, "old width runs do not reproduce (1 % of peak, ASSUMED) with the current code; run new samples instead"
    fixed = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    print(f"found earlier width run: {OLD} ({len(old)} runs, width {old.width.min():.3f}-{old.width.max():.3f}); "
          f"verification: {len(ver)} reruns, max |diff| {ver.diff_MPa.abs().max():.3f} MPa ({100 * (ver.diff_MPa.abs() / ver.peak_old).max():.2f} % of peak), "
          f"RMS diff {np.sqrt((ver.diff_MPa ** 2).mean()):.3f} MPa, max stiffness diff {(ver.stiffness_rerun - ver.stiffness_old).abs().max():.4f} N/mm; fixed-width runs {len(fixed)}")
    xf, yf = fixed[FE5].to_numpy(), fixed["peak_vm"].to_numpy()
    xo, wo, yo = old[FE5].to_numpy(), old.width.to_numpy(), old.peak_vm.to_numpy()

    # ---------- 2. CV on width-varying points
    pred = np.empty(len(old)); fold = np.empty(len(old), int)
    t0 = time.perf_counter()
    for k, (tr, te) in enumerate(KFold(5, shuffle=True, random_state=0).split(xo)):
        gp = fit_gp6(np.vstack([xf, xo[tr]]), np.r_[np.full(len(xf), W0), wo[tr]], np.r_[yf, yo[tr]])
        pred[te] = predict6(gp, xo[te], wo[te]); fold[te] = k
        print(f"  fold {k}: RMSE {np.sqrt(np.mean((pred[te] - yo[te]) ** 2)):.3f} MPa ({len(te)} held-out width-varying points)")
    e = pred - yo
    cv = {"n_holdout_width_varying": len(old), "RMSE_MPa": float(np.sqrt(np.mean(e**2))), "bias_MPa": float(e.mean()),
          "Q2": float(1 - np.sum(e**2) / np.sum((yo - yo.mean()) ** 2)), "max_abs_MPa": float(np.abs(e).max()),
          "max_abs_pct_of_peak": float(100 * np.max(np.abs(e) / yo)), "mean_peak_MPa": float(yo.mean()),
          "cv_wall_s": time.perf_counter() - t0}
    edge = np.abs(wo - W0) > 0.75
    cv["RMSE_MPa_edge_|w-16.365|>0.75"] = float(np.sqrt(np.mean(e[edge] ** 2)))
    cv["RMSE_MPa_inner"] = float(np.sqrt(np.mean(e[~edge] ** 2)))
    pd.DataFrame({"old_row": np.arange(len(old)), "fold": fold, "width": wo, "peak_fe": yo, "peak_gp_cv": pred}).to_csv(OUT / "cv_predictions.csv", index=False)
    pd.Series(cv).to_csv(OUT / "cv_scores.csv", header=["value"]); print(pd.Series(cv).to_string())
    poor = not (cv["Q2"] >= Q2_MIN and cv["RMSE_MPa"] <= RMSE_MAX_FRAC * cv["mean_peak_MPa"])
    cv["poor_fit_flag"] = poor
    if poor:
        print(f"STOP: fit poor (Q2 {cv['Q2']:.4f}, RMSE {cv['RMSE_MPa']:.2f} MPa)"); return

    t0 = time.perf_counter()
    gp = fit_gp6(np.vstack([xf, xo]), np.r_[np.full(len(xf), W0), wo], np.r_[yf, yo])
    print("final GP:", gp.kernel_, f"({time.perf_counter() - t0:.1f} s)")
    ls = gp.kernel_.k1.k2.length_scale
    pd.DataFrame({"input": FE5 + ["width"], "length_scale_scaled": ls}).to_csv(OUT / "gp6_length_scales.csv", index=False)

    # ---------- 5. nominal peak vs width (GP) and FE check
    nom5 = np.array([[d.mean() for d in unc.FE_INPUTS.values()]])
    jobs = [(f"nominal_w{w:.3f}", w, list(nom5[0])) for w in WIDTHS_TABLE]
    with Pool(4) as pool:
        fe_nom = pool.map(run_fe, jobs)
    rows = []
    for w, r in zip(WIDTHS_TABLE, fe_nom):
        g = float(predict6(gp, nom5, [w])[0])
        rows.append({"width_mm": w, "GP_peak_MPa": g, "FE_peak_MPa": r[0], "GP_minus_FE_MPa": g - r[0],
                     "GP_vs_published_857.8_pct": 100 * (g / p.PAPER_PEAK_FE - 1), "FE_vs_published_857.8_pct": 100 * (r[0] / p.PAPER_PEAK_FE - 1),
                     "FE_stiffness_N_per_mm": r[1], "solver_s": r[2], "inside_old_FE_width_range": bool(old.width.min() <= w <= old.width.max())})
    nw = pd.DataFrame(rows); nw.to_csv(OUT / "nominal_peak_vs_width.csv", index=False); print(nw.round(3).to_string(index=False))
    wg = np.linspace(15.2, 17.5, 60)
    gp_curve = predict6(gp, np.repeat(nom5, len(wg), 0), wg)
    slope = np.polyfit(wg, gp_curve, 1)[0]
    print(f"GP nominal peak slope vs width {slope:.2f} MPa/mm")

    # ---------- 3/4. cases
    rng = np.random.default_rng(1)
    u8 = rng.random((N_MC, 8))                       # identical stream to Stage D (so W0 is comparable)
    uw = np.random.default_rng(101).random(N_MC)     # width driver, common to W1 and W2
    u9 = np.column_stack([u8, uw])
    prob = {"num_vars": 9, "names": ft.NAMES + ["width"], "bounds": [[0, 1]] * 9}
    us = sobol_sample.sample(prob, N_SOBOL, calc_second_order=False, seed=1)
    gp5 = None
    rows, sob, base = [], [], {}
    mean8 = np.array([[d.mean() for d in ft.DISTS]])
    for cname, half in CASES.items():
        x, w = build_x(u9, half)
        pk = peak_at(gp, x, w); base[cname] = pk
        uts, ratio = x[:, -2], x[:, -1]
        xs, ws = build_x(us, half); pks = peak_at(gp, xs, ws)
        for c in CRITS:
            f = ft.g_fun(pk, uts, ratio, c) < 0
            pf = f.mean()
            pff, beta, ok = form9(gp, c, half)
            mx = mean8.copy()
            L = brentq(lambda ld: float(ft.g_fun(predict6(gp, mx[:, :5], [W0]) * ld / p.RUNOUT_FORCE, mx[:, -2], mx[:, -1], c)[0]), 100, 10000)
            punit = pk / (x[:, 5] / p.RUNOUT_FORCE)
            L50 = brentq(lambda ld: np.mean(ft.g_fun(punit * ld * (x[:, 5] / p.RUNOUT_FORCE) / p.RUNOUT_FORCE, uts, ratio, c) < 0) - 0.5, 300, 5000)
            rows.append({"case": cname, "half_width_mm": half, "criterion": c, "Pf": pf, "mc_se": np.sqrt(pf * (1 - pf) / N_MC),
                         "Pf_FORM": pff, "beta": beta, "FORM_converged": ok, "runout_load_nominal_N": L,
                         "runout_vs_1200_pct": 100 * (L / p.RUNOUT_FORCE - 1), "load_at_Pf_0.5_N": L50,
                         "load_at_Pf_0.5_vs_1200_pct": 100 * (L50 / p.RUNOUT_FORCE - 1)})
            s = sobol.analyze(prob, ft.g_fun(pks, xs[:, -2], xs[:, -1], c), calc_second_order=False, seed=1)
            for i, n in enumerate(prob["names"]):
                sob.append({"case": cname, "output": f"g_{c}", "variable": n, "S1": s["S1"][i], "ST": s["ST"][i], "ST_conf": s["ST_conf"][i]})
        s = sobol.analyze(prob, pks, calc_second_order=False, seed=1)
        for i, n in enumerate(prob["names"]):
            sob.append({"case": cname, "output": "peak_stress", "variable": n, "S1": s["S1"][i], "ST": s["ST"][i], "ST_conf": s["ST_conf"][i]})
        base[cname + "_sd"] = float(pk.std())
    res = pd.DataFrame(rows); res.to_csv(OUT / "pf_form_runout_by_case.csv", index=False)
    so = pd.DataFrame(sob); so.to_csv(OUT / "sobol_by_case.csv", index=False)
    print(res.to_string(index=False, float_format=lambda v: f"{v:.5g}"))
    print(so[so.output.isin(["peak_stress", "g_goodman", "g_gerber"])].pivot_table(index=["output", "variable"], columns="case", values="ST").round(3).to_string())
    top = so[so.output.str.startswith("g_")].groupby(["case", "output"]).apply(lambda d: " > ".join(d.sort_values("ST", ascending=False).variable.head(3))).rename("top3_ST")
    top.to_csv(OUT / "sobol_top3_g.csv"); print(top.to_string())
    # paired differences vs W0
    rows = []
    xw0, ww0 = build_x(u9, 0.0)
    for c in CRITS:
        f0 = ft.g_fun(base["W0"], xw0[:, -2], xw0[:, -1], c) < 0
        for cname in ("W1", "W2"):
            x, w = build_x(u9, CASES[cname])
            f = ft.g_fun(base[cname], x[:, -2], x[:, -1], c) < 0
            rows.append({"criterion": c, "case": cname, "Pf": f.mean(), "Pf_W0": f0.mean(), "diff": f.mean() - f0.mean(),
                         "se_paired": np.std(f.astype(float) - f0.astype(float), ddof=1) / np.sqrt(N_MC)})
    pd.DataFrame(rows).to_csv(OUT / "pf_differences_vs_W0.csv", index=False); print(pd.DataFrame(rows).round(4).to_string(index=False))
    # W0 vs the Stage D 5-input GP
    df5 = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    gp5, q2_5, rmse_5 = sens.fit_surrogate(df5[FE5].to_numpy(), df5["peak_vm"].to_numpy())
    pk5 = ft.peak_stress(gp5, ft.to_physical(u8))
    xx = ft.to_physical(u8)
    cmp_rows = [{"criterion": c, "Pf_5input_GP_StageD": (ft.g_fun(pk5, xx[:, -2], xx[:, -1], c) < 0).mean(),
                 "Pf_W0_6input_GP": (ft.g_fun(base["W0"], xx[:, -2], xx[:, -1], c) < 0).mean(),
                 "mean_peak_5input": pk5.mean(), "mean_peak_W0_6input": base["W0"].mean(), "rms_peak_diff_MPa": float(np.sqrt(np.mean((pk5 - base["W0"]) ** 2)))}
                for c in CRITS]
    pd.DataFrame(cmp_rows).to_csv(OUT / "W0_vs_5input_GP.csv", index=False); print(pd.DataFrame(cmp_rows).round(5).to_string(index=False))
    # stress variance shares from the peak-stress Sobol
    json.dump({"cv": cv, "stress_sd_W0": base["W0_sd"], "stress_sd_W1": base["W1_sd"], "stress_sd_W2": base["W2_sd"],
               "gp_slope_MPa_per_mm": slope, "n_old_runs": len(old), "n_new_FE_runs": len(jobs) + 8, "kernel": str(gp.kernel_)},
              open(OUT / "meta.json", "w"), indent=2, default=str)

    # ---------- figures
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.scatter(wo, yo, s=8, color="0.6", label="FE runs (256, other inputs vary)")
    ax.plot(wg, gp_curve, color="#4c78a8", label="GP, other inputs at nominal")
    ax.plot(nw.width_mm, nw.FE_peak_MPa, "o", color="#f58518", label="FE, nominal inputs")
    ax.axhline(p.PAPER_PEAK_FE, color="k", ls="--", lw=0.8, label="published 857.8 MPa")
    ax.axvline(W0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("plate width (mm)"); ax.set_ylabel("peak von Mises at 1200 N (MPa)"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig_peak_vs_width.png", dpi=600); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    wd = 0.35
    for j, c in enumerate(("goodman", "gerber")):
        d = res[res.criterion == c]
        ax.bar(np.arange(3) + (j - 0.5) * wd, d.Pf, wd, yerr=1.96 * d.mc_se, capsize=3, label=c.capitalize())
    ax.set_xticks(range(3)); ax.set_xticklabels(["W0 (fixed)", "W1 (+-0.5 mm, ASSUMED)", "W2 (+-1.0 mm, ASSUMED)"], fontsize=8)
    ax.set_ylabel("P_f at 1200 N (95 % MC interval)"); ax.legend()
    fig.tight_layout(); fig.savefig(OUT / "fig_pf_by_width_case.png", dpi=600); plt.close(fig)


if __name__ == "__main__":
    main()
