"""Stage F: surrogate benchmark for the peak von Mises stress, and its effect on the failure probability.

Data: results/stageA/fe_samples_al.csv (416 FE runs: 256 Sobol DOE + 160 active-learning boundary points), 5 FE inputs.
Split: random 80/20, stratified by source (DOE / active learning), seed 0. All tuning uses 5-fold CV on the training
part only; the 20 % test set is touched once for the scores.
Models: GP (as sensitivity.py), neural network (sklearn MLP), gradient boosting (sklearn), polynomial chaos
(OpenTURNS, Hermite/Legendre basis on the input distributions, least squares or LARS, degree chosen by CV).
P_f: each model replaces the surrogate in the Stage D "assumed" limit state (Goodman, Gerber, SWT; Kf = Kt), same
400,000 MC inputs for every model. GP uncertainty: MC error + epistemic error from joint posterior draws of the
LATENT GP (white-noise kernel excluded) on a 5,000-point subset of the MC pool.
Usage: python stage_f.py
"""
import json
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openturns as ot
import pandas as pd
from scipy import stats
from scipy.linalg import solve_triangular
from sklearn.compose import TransformedTargetRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GridSearchCV, KFold, train_test_split
from sklearn.neural_network import MLPRegressor

import fatigue as ft
import params as p
import sensitivity as sens
import uncertainty as unc

warnings.filterwarnings("ignore")
OUT = p.RESULTS_DIR / "stageF"
N_MC, M_SUB, S_DRAWS = 400_000, 5000, 300
CRITERIA = ("goodman", "gerber", "swt")
NAMES = unc.FE_NAMES


def scale(x):
    return (x - sens.MU) / sens.SD


# ---------------- models ----------------
def fit_gp(x, y):
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
    k = ConstantKernel(1.0) * RBF(np.ones(x.shape[1]), (1e-2, 1e3)) + WhiteKernel(1e-4, (1e-8, 1e-1))
    return GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=3, random_state=0).fit(scale(x), y)


def fit_nn(x, y):
    base = TransformedTargetRegressor(MLPRegressor(solver="lbfgs", activation="tanh", max_iter=5000, random_state=0),
                                      transformer=__import__("sklearn.preprocessing", fromlist=["x"]).StandardScaler())
    grid = {"regressor__hidden_layer_sizes": [(16,), (32, 32), (64, 64)], "regressor__alpha": [1e-4, 1e-2, 1.0]}
    gs = GridSearchCV(base, grid, cv=KFold(5, shuffle=True, random_state=0), scoring="neg_root_mean_squared_error").fit(scale(x), y)
    return gs


def fit_gbm(x, y):
    grid = {"n_estimators": [200, 600], "max_depth": [2, 3], "learning_rate": [0.05, 0.1], "subsample": [1.0, 0.7]}
    gs = GridSearchCV(GradientBoostingRegressor(random_state=0), grid, cv=KFold(5, shuffle=True, random_state=0),
                      scoring="neg_root_mean_squared_error").fit(scale(x), y)
    return gs


def ot_input():
    marg = []
    for k in NAMES:
        d = unc.FE_INPUTS[k]
        if d.dist.name == "norm":
            marg.append(ot.Normal(d.mean(), d.std()))
        else:
            marg.append(ot.Uniform(float(d.ppf(0)), float(d.ppf(1))))
    return marg


def fit_pce_cfg(x, y, degree, sparse):
    marg = ot_input()
    dist = ot.JointDistribution(marg)
    polys = [ot.StandardDistributionPolynomialFactory(m) for m in marg]
    enum = ot.LinearEnumerateFunction(len(marg))
    basis = ot.OrthogonalProductPolynomialFactory(polys, enum)
    size = enum.getStrataCumulatedCardinal(degree)
    if sparse:
        strat = ot.LeastSquaresMetaModelSelectionFactory(ot.LARS(), ot.CorrectedLeaveOneOut())
        adaptive = ot.FixedStrategy(basis, size)
        proj = ot.LeastSquaresStrategy(strat)
    else:
        adaptive = ot.FixedStrategy(basis, size)
        proj = ot.LeastSquaresStrategy()
    algo = ot.FunctionalChaosAlgorithm(ot.Sample(x), ot.Sample(y.reshape(-1, 1)), dist, adaptive, proj)
    algo.run()
    mm = algo.getResult().getMetaModel()
    return lambda xn: np.array(mm(ot.Sample(xn))).ravel()


def fit_pce(x, y):
    best, rows = None, []
    for sparse, degs in ((False, (2, 3, 4)), (True, (3, 4, 5, 6))):
        for d in degs:
            errs = []
            for tr, va in KFold(5, shuffle=True, random_state=0).split(x):
                f = fit_pce_cfg(x[tr], y[tr], d, sparse)
                errs.append(np.sqrt(np.mean((f(x[va]) - y[va]) ** 2)))
            rows.append({"sparse_LARS": sparse, "degree": d, "cv_rmse": float(np.mean(errs))})
    cv = pd.DataFrame(rows)
    b = cv.loc[cv.cv_rmse.idxmin()]
    return fit_pce_cfg(x, y, int(b.degree), bool(b.sparse_LARS)), cv, b


# ---------------- GP latent posterior ----------------
def gp_latent(gp, xs_fe):
    """Mean and covariance (MPa) of the latent GP at points xs_fe (physical FE inputs), white-noise kernel excluded."""
    k1 = gp.kernel_.k1
    xs = scale(xs_fe)
    ks = k1(xs, gp.X_train_)
    mean = ks @ gp.alpha_ * gp._y_train_std + gp._y_train_mean
    v = solve_triangular(gp.L_, ks.T, lower=True)
    cov = (k1(xs) - v.T @ v) * gp._y_train_std**2
    return mean, cov


def pf_of(peak_model, force, uts, ratio, crit):
    g = ft.g_fun(peak_model * force / p.RUNOUT_FORCE, uts, ratio, crit)
    v = np.mean(g < 0)
    return v, np.sqrt(v * (1 - v) / len(g))


def main():
    OUT.mkdir(exist_ok=True)
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    x, y = df[NAMES].to_numpy(), df["peak_vm"].to_numpy()
    src = np.r_[np.zeros(256), np.ones(len(df) - 256)]
    idx = np.arange(len(df))
    tr, te = train_test_split(idx, test_size=0.2, random_state=0, stratify=src)
    print(f"train {len(tr)}, test {len(te)} (DOE {int((src[te] == 0).sum())}, AL {int((src[te] == 1).sum())})")
    xt, yt, xe, ye = x[tr], y[tr], x[te], y[te]

    models, info = {}, {}
    gp = fit_gp(xt, yt); models["GP"] = lambda a, g=gp: g.predict(scale(a))
    info["GP"] = {"kernel": str(gp.kernel_)}
    nn = fit_nn(xt, yt); models["NeuralNet"] = lambda a, g=nn: g.predict(scale(a))
    info["NeuralNet"] = {"best": {k: str(v) for k, v in nn.best_params_.items()}, "cv_rmse": -nn.best_score_}
    gb = fit_gbm(xt, yt); models["GradBoost"] = lambda a, g=gb: g.predict(scale(a))
    info["GradBoost"] = {"best": {k: str(v) for k, v in gb.best_params_.items()}, "cv_rmse": -gb.best_score_}
    pce, cv, b = fit_pce(xt, yt); models["PCE"] = pce
    info["PCE"] = {"best": {"degree": int(b.degree), "sparse_LARS": bool(b.sparse_LARS)}, "cv_rmse": float(b.cv_rmse)}
    cv.to_csv(OUT / "pce_cv.csv", index=False)

    # ---- test-set scores ----
    rows, preds = [], {}
    for name, f in models.items():
        pr = f(xe); preds[name] = pr
        for sub, m in (("all", np.ones(len(te), bool)), ("DOE points", src[te] == 0), ("active-learning points", src[te] == 1)):
            e = pr[m] - ye[m]
            rows.append({"model": name, "test_subset": sub, "n": int(m.sum()),
                         "R2": 1 - np.sum(e**2) / np.sum((ye[m] - ye[m].mean()) ** 2), "RMSE_MPa": np.sqrt(np.mean(e**2)),
                         "MAE_MPa": np.mean(np.abs(e)), "max_abs_MPa": np.abs(e).max(),
                         "max_abs_pct_of_peak": 100 * np.max(np.abs(e) / ye[m]),
                         "cv_rmse_train": info[name].get("cv_rmse", np.nan)})
    sc = pd.DataFrame(rows); sc.to_csv(OUT / "test_scores.csv", index=False)
    pd.DataFrame({"test_index": te, "fe_peak": ye, **{k: v for k, v in preds.items()}}).to_csv(OUT / "test_predictions.csv", index=False)
    # paired bootstrap (2000 resamples of the 84 test points): RMSE difference vs GP, 95 % interval
    brng, bs = np.random.default_rng(0), []
    ids = brng.integers(0, len(te), (2000, len(te)))
    rm = lambda pr: np.sqrt(np.mean((pr[ids] - ye[ids]) ** 2, axis=1))
    for name in preds:
        if name != "GP":
            d = rm(preds[name]) - rm(preds["GP"])
            bs.append({"model": name, "rmse_minus_GP_mean": d.mean(), "ci95_low": np.percentile(d, 2.5), "ci95_high": np.percentile(d, 97.5)})
    pd.DataFrame(bs).to_csv(OUT / "test_rmse_vs_gp_bootstrap.csv", index=False); print(pd.DataFrame(bs).round(3).to_string(index=False))
    print(sc.round(4).to_string(index=False))

    # ---- failure probability with each surrogate (common random numbers) ----
    rng = np.random.default_rng(1)
    xp = ft.to_physical(rng.random((N_MC, len(ft.NAMES))))
    xfe, force, uts, ratio = xp[:, :len(NAMES)], xp[:, len(NAMES)], xp[:, -2], xp[:, -1]
    rows = []
    peak_by = {}
    for name, f in models.items():
        pk = np.concatenate([f(xfe[i:i + 100_000]) for i in range(0, N_MC, 100_000)]); peak_by[name] = pk
        for crit in CRITERIA:
            v, se = pf_of(pk, force, uts, ratio, crit)
            rows.append({"model": name, "criterion": crit, "Pf": v, "mc_se": se, "peak_mean_MPa": pk.mean()})
    pfdf = pd.DataFrame(rows); pfdf.to_csv(OUT / "pf_by_model.csv", index=False)
    print(pfdf.pivot(index="model", columns="criterion", values="Pf").round(4).to_string())

    # ---- GP uncertainty in P_f (train-only GP and GP refitted on all 416 runs) ----
    gp_all = fit_gp(x, y)
    rows = []
    sub = slice(0, M_SUB)
    for lab, g in (("GP (train 332)", gp), ("GP (all 416)", gp_all)):
        mean, cov = gp_latent(g, xfe[sub])
        assert np.allclose(mean, g.predict(scale(xfe[sub])), atol=1e-6), "latent mean mismatch"
        L = np.linalg.cholesky(cov + 1e-6 * np.eye(M_SUB))
        draws = mean[:, None] + L @ np.random.default_rng(2).standard_normal((M_SUB, S_DRAWS))
        pm, sd_tot = g.predict(scale(xfe), return_std=True)
        for crit in CRITERIA:
            v, se = pf_of(pm, force, uts, ratio, crit)
            pf_draw = np.array([np.mean(ft.g_fun(draws[:, s] * force[sub] / p.RUNOUT_FORCE, uts[sub], ratio[sub], crit) < 0)
                                for s in range(S_DRAWS)])
            pf_sub_plug = np.mean(ft.g_fun(mean * force[sub] / p.RUNOUT_FORCE, uts[sub], ratio[sub], crit) < 0)
            # per-point probabilistic classification with the total predictive sd (includes FE noise, 2-sided Phi)
            gm = ft.g_fun(pm * force / p.RUNOUT_FORCE, uts, ratio, crit)
            slope = ft.g_fun(pm * force / p.RUNOUT_FORCE + 1.0, uts, ratio, crit) - gm
            pf_prob = np.mean(stats.norm.cdf(-gm / np.maximum(np.abs(slope) * sd_tot * force / p.RUNOUT_FORCE, 1e-9)))
            epi = pf_draw.std(ddof=1)
            tot = np.sqrt(se**2 + epi**2)
            rows.append({"gp": lab, "criterion": crit, "Pf_plugin": v, "mc_se": se, "epistemic_sd_latent": epi,
                         "ci95_low": v - 1.96 * tot, "ci95_high": v + 1.96 * tot,
                         "draw_mean_minus_plugin_on_subset": pf_draw.mean() - pf_sub_plug,
                         "Pf_pointwise_prob_total_sd": pf_prob, "mean_pred_sd_MPa": sd_tot.mean(),
                         "mean_latent_sd_MPa": float(np.sqrt(np.diag(cov)).mean())})
    gu = pd.DataFrame(rows); gu.to_csv(OUT / "gp_pf_uncertainty.csv", index=False)
    print(gu.round(4).to_string(index=False))
    json.dump(info, open(OUT / "model_selection.json", "w"), indent=2, default=str)

    fig, ax = plt.subplots(1, 5, figsize=(20, 4))
    for a, (name, pr) in zip(ax, preds.items()):
        a.scatter(ye, pr, s=12, c=src[te], cmap="coolwarm"); lim = [ye.min(), ye.max()]
        a.plot(lim, lim, "k--", lw=1); a.set_title(f"{name}: RMSE {np.sqrt(np.mean((pr - ye) ** 2)):.2f} MPa")
        a.set_xlabel("FE peak (MPa)"); a.set_ylabel("surrogate (MPa)")
    w = 0.2
    for i, crit in enumerate(CRITERIA):
        d = pfdf[pfdf.criterion == crit]
        ax[4].bar(np.arange(len(d)) + (i - 1) * w, d.Pf, w, yerr=1.96 * d.mc_se, label=crit)
    ax[4].set_xticks(range(len(d))); ax[4].set_xticklabels(d.model, rotation=20); ax[4].set_ylabel("P_f"); ax[4].legend()
    fig.tight_layout(); fig.savefig(OUT / "benchmark.png", dpi=150)


if __name__ == "__main__":
    main()
