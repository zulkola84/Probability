"""Phase 2, step 3: closed-form baseline vs the Gaussian-process surrogate for the peak stress.

Closed form: sigma_pk = Kt_net * sigma_nom,  sigma_nom = 6 M / ((W - d) t^2),  M = F h / 2
(F total force over both loading rollers, h = lever arm, W = 16.365 mm CALIBRATED width, never sampled).
E and nu do not enter the closed form.
Kt_net versions:
  (a) fixed = Kt_rect_net from results/stress_field/summary.json (FE, rectangle, net section);
  (b) Howland / Heywood finite-width fit for a circular hole in IN-PLANE TENSION (Peterson, Stress Concentration
      Factors, Chart 4.1): Kt_net = 3.00 - 3.13 (d/W) + 3.66 (d/W)^2 - 1.53 (d/W)^3. ASSUMED: the source for
      a bending-specific formula for this out-of-plane-bending geometry could not be confirmed, and the formula is
      quoted from memory (not checked against the book). It is a different loading case.
Comparison: the 84-point hold-out set of the Stage F benchmark (same split), GP trained on the other 332 runs.
Limit state, input distributions and the 400,000 common-random-number MC samples as in Stage D ("assumed" case,
Kf = Kt): P_f (Goodman, Gerber, SWT), FORM, total-order Sobol indices (base sample 8192).
All outputs go to OUT (outside the project's results/ folder).
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SALib.analyze import sobol
from SALib.sample import sobol as sobol_sample
from scipy import stats
from sklearn.model_selection import train_test_split

import fatigue as ft
import params as p
import sensitivity as sens
import stage_f
import uncertainty as unc

OUT = Path(__file__).resolve().parents[1] / "results" / "phase2" / "step3_baseline"
N_MC, N_SOBOL = 400_000, 2**13
CRITERIA = ("goodman", "gerber", "swt")
NAMES = ft.NAMES                          # thickness, hole_d, E, nu, lever_arm, force, uts, endurance_ratio
FE_NAMES = unc.FE_NAMES
IT, ID, IE, INU, IH, IF = range(6)
W = p.WIDTH


def sigma_nom(t, d, h, f):
    return 6 * (f * h / 2) / ((W - d) * t**2)


def kt_howland(d):
    r = d / W
    return 3.00 - 3.13 * r + 3.66 * r**2 - 1.53 * r**3


def make_models(gp, kt_fixed):
    def gp_peak(x):
        return ft.peak_stress(gp, x)

    def cf_a(x):
        return kt_fixed * sigma_nom(x[:, IT], x[:, ID], x[:, IH], x[:, IF])

    def cf_b(x):
        return kt_howland(x[:, ID]) * sigma_nom(x[:, IT], x[:, ID], x[:, IH], x[:, IF])
    return {"GP": gp_peak, "closed_form_a_Kt_FE": cf_a, "closed_form_b_Kt_Howland": cf_b}


def scores(pred, fe, mask=None):
    mask = np.ones(len(fe), bool) if mask is None else mask
    e, y = pred[mask] - fe[mask], fe[mask]
    return {"n": int(mask.sum()), "RMSE_MPa": np.sqrt(np.mean(e**2)), "bias_MPa": e.mean(),
            "max_abs_err_MPa": np.abs(e).max(), "max_abs_err_pct_of_peak": 100 * np.max(np.abs(e) / y),
            "Q2": 1 - np.sum(e**2) / np.sum((y - y.mean()) ** 2)}


def form(peak_fn, crit, iters=100):
    """HL-RF in standard-normal space (same algorithm and convention as fatigue.form). Returns Pf, beta, converged."""
    def g_u(z):
        x = ft.to_physical(stats.norm.cdf(z)[None, :])
        return float(ft.g_fun(peak_fn(x), x[:, -2], x[:, -1], crit)[0])
    z, done = np.zeros(len(NAMES)), False
    for _ in range(iters):
        g0 = g_u(z)
        grad = np.array([(g_u(z + 1e-4 * e) - g0) / 1e-4 for e in np.eye(len(z))])
        z_new = (grad @ z - g0) / (grad @ grad) * grad
        done = np.linalg.norm(z_new - z) < 1e-6
        z = z_new
        if done:
            break
    beta = np.linalg.norm(z)
    if g_u(np.zeros(len(NAMES))) < 0:       # nominal plate fails: origin in the failure region (beta reported negative)
        return stats.norm.cdf(beta), -beta, done
    return stats.norm.cdf(-beta), beta, done


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summ = json.load(open(p.RESULTS_DIR / "stress_field" / "summary.json"))
    kt_fixed = summ["Kt_rect_net"]

    # ---- step 1: nominal check ----
    sn0 = float(sigma_nom(p.THICKNESS, p.HOLE_DIAMETER, p.LEVER_ARM, p.RUNOUT_FORCE))
    print(f"nominal sigma_nom = {sn0:.2f} MPa (summary.json nominal_rect_net_MPa {summ['nominal_rect_net_MPa']:.2f})")
    if abs(sn0 - 609) > 0.01 * 609:
        raise SystemExit(f"STOP: sigma_nom {sn0:.1f} MPa is not about 609 MPa")
    kt_b0 = float(kt_howland(p.HOLE_DIAMETER))
    nominal = {"sigma_nom_MPa": sn0, "Kt_a_fixed_FE": kt_fixed, "Kt_b_Howland_at_d7.8": kt_b0,
               "peak_a_MPa": kt_fixed * sn0, "peak_b_MPa": kt_b0 * sn0, "FE_peak_vm_stress_field_MPa": summ["overall_peak_von_mises_MPa"],
               "d_over_W": p.HOLE_DIAMETER / W}
    print(nominal)

    # ---- data, split, GPs ----
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    x, y = df[FE_NAMES].to_numpy(), df["peak_vm"].to_numpy()
    src = np.r_[np.zeros(256), np.ones(len(df) - 256)]
    tr, te = train_test_split(np.arange(len(df)), test_size=0.2, random_state=0, stratify=src)
    gp_tr = stage_f.fit_gp(x[tr], y[tr])
    gp_all, q2_all, rmse_all = sens.fit_surrogate(x, y)      # as Stage D (used for P_f, FORM, Sobol)
    pred_te = {"GP": gp_tr.predict(stage_f.scale(x[te]))}
    # closed forms on the test inputs (force fixed at 1200 N, the FE runs' reference load)
    xt = np.zeros((len(te), 6))
    xt[:, :5] = x[te]; xt[:, IF] = p.RUNOUT_FORCE
    mods_te = make_models(gp_all, kt_fixed)
    pred_te["closed_form_a_Kt_FE"] = mods_te["closed_form_a_Kt_FE"](xt)
    pred_te["closed_form_b_Kt_Howland"] = mods_te["closed_form_b_Kt_Howland"](xt)
    ye = y[te]
    rows = []
    for name, pr in pred_te.items():
        for sub, m in (("all", None), ("DOE points", src[te] == 0), ("active-learning points", src[te] == 1)):
            rows.append({"model": name, "test_subset": sub, **scores(pr, ye, m)})
    sc = pd.DataFrame(rows); sc.to_csv(OUT / "holdout_scores.csv", index=False)
    print(sc.round(4).to_string(index=False))
    gp_row = sc[(sc.model == "GP") & (sc.test_subset == "all")].iloc[0]
    print(f"GP check vs Stage F: RMSE {gp_row.RMSE_MPa:.3f} (Stage F 2.89), Q2 {gp_row.Q2:.4f} (0.9735)")
    pd.DataFrame({"row_in_fe_samples_al": te, "source": np.where(src[te] == 0, "DOE", "AL"), "fe_peak_MPa": ye,
                  **{k: v for k, v in pred_te.items()}, **dict(zip(FE_NAMES, x[te].T))}).to_csv(OUT / "holdout_predictions.csv", index=False)

    # residual correlation with inputs (where does the closed form fail?)
    rr = []
    for name in ("closed_form_a_Kt_FE", "closed_form_b_Kt_Howland", "GP"):
        res = pred_te[name] - ye
        for k, col in zip(FE_NAMES, x[te].T):
            rr.append({"model": name, "input": k, "pearson_r_residual_vs_input": np.corrcoef(res, col)[0, 1]})
        rr.append({"model": name, "input": "fe_peak", "pearson_r_residual_vs_input": np.corrcoef(res, ye)[0, 1]})
    pd.DataFrame(rr).to_csv(OUT / "holdout_residual_correlation.csv", index=False)

    # local log-elasticities d ln(peak)/d ln(x) at the nominal point (central difference, +-1 %)
    mods = make_models(gp_all, kt_fixed)
    nom = np.array([[d.mean() for d in ft.DISTS[:6]]])
    el = []
    for i, k in enumerate(NAMES[:6]):
        row = {"input": k}
        for name, f in mods.items():
            xp_, xm_ = nom.copy(), nom.copy(); xp_[0, i] *= 1.01; xm_[0, i] *= 0.99
            row[name] = float((np.log(f(xp_)[0]) - np.log(f(xm_)[0])) / (np.log(1.01) - np.log(0.99)))
        el.append(row)
    el = pd.DataFrame(el); el.to_csv(OUT / "nominal_elasticities.csv", index=False); print(el.round(3).to_string(index=False))

    # ---- Monte Carlo (common random numbers; same seed/construction as Stage D) ----
    rng = np.random.default_rng(1)
    xp = ft.to_physical(rng.random((N_MC, len(NAMES))))
    uts, ratio = xp[:, -2], xp[:, -1]
    peaks = {k: np.concatenate([f(xp[i:i + 100_000]) for i in range(0, N_MC, 100_000)]) for k, f in mods.items()}
    ind, rows = {}, []
    for name, pk in peaks.items():
        for crit in CRITERIA:
            f = ft.g_fun(pk, uts, ratio, crit) < 0
            ind[(name, crit)] = f
            v = f.mean()
            rows.append({"method": name, "criterion": crit, "Pf": v, "mc_se": np.sqrt(v * (1 - v) / N_MC),
                         "peak_mean_MPa": pk.mean(), "peak_sd_MPa": pk.std()})
    pf = pd.DataFrame(rows); pf.to_csv(OUT / "pf_by_method.csv", index=False)
    print(pf.round(4).to_string(index=False))

    rows = []
    for name in mods:
        if name == "GP":
            continue
        for crit in CRITERIA:
            a, b = ind[(name, crit)], ind[("GP", crit)]
            d = a.mean() - b.mean()
            se_ind = np.sqrt(a.mean() * (1 - a.mean()) / N_MC + b.mean() * (1 - b.mean()) / N_MC)
            se_pair = np.std(a.astype(float) - b.astype(float), ddof=1) / np.sqrt(N_MC)
            rows.append({"method": name, "criterion": crit, "Pf_method": a.mean(), "Pf_GP": b.mean(), "diff": d,
                         "se_independent": se_ind, "se_paired_CRN": se_pair, "abs_diff_over_se_independent": abs(d) / se_ind,
                         "within_1.96se_independent": bool(abs(d) <= 1.96 * se_ind),
                         "within_1.96se_paired": bool(abs(d) <= 1.96 * se_pair),
                         "n_classified_differently": int((a != b).sum()),
                         "mean_peak_diff_MPa": peaks[name].mean() - peaks["GP"].mean(),
                         "rms_peak_diff_MPa": float(np.sqrt(np.mean((peaks[name] - peaks["GP"]) ** 2)))})
    pdiff = pd.DataFrame(rows); pdiff.to_csv(OUT / "pf_differences.csv", index=False); print(pdiff.round(4).to_string(index=False))

    # DIAGNOSTIC only (not a model): closed form (a) with its mean offset to the GP over the MC pool removed, to see
    # whether the P_f difference is a pure mean shift or a shape (input-dependence) difference.
    off = (peaks["closed_form_a_Kt_FE"] - peaks["GP"]).mean()
    rows = []
    for crit in CRITERIA:
        a = ft.g_fun(peaks["closed_form_a_Kt_FE"] - off, uts, ratio, crit) < 0
        b = ind[("GP", crit)]
        d = a.mean() - b.mean()
        se_ind = np.sqrt(a.mean() * (1 - a.mean()) / N_MC + b.mean() * (1 - b.mean()) / N_MC)
        rows.append({"criterion": crit, "offset_removed_MPa": off, "Pf_cf_a_shifted": a.mean(), "Pf_GP": b.mean(), "diff": d,
                     "se_independent": se_ind, "within_1.96se_independent": bool(abs(d) <= 1.96 * se_ind)})
    pd.DataFrame(rows).to_csv(OUT / "diagnostic_cf_a_offset_removed.csv", index=False)
    print(pd.DataFrame(rows).round(4).to_string(index=False))

    # ---- FORM ----
    rows = []
    for name, f in mods.items():
        for crit in CRITERIA:
            pff, beta, ok = form(f, crit)
            rows.append({"method": name, "criterion": crit, "Pf_FORM": pff, "beta": beta, "converged": ok})
    fm = pd.DataFrame(rows); fm.to_csv(OUT / "form_by_method.csv", index=False); print(fm.round(4).to_string(index=False))

    # ---- Sobol total order of g ----
    prob = {"num_vars": len(NAMES), "names": NAMES, "bounds": [[0, 1]] * len(NAMES)}
    us = sobol_sample.sample(prob, N_SOBOL, calc_second_order=False, seed=1)
    xs = ft.to_physical(us)
    rows = []
    for name, f in mods.items():
        pk = f(xs)
        for crit in CRITERIA:
            s = sobol.analyze(prob, ft.g_fun(pk, xs[:, -2], xs[:, -1], crit), calc_second_order=False, seed=1)
            for i, n in enumerate(NAMES):
                rows.append({"method": name, "criterion": crit, "variable": n, "S1": s["S1"][i], "ST": s["ST"][i], "ST_conf": s["ST_conf"][i]})
    so = pd.DataFrame(rows); so.to_csv(OUT / "sobol_by_method.csv", index=False)
    gpst = so[so.method == "GP"].set_index(["criterion", "variable"])
    drows = []
    for name in mods:
        if name == "GP":
            continue
        t = so[so.method == name].set_index(["criterion", "variable"])
        for (crit, v), r in t.iterrows():
            g = gpst.loc[(crit, v)]
            drows.append({"method": name, "criterion": crit, "variable": v, "ST_method": r.ST, "ST_GP": g.ST,
                          "diff": r.ST - g.ST, "conf95_combined": float(np.hypot(r.ST_conf, g.ST_conf))})
    sd = pd.DataFrame(drows); sd.to_csv(OUT / "sobol_differences.csv", index=False)
    print(so.pivot_table(index=["criterion", "variable"], columns="method", values="ST").round(3).to_string())
    rank = so.groupby(["method", "criterion"]).apply(lambda d: " > ".join(d.sort_values("ST", ascending=False).variable.head(3))).rename("top3_ST")
    rank.to_csv(OUT / "sobol_top3.csv"); print(rank.to_string())

    # ---- figure ----
    cols = {"GP": "#222222", "closed_form_a_Kt_FE": "#4c78a8", "closed_form_b_Kt_Howland": "#f58518"}
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    for name, pr in pred_te.items():
        ax[0].scatter(ye, pr, s=16, color=cols[name], alpha=.8,
                      label=f"{name} (RMSE {scores(pr, ye)['RMSE_MPa']:.1f} MPa)")
    lim = [ye.min() - 20, ye.max() + 20]
    ax[0].plot(lim, lim, "k--", lw=1); ax[0].set_xlabel("FE peak von Mises at 1200 N (MPa), 84 hold-out points")
    ax[0].set_ylabel("predicted peak stress (MPa)"); ax[0].legend(fontsize=8)
    w = 0.26
    for j, (name, c) in enumerate(cols.items()):
        d = pf[pf.method == name].set_index("criterion").loc[list(CRITERIA)]
        ax[1].bar(np.arange(3) + (j - 1) * w, d.Pf, w, yerr=1.96 * d.mc_se, color=c, label=name, capsize=2)
    ax[1].set_xticks(range(3)); ax[1].set_xticklabels([c.capitalize() if c != "swt" else "SWT" for c in CRITERIA])
    ax[1].set_ylabel("P_f at 1200 N (95 % MC interval)"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "comparison.png", dpi=600)

    json.dump({"nominal": nominal, "gp_all_q2": q2_all, "gp_all_cv_rmse": rmse_all, "n_mc": N_MC, "n_sobol_base": N_SOBOL,
               "note": "Kt version (b) is ASSUMED (in-plane tension formula, source not confirmed for this loading)."},
              open(OUT / "meta.json", "w"), indent=2)


if __name__ == "__main__":
    main()
