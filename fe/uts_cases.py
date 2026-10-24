"""Phase 2, step 5: tensile strength 900 vs 991 MPa.

Cases (GP on all 416 runs, Stage D "assumed" scatter for everything else, Kf = Kt = 1.66, 400,000 CRN samples, seed 1, as Stage D):
  A  UTS mean 900 MPa, CV 3 %           (control: must reproduce Stage D 0.9800 / 0.3365 / 0.9932)
  B  UTS mean 991 MPa, CV 3 %
  C  UTS mean 945.5 MPa (midpoint), CV 5 % (ASSUMED, stands for not knowing which value applies)
Endurance limit = ratio * UTS / Kf with ratio ~ N(0.5, 0.05) (ASSUMED scatter, mean 0.5 from the paper). UTS is lognormal.
Load-line runout amplitude sa* and the V2' life line from goodman_gap_v2prime.py (imported unchanged).
Outputs go to OUT (outside results/).
"""
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

import fatigue as ft
import goodman_gap_v2prime as v2
import params as p
import sensitivity as sens
import stage_e as se
import uncertainty as unc

OUT = Path(__file__).resolve().parents[1] / "results" / "phase2" / "step5_uts991"
NMC, N_SOBOL = 400_000, 2**13
CRITS = ("goodman", "gerber", "swt")
CASES = {"A_900_cv3": (900.0, 0.03), "B_991_cv3": (991.0, 0.03), "C_945.5_cv5_ASSUMED": (945.5, 0.05)}
N_FE = len(unc.FE_NAMES)


def physical(u, mean, cv):
    """Stage D physical samples with the UTS replaced by lognormal(mean, cv) driven by the same uniform column."""
    x = ft.to_physical(u).copy()
    x[:, -2] = ft.uts_dist(cv, mean).ppf(np.clip(u[:, -2], 1e-9, 1 - 1e-9))
    return x


def form(gp, crit, mean, cv, iters=100):
    def g_u(z):
        x = physical(stats.norm.cdf(z)[None, :], mean, cv)
        return float(ft.g_fun(ft.peak_stress(gp, x), x[:, -2], x[:, -1], crit)[0])
    z, done = np.zeros(len(ft.NAMES)), False
    for _ in range(iters):
        g0 = g_u(z)
        grad = np.array([(g_u(z + 1e-4 * e) - g0) / 1e-4 for e in np.eye(len(z))])
        z_new = (grad @ z - g0) / (grad @ grad) * grad
        done = np.linalg.norm(z_new - z) < 1e-6
        z = z_new
        if done:
            break
    beta = np.linalg.norm(z)
    if g_u(np.zeros(len(z))) < 0:
        return stats.norm.cdf(beta), -beta, done
    return stats.norm.cdf(-beta), beta, done


def runout_load(gp, crit, uts_mean):
    mean = np.array([[d.mean() for d in ft.DISTS]]); mean[0, -2] = uts_mean

    def f(load):
        x = mean.copy(); x[0, N_FE] = load
        return float(ft.g_fun(ft.peak_stress(gp, x), x[:, -2], x[:, -1], crit)[0])
    return brentq(f, 100, 10000)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    gp, q2, rmse = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
    print(f"GP on {len(df)} runs: Q2 {q2:.4f}, CV RMSE {rmse:.2f} MPa")

    # ---------- 3. which UTS reproduces Table 2 ----------
    print("\n[3] Table 2 amplitudes (198.14 Goodman, 241.84 Gerber; 1e3 anchor 364.50) and the UTS they imply (ratio 0.5, Kf 1.66)")
    rows3 = []
    for c, t2 in (("goodman", 198.14), ("gerber", 241.84)):
        uts_imp = brentq(lambda u: v2.sa_star(c, u, 0.5) - t2, 300, 2000, xtol=1e-10)
        r991 = brentq(lambda r: v2.sa_star(c, 991.0, r) - t2, 0.05, 1.5, xtol=1e-12)
        rows3.append({"criterion": c, "table2_1e6_MPa": t2, "sa_star_at_UTS900": float(v2.sa_star(c, 900.0, 0.5)),
                      "sa_star_at_UTS991": float(v2.sa_star(c, 991.0, 0.5)), "UTS_implied_by_table2": uts_imp,
                      "ratio_needed_at_UTS991_to_give_table2": r991})
    uts_anchor = 364.50 / (0.9 * 0.45)
    t3 = pd.DataFrame(rows3); t3["UTS_implied_by_1e3_anchor_364.50"] = uts_anchor
    t3.to_csv(OUT / "table2_vs_uts.csv", index=False); print(t3.T.to_string())
    print("  sa* SWT at 900 / 991:", float(v2.sa_star("swt", 900.0, 0.5)), float(v2.sa_star("swt", 991.0, 0.5)))

    # ---------- 1/2. MC, FORM, runout, Sobol, life ----------
    rng = np.random.default_rng(1)
    u = rng.random((NMC, len(ft.NAMES)))
    prob = {"num_vars": len(ft.NAMES), "names": ft.NAMES, "bounds": [[0, 1]] * len(ft.NAMES)}
    us = sobol_sample.sample(prob, N_SOBOL, calc_second_order=False, seed=1)
    rows, sob, life = [], [], []
    for cname, (mean, cv) in CASES.items():
        x = physical(u, mean, cv)
        peak = ft.peak_stress(gp, x)
        uts, ratio = x[:, -2], x[:, -1]
        sa = (1 - ft.R) / 2 * peak / ft.KT_B
        xs = physical(us, mean, cv); pks = ft.peak_stress(gp, xs)
        for c in CRITS:
            f = ft.g_fun(peak, uts, ratio, c) < 0
            pf = f.mean()
            pff, beta, ok = form(gp, c, mean, cv)
            L = runout_load(gp, c, mean)
            rows.append({"case": cname, "uts_mean": mean, "uts_cv": cv, "criterion": c, "Pf": pf, "mc_se": np.sqrt(pf * (1 - pf) / NMC),
                         "Pf_FORM": pff, "beta": beta, "FORM_converged": ok, "runout_load_N": L,
                         "runout_vs_1200_pct": 100 * (L / p.RUNOUT_FORCE - 1),
                         "sa_star_nominal_at_case_mean": float(v2.sa_star(c, mean, 0.5)),
                         "sa_star_nominal_UTS900": float(v2.sa_star(c, 900.0, 0.5)),
                         "sa_star_nominal_UTS991": float(v2.sa_star(c, 991.0, 0.5))})
            s = sobol.analyze(prob, ft.g_fun(pks, xs[:, -2], xs[:, -1], c), calc_second_order=False, seed=1)
            for i, n in enumerate(ft.NAMES):
                sob.append({"case": cname, "criterion": c, "variable": n, "S1": s["S1"][i], "ST": s["ST"][i], "ST_conf": s["ST_conf"][i]})
            if c != "swt":
                b = se.TABLE2[c][1]
                n2 = 1e6 * (sa / v2.sa_star(c, uts, ratio)) ** (1 / b)
                q = np.percentile(n2, [5, 50, 95])
                life.append({"case": cname, "criterion": c, "V2prime_N_p05": q[0], "V2prime_N_median": q[1], "V2prime_N_p95": q[2],
                             "P_N_below_1e6": np.mean(n2 < 1e6), "disagreements_with_g<0": int(((n2 < 1e6) != f).sum())})
    t1 = pd.DataFrame(rows); t1.to_csv(OUT / "pf_form_runout.csv", index=False)
    so = pd.DataFrame(sob); so.to_csv(OUT / "sobol_by_case.csv", index=False)
    tl = pd.DataFrame(life); tl.to_csv(OUT / "life_v2prime_by_case.csv", index=False)
    top = so.groupby(["case", "criterion"]).apply(lambda d: " > ".join(d.sort_values("ST", ascending=False).variable.head(3))).rename("top3_ST")
    top.to_csv(OUT / "sobol_top3.csv")
    print("\n[1,2] Pf / FORM / runout / sa*"); print(t1.to_string(index=False, float_format=lambda v: f"{v:.5g}"))
    print("\nSobol ST"); print(so.pivot_table(index=["criterion", "variable"], columns="case", values="ST").round(3).to_string())
    print(top.to_string())
    print("\n[4] V2' life at 1200 N"); print(tl.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    # ---------- 5 helper: ranking of Pf across criteria ----------
    rk = t1.pivot(index="case", columns="criterion", values="Pf")
    print("\nPf ordering per case:", {k: " > ".join(r.sort_values(ascending=False).index) for k, r in rk.iterrows()})

    # ---------- 6. Pf vs load ----------
    loads = np.arange(800, 1701, 50)
    n_l = 100_000
    ul = np.random.default_rng(11).random((n_l, len(ft.NAMES)))
    rows = []
    cache = {}
    for cname in ("A_900_cv3", "B_991_cv3"):
        mean, cv = CASES[cname]
        x = physical(ul, mean, cv); x[:, N_FE] = p.RUNOUT_FORCE
        pk1200 = ft.peak_stress(gp, x)
        for L in loads:
            for c in ("goodman", "gerber"):
                f = (ft.g_fun(pk1200 * L / p.RUNOUT_FORCE, x[:, -2], x[:, -1], c) < 0).mean()
                rows.append({"case": cname, "criterion": c, "load_N": L, "Pf": f, "mc_se": np.sqrt(f * (1 - f) / n_l)})
    pl = pd.DataFrame(rows); pl.to_csv(OUT / "pf_vs_load.csv", index=False)
    # load at which Pf = 0.5 (median runout), per case, by interpolation on the curve
    med = []
    for (cn, c), d in pl.groupby(["case", "criterion"]):
        med.append({"case": cn, "criterion": c, "load_at_Pf_0.5_N": float(np.interp(0.5, d.Pf.values, d.load_N.values)),
                    "Pf_at_1200N": float(np.interp(1200, d.load_N.values, d.Pf.values))})
    pd.DataFrame(med).to_csv(OUT / "pf_vs_load_summary.csv", index=False); print("\n[6]"); print(pd.DataFrame(med).to_string(index=False))
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    sty = {("goodman", "A_900_cv3"): ("#4c78a8", "-"), ("goodman", "B_991_cv3"): ("#4c78a8", "--"),
           ("gerber", "A_900_cv3"): ("#f58518", "-"), ("gerber", "B_991_cv3"): ("#f58518", "--")}
    for (c, cn), (col, ls) in sty.items():
        d = pl[(pl.criterion == c) & (pl.case == cn)]
        ax.plot(d.load_N, d.Pf, color=col, ls=ls, marker="o", ms=3, label=f"{c.capitalize()}, UTS {int(CASES[cn][0])} MPa")
    ax.axvline(1200, color="k", lw=1.2); ax.set_xlabel("applied load (N)"); ax.set_ylabel("P_f (no runout at 1e6 cycles)")
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "pf_vs_load_uts.png", dpi=600)


if __name__ == "__main__":
    main()
