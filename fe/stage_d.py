"""Stage D: material scatter and model-form uncertainty for the fatigue limit state.

Surrogate: GP on results/stageA/fe_samples_al.csv (416 FE runs; Stage A). The surrogate peak stress does not depend on
UTS or the endurance ratio, so one set of FE-input/force samples is evaluated once and reused for every case
(common random numbers); only g changes. Writes results/stageD/*.csv.

Parts: (A) literature vs assumed scatter; (B) criteria Goodman/Gerber/SWT x Kf = Kt, Peterson, Okazaki q;
(C) robustness grid: endurance-ratio SD 0.03/0.05/0.08 x UTS CV 2/3/5 % x UTS-ratio correlation 0 / 0.7,
P_f and Sobol ST ranking. Correlation is a Gaussian copula between the UTS and ratio drivers (ASSUMED rho = 0.7, brief).
With correlation, Sobol indices refer to the independent drivers z_uts and z_ratio_residual, not the marginals.
"""
import itertools

import numpy as np
import pandas as pd
from SALib.analyze import sobol
from SALib.sample import sobol as sobol_sample
from scipy import stats
from scipy.special import ndtri

import fatigue as ft
import params as p
import sensitivity as sens
import uncertainty as unc

OUT = p.RESULTS_DIR / "stageD"
N_MC = 400_000
N_SOBOL = 2**13
N_FE = len(unc.FE_NAMES)
NAMES = ft.NAMES          # 5 FE inputs, force, uts, endurance_ratio
CRITERIA = ("goodman", "gerber", "swt")
Z = ndtri


def kf_options():
    """Fatigue notch factors applied to the endurance limit (nominal stress unchanged)."""
    kt = ft.KT_B
    a = 0.0254 * (2070 / 900.0) ** 1.8                   # mm; Peterson constant, a STEEL correlation applied to Ti (ASSUMED)
    r = p.HOLE_DIAMETER / 2
    q_pet = 1 / (1 + a / r)
    return {"kf=Kt": kt, "peterson": 1 + q_pet * (kt - 1),
            "okazaki_q0.25": 1 + 0.25 * (kt - 1)}, {"peterson_a_mm": a, "peterson_q": q_pet}


def physical(u, cv, sd, rho):
    """Map unit-cube columns to physical variables for a given UTS CV, ratio SD and correlation."""
    x = np.column_stack([unc.FE_INPUTS[k].ppf(np.clip(u[:, i], 1e-9, 1 - 1e-9)) for i, k in enumerate(unc.FE_NAMES)]
                        + [unc.FORCE.ppf(np.clip(u[:, N_FE], 1e-9, 1 - 1e-9))])
    u6, u7 = np.clip(u[:, N_FE + 1], 1e-9, 1 - 1e-9), np.clip(u[:, N_FE + 2], 1e-9, 1 - 1e-9)
    uts = ft.uts_dist(cv).ppf(u6)
    ratio = ft.ratio_dist(sd).mean() + sd * (rho * Z(u6) + np.sqrt(1 - rho**2) * Z(u7))
    return x, uts, ratio


def peak_of(gp, u):
    x, _, _ = physical(u, 0.03, 0.05, 0.0)
    return ft.peak_stress(gp, x)


def pf(peak, uts, ratio, crit, kf=ft.KT_B):
    g = ft.g_fun(peak, uts, ratio, crit, kf)
    p_ = np.mean(g < 0)
    return p_, np.sqrt(p_ * (1 - p_) / len(g))


def main():
    OUT.mkdir(exist_ok=True)
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    gp, q2, rmse = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
    print(f"surrogate: {len(df)} runs, Q2 {q2:.4f}, CV RMSE {rmse:.2f} MPa")
    rng = np.random.default_rng(1)
    u = rng.random((N_MC, len(NAMES)))
    peak = peak_of(gp, u)
    kfs, kinfo = kf_options()
    print("Kf options", {k: round(v, 4) for k, v in kfs.items()}, kinfo)

    # (A) literature vs assumed scatter, Kf = Kt
    rows = []
    for name, c in ft.CASES.items():
        _, uts, ratio = physical(u, c["uts_cv"], c["ratio_sd"], 0.0)
        for crit in CRITERIA:
            v, se = pf(peak, uts, ratio, crit)
            rows.append({"case": name, **c, "criterion": crit, "Pf": v, "mc_se": se})
    a = pd.DataFrame(rows); a.to_csv(OUT / "A_literature_cases.csv", index=False); print(a.round(4).to_string(index=False))

    # (B) criterion x Kf, original (assumed) scatter
    rows = []
    _, uts, ratio = physical(u, 0.03, 0.05, 0.0)
    for crit, (kn, kv) in itertools.product(CRITERIA, kfs.items()):
        v, se = pf(peak, uts, ratio, crit, kv)
        mean_x = np.array([[d.mean() for d in unc.FE_INPUTS.values()] + [p.RUNOUT_FORCE]])
        pk0 = ft.peak_stress(gp, mean_x)
        rows.append({"criterion": crit, "kf_option": kn, "Kf": kv, "Pf": v, "mc_se": se,
                     "g_nominal_MPa": float(ft.g_fun(pk0, 900.0, 0.5, crit, kv)[0])})
    b = pd.DataFrame(rows); b.to_csv(OUT / "B_criteria_kf.csv", index=False); print(b.round(4).to_string(index=False))

    # (C) robustness grid
    prob = {"num_vars": len(NAMES), "names": NAMES, "bounds": [[0, 1]] * len(NAMES)}
    us = sobol_sample.sample(prob, N_SOBOL, calc_second_order=False, seed=1)
    peak_s = peak_of(gp, us)
    rows, srows = [], []
    for sd, cv, rho in itertools.product((0.03, 0.05, 0.08), (0.02, 0.03, 0.05), (0.0, 0.7)):
        _, uts, ratio = physical(u, cv, sd, rho)
        _, uts_s, ratio_s = physical(us, cv, sd, rho)
        for crit in CRITERIA:
            v, se = pf(peak, uts, ratio, crit)
            s = sobol.analyze(prob, ft.g_fun(peak_s, uts_s, ratio_s, crit), calc_second_order=False, seed=1)
            st = dict(zip(NAMES, s["ST"]))
            rank = sorted(st, key=st.get, reverse=True)
            rows.append({"ratio_sd": sd, "uts_cv": cv, "rho": rho, "criterion": crit, "Pf": v, "mc_se": se,
                         "rank1": rank[0], "rank2": rank[1], "rank3": rank[2], "ST_ratio": st["endurance_ratio"],
                         "ST_uts": st["uts"], "ST_force": st["force"], "ST_thickness": st["thickness"]})
            for n in NAMES:
                srows.append({"ratio_sd": sd, "uts_cv": cv, "rho": rho, "criterion": crit, "variable": n,
                              "S1": s["S1"][NAMES.index(n)], "ST": st[n]})
    c = pd.DataFrame(rows); c.to_csv(OUT / "C_robustness.csv", index=False)
    pd.DataFrame(srows).to_csv(OUT / "C_sobol_all.csv", index=False)
    print(c.round(3).to_string(index=False))
    pd.Series({**kinfo, "surrogate_q2": q2, "surrogate_cv_rmse": rmse, "n_fe": len(df), "n_mc": N_MC,
               "n_sobol_base": N_SOBOL}).to_json(OUT / "meta.json", indent=2)


if __name__ == "__main__":
    main()
