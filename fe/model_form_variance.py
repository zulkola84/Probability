"""Phase 3: how much of the failure-probability uncertainty is model form, and what the survival at 1200 N says.

Epistemic factors (each a discrete set of equally weighted levels, ASSUMED prior):
  criterion        Goodman, Gerber, SWT
  notch factor Kf  Kt = 1.66 (Bologna et al.), Peterson 1.641, Okazaki q = 0.25 -> 1.165
  mean UTS         900, 945.5, 991 MPa (certificate values of the two plates in the cited studies)
  scatter set      low (ratio SD 0.03, UTS CV 2 %), base (0.05, 3 %), high (0.08, 5 %)
Aleatory inputs: thickness, hole diameter, E, nu, lever arm, force, UTS, endurance ratio (Table 2 of the paper).
For every configuration e, Pf(e) is estimated by Monte Carlo with common random numbers on the 416-run GP.

Law of total variance of the failure indicator I = 1{g < 0}:
  Var(I) = Var_e[Pf(e)] (epistemic) + E_e[Pf(e)(1 - Pf(e))] (aleatory, within a configuration).
Var_e[Pf(e)] is split by exact ANOVA over the factor grid (first-order and total indices).
Survival at 1200 N: likelihood of a plate surviving is 1 - Pf(e); Bayes weights with equal priors.
"""
import itertools
import json

import numpy as np
import pandas as pd

import fatigue as ft
import params as p
import sensitivity as sens
import stage_d as sd
import uncertainty as unc

OUT = p.RESULTS_DIR / "phase3"
N_MC = 400_000
CRIT = ("goodman", "gerber", "swt")
KF = {"Kt=1.66": ft.KT_B, "Peterson": 1 + 0.9 * (ft.KT_B - 1) * 0 + (sd.kf_options()[0]["peterson"] - 1), "q=0.25": 1 + 0.25 * (ft.KT_B - 1)}
KF["Peterson"] = sd.kf_options()[0]["peterson"]
UTS_MEAN = (900.0, 945.5, 991.0)
SCATTER = {"low": (0.03, 0.02), "base": (0.05, 0.03), "high": (0.08, 0.05)}   # (ratio SD, UTS CV)


def configs_pf(gp):
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    rng = np.random.default_rng(7)
    u = rng.random((N_MC, len(sd.NAMES)))
    peak = sd.peak_of(gp, u)
    rows = []
    for (sname, (rsd, cv)) in SCATTER.items():
        _, uts0, ratio = sd.physical(u, cv, rsd, 0.0)      # UTS with mean 900
        for m in UTS_MEAN:
            uts = uts0 * m / 900.0
            for crit, (kn, kv) in itertools.product(CRIT, KF.items()):
                g = ft.g_fun(peak, uts, ratio, crit, kv)
                pf = float(np.mean(g < 0))
                rows.append({"scatter": sname, "uts_mean": m, "criterion": crit, "kf": kn, "Pf": pf})
    return pd.DataFrame(rows)


def anova(d, factors):
    """Exact ANOVA of Pf over a full grid: first-order and total indices (share of Var_e[Pf])."""
    y = d["Pf"].to_numpy()
    var = y.var()
    out = {}
    for f in factors:
        m = d.groupby(f)["Pf"].mean()
        s1 = ((m - y.mean()) ** 2).mean() / var
        rest = [g for g in factors if g != f]
        mm = d.groupby(rest)["Pf"].mean()                    # variance explained by all other factors
        st = 1 - ((mm - y.mean()) ** 2).mean() / var
        out[f] = {"S1": float(s1), "ST": float(st)}
    return float(var), out


def summarise(d, label, factors):
    pf = d["Pf"].to_numpy()
    epi = pf.var()
    ale = np.mean(pf * (1 - pf))
    tot = epi + ale
    var, idx = anova(d, factors)
    return {"label": label, "n_configs": len(d), "mean_Pf": float(pf.mean()), "min_Pf": float(pf.min()),
            "max_Pf": float(pf.max()), "epistemic_share": float(epi / tot), "aleatory_share": float(ale / tot),
            "anova": idx}


def main():
    OUT.mkdir(exist_ok=True)
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    gp, q2, rmse = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
    d = configs_pf(gp)
    d.to_csv(OUT / "Pf_by_configuration.csv", index=False)
    fac = ["criterion", "kf", "uts_mean", "scatter"]
    res = [summarise(d, "all 81 configurations", fac),
           summarise(d[d.kf != "q=0.25"], "without q = 0.25 (54)", fac),
           summarise(d[d.criterion != "swt"], "Goodman and Gerber only (54)", fac),
           summarise(d[(d.criterion != "swt") & (d.kf != "q=0.25")], "Goodman, Gerber; Kf = Kt or Peterson (36)", fac),
           summarise(d[(d.uts_mean == 900.0)], "mean UTS fixed at 900 (27)", ["criterion", "kf", "scatter"])]
    # survival likelihood at 1200 N, equal priors
    d["like"] = 1 - d["Pf"]
    post = {}
    for f in ("criterion", "kf", "uts_mean", "scatter"):
        w = d.groupby(f)["like"].mean()
        post[f] = (w / w.sum()).round(4).to_dict()
    base = d[(d.scatter == "base") & (d.uts_mean == 900.0) & (d.kf == "Kt=1.66")].set_index("criterion")["like"]
    lr = {"Gerber_over_Goodman_baseline": float(base["gerber"] / base["goodman"]),
          "Gerber_over_SWT_baseline": float(base["gerber"] / base["swt"]),
          "like_goodman": float(base["goodman"]), "like_gerber": float(base["gerber"]), "like_swt": float(base["swt"])}
    marg = d.groupby("criterion")["like"].mean()
    lr["marginal_like_by_criterion"] = marg.round(4).to_dict()
    lr["marginal_Gerber_over_Goodman"] = float(marg["gerber"] / marg["goodman"])
    json.dump({"summaries": res, "posterior_weights_equal_priors": post, "likelihood": lr,
               "surrogate_q2": q2, "surrogate_rmse": rmse}, open(OUT / "model_form_variance.json", "w"), indent=2)
    print(json.dumps({"summaries": res, "posterior": post, "likelihood": lr}, indent=1))


if __name__ == "__main__":
    main()
