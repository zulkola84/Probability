"""Fatigue failure probability of the Bologna et al. 2023 plate at the 1e6-cycle runout condition.

Limit state (the paper's rules, Bologna et al. 2023):
    nominal max stress  sn_max = peak_vm / KT_B          (peak from the FE surrogate, scaled by force/1200)
    R = 0.1:            sa = (1-R)/2 * sn_max ,  sm = (1+R)/2 * sn_max
    fatigue limit       sD = ENDURANCE_RATIO * UTS / KF   (paper: 0.5 UTS, Kf = Kt,b = 1.66)
    Goodman             allowable = sD * (1 - sm/UTS)          (lower bound)
    Gerber              allowable = sD * (1 - (sm/UTS)**2)     (upper bound)
    g = allowable - sa ;  failure (no runout at 1e6 cycles) when g < 0.

Random inputs: the five FE inputs and the force from uncertainty.py, plus UTS and the endurance ratio
(extra spreads below are ASSUMED; the paper gives none).
Needs results/fe_samples.csv. Usage: python fatigue.py [N_MC, default 400000]
Writes results/fatigue_summary.json, fatigue_sobol.csv, fatigue.png
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq
from SALib.analyze import sobol
from SALib.sample import sobol as sobol_sample

import params as p
import sensitivity as sens
import uncertainty as unc

R = 0.1
KT_B = 1.66                       # paper: Kf = Kt,b = 1.66
UTS = stats.lognorm(s=np.sqrt(np.log(1 + 0.03**2)), scale=900 * np.exp(-np.log(1 + 0.03**2) / 2))  # ASSUMED 3 % CV, mean from certificate
RATIO = stats.norm(0.5, 0.05)     # ASSUMED scatter around the paper's 0.5 rule
NAMES = unc.FE_NAMES + ["force", "uts", "endurance_ratio"]
DISTS = [unc.FE_INPUTS[k] for k in unc.FE_NAMES] + [unc.FORCE, UTS, RATIO]


def uts_dist(cv, mean=900.0):
    """Lognormal UTS with the given CV (mean 900 MPa from the certificate, Bologna et al. 2023)."""
    s2 = np.log(1 + cv**2)
    return stats.lognorm(s=np.sqrt(s2), scale=mean * np.exp(-s2 / 2))


def ratio_dist(sd, mean=0.5):
    return stats.norm(mean, sd)


# Stage D cases. "assumed" is the original case above (the default). Literature values (Stage D REPORT.md):
#   CV 1.1 %: Okazaki 2012, grade 4 Ti (ISO 5832-2) annealed UTS 700 +- 8 MPa, n = 3 (doi 10.3390/ma5122981)
#   CV 8.5 %: Suzuki 2018, grade 4 machined implant specimens, static yield FORCE 1088 +- 93 N (doi 10.4012/dmj.2017-125),
#             a component-level proxy for strength scatter (undefined +-, n unknown)
# No published fatigue-strength SD for grade 4 was found: the endurance-ratio SD stays ASSUMED in every case.
CASES = {"assumed": dict(uts_cv=0.03, ratio_sd=0.05), "lit_material": dict(uts_cv=0.011, ratio_sd=0.05),
         "lit_component": dict(uts_cv=0.085, ratio_sd=0.05)}


def to_physical(u):
    u = np.clip(u, 1e-9, 1 - 1e-9)
    return np.column_stack([d.ppf(u[:, i]) for i, d in enumerate(DISTS)])


def peak_stress(gp, x):
    """Peak von Mises (MPa) for physical samples x (n, 8)."""
    n_fe = len(unc.FE_NAMES)
    xs = (x[:, :n_fe] - sens.MU) / sens.SD
    out = np.empty(len(x))
    for i in range(0, len(x), 100_000):
        out[i:i + 100_000] = gp.predict(xs[i:i + 100_000])
    return out * x[:, n_fe] / p.RUNOUT_FORCE


def g_fun(peak, uts, ratio, criterion, kf=KT_B):
    """kf is the fatigue notch factor applied to the endurance limit (default Kf = Kt,b); nominal stress
    is always peak / KT_B. criterion: goodman | gerber | swt (Smith-Watson-Topper, R = -1 limit sd)."""
    sn = peak / KT_B
    sa, sm = (1 - R) / 2 * sn, (1 + R) / 2 * sn
    sd = ratio * uts / kf
    if criterion == "swt":
        return sd - np.sqrt(sn * sa)          # sqrt(smax * sa) compared with the R = -1 limit
    k = sm / uts
    allow = sd * (1 - k) if criterion == "goodman" else sd * (1 - k**2)
    return allow - sa


def g_x(gp, x, criterion):
    return g_fun(peak_stress(gp, x), x[:, -2], x[:, -1], criterion)


def mc(gp, n, criterion, rng, batches=20):
    pf = []
    for _ in range(batches):
        x = to_physical(rng.random((n // batches, len(NAMES))))
        pf.append(np.mean(g_x(gp, x, criterion) < 0))
    pf = np.array(pf)
    boot = [rng.choice(pf, batches).mean() for _ in range(2000)]
    return pf.mean(), np.percentile(boot, [2.5, 97.5])


def form(gp, criterion, iters=60):
    """First-order reliability (HL-RF) in standard-normal space: a check independent of random sampling."""
    def g_u(z):
        u = stats.norm.cdf(z)[None, :]
        return float(g_x(gp, to_physical(u), criterion)[0])
    z = np.zeros(len(NAMES))
    for _ in range(iters):
        g0 = g_u(z)
        grad = np.array([(g_u(z + 1e-4 * e) - g0) / 1e-4 for e in np.eye(len(z))])
        z_new = (grad @ z - g0) / (grad @ grad) * grad
        done = np.linalg.norm(z_new - z) < 1e-6
        z = z_new
        if done:
            break
    beta = np.linalg.norm(z)
    if g_u(np.zeros(len(NAMES))) < 0:      # nominal plate already fails: the failure region contains the origin
        return stats.norm.cdf(beta), -beta
    return stats.norm.cdf(-beta), beta


def runout_load(gp, criterion):
    """Load at which the nominal plate (every input at its mean) just reaches g = 0."""
    mean = np.array([[d.mean() for d in DISTS]])
    n_fe = len(unc.FE_NAMES)

    def f(load):
        x = mean.copy(); x[0, n_fe] = load
        return g_x(gp, x, criterion)[0]
    return brentq(f, 100, 10000)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 400_000
    rng = np.random.default_rng(1)
    df = pd.read_csv(p.RESULTS_DIR / "fe_samples.csv")
    gp, q2, rmse = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
    print(f"Surrogate for peak stress: Q2={q2:.4f}, CV RMSE={rmse:.2f} MPa")

    res = {"surrogate_q2": q2, "surrogate_cv_rmse_MPa": rmse, "n_mc": n, "criteria": {}}
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    loads = np.arange(800, 1701, 100)
    for crit, col in (("goodman", "#4c78a8"), ("gerber", "#f58518")):
        pf, ci = mc(gp, n, crit, rng)
        pf_f, beta = form(gp, crit)
        L = runout_load(gp, crit)
        # sampled P_f versus load level (other inputs random, force fixed)
        curve = []
        for L_i in loads:
            x = to_physical(rng.random((50_000, len(NAMES)))); x[:, len(unc.FE_NAMES)] = L_i
            curve.append(np.mean(g_x(gp, x, crit) < 0))
        # Sobol on g
        prob = {"num_vars": len(NAMES), "names": NAMES, "bounds": [[0, 1]] * len(NAMES)}
        u = sobol_sample.sample(prob, 2**13, calc_second_order=False, seed=1)
        s = sobol.analyze(prob, g_x(gp, to_physical(u), crit), calc_second_order=False, seed=1)
        sob = pd.DataFrame({"variable": NAMES, "S1": s["S1"], "S1_conf": s["S1_conf"], "ST": s["ST"], "ST_conf": s["ST_conf"]})
        sob.insert(0, "criterion", crit)
        res["criteria"][crit] = {"Pf_at_1200N": pf, "Pf_ci95": list(ci), "Pf_FORM": pf_f, "beta": beta,
                                 "predicted_runout_N": L, "runout_error_pct": (L / p.RUNOUT_FORCE - 1) * 100}
        print(f"\n{crit.upper()}: P_f(1200 N) = {pf:.4f} (95% CI {ci[0]:.4f}-{ci[1]:.4f}), FORM {pf_f:.4f} (beta {beta:.2f});"
              f" predicted runout {L:.0f} N ({(L / p.RUNOUT_FORCE - 1) * 100:+.1f}% vs 1200 N)")
        print(sob.round(3).to_string(index=False))
        sob.to_csv(p.RESULTS_DIR / f"fatigue_sobol_{crit}.csv", index=False)
        x = to_physical(rng.random((100_000, len(NAMES))))
        ax[0].hist(g_x(gp, x, crit), 80, alpha=.6, color=col, label=crit)
        ax[1].plot(loads, curve, "-o", color=col, label=crit)
        d = sob.set_index("variable")["ST"]
        ax[2].bar(np.arange(len(d)) + (0.2 if crit == "gerber" else -0.2), d.values, 0.4, color=col, label=crit)
    ax[0].axvline(0, color="k"); ax[0].set_xlabel("g (MPa)"); ax[0].legend()
    ax[1].axvline(p.RUNOUT_FORCE, color="k", ls="--"); ax[1].set_xlabel("load (N)"); ax[1].set_ylabel("P(no runout)"); ax[1].legend()
    ax[2].set_xticks(range(len(NAMES))); ax[2].set_xticklabels(NAMES, rotation=40, ha="right"); ax[2].set_ylabel("Sobol ST of g"); ax[2].legend()
    fig.tight_layout(); fig.savefig(p.RESULTS_DIR / "fatigue.png", dpi=150)
    res["note"] = ("UTS CV 3%, endurance-ratio SD 0.05, force CV 2% and the FE-input spreads are assumptions; "
                   "the mean values follow Bologna et al. 2023. Predicted runout is a validation check against 1200 N.")
    json.dump(res, open(p.RESULTS_DIR / "fatigue_summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
