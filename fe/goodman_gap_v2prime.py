"""Phase 2, step 6b: load-independent per-sample S-N line V2' and a direct test of the Goodman-gap explanation.

Runout amplitude on the R = 0.1 load line (sm = r*sa, r = (1+R)/(1-R) = 0.55/0.45), sample i with sD_i = ratio_i*UTS_i/Kf, Kf = 1.66:
  Goodman  sa = sD (1 - r sa/UTS)          ->  sa* = sD / (1 + r sD/UTS)  = (ratio UTS/Kf) / (1 + r ratio/Kf)
  Gerber   sa = sD (1 - (r sa/UTS)^2)      ->  k sa^2 + sa - sD = 0, k = sD r^2/UTS^2  ->  sa* = (-1 + sqrt(1 + 4 k sD)) / (2 k)
  SWT      sa = sD sqrt(sa * smax/..)      ->  sqrt(sn*sa) = sa/sqrt(0.45) = sD  ->  sa* = sD sqrt(0.45)
(sn = sa/0.45 is the nominal maximum stress.) The closed forms are verified against a numerical root of ft.g_fun.
V2': N = 1e6 (sa / sa*_i)^(1/b), b = Table 2 slope (fixed), line independent of the applied load.
V1 (Stage E / step 6): N = (sa / (A scale))^(1/b), A scale 1e6^b = A 1e6^b * ratio*UTS/450 (endpoint linear in ratio and UTS).
Samples as in step 6 (GP on all 416 runs, seed 3, 400,000, same life-scatter draw). Life scatter log10 SD 0.2 is ASSUMED.
Outputs go to OUT.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

import fatigue as ft
import params as p
import sensitivity as sens
import stage_e as se
import uncertainty as unc

OUT = Path(__file__).resolve().parents[1] / "results" / "phase2" / "step6_goodman_gap/v2prime"
NMC = 400_000
KF = ft.KT_B
RR = (1 + ft.R) / (1 - ft.R)            # sm / sa on the R = 0.1 load line = 0.55 / 0.45
FORMULAS = {
    "goodman": "sa* = (ratio*UTS/Kf) / (1 + (0.55/0.45)*ratio/Kf)",
    "gerber": "sa* = (-1 + sqrt(1 + 4*k*sD)) / (2*k), sD = ratio*UTS/Kf, k = sD*(0.55/0.45)^2/UTS^2",
    "swt": "sa* = sqrt(0.45) * ratio*UTS/Kf",
}


def sa_star(crit, uts, ratio):
    sd = ratio * uts / KF
    if crit == "goodman":
        return sd / (1 + RR * sd / uts)
    if crit == "gerber":
        k = sd * RR**2 / uts**2
        return (-1 + np.sqrt(1 + 4 * k * sd)) / (2 * k)
    return np.sqrt(0.45) * sd


def sa_star_numeric(crit, uts, ratio):
    """Amplitude at which g = 0 on the R = 0.1 load line, from ft.g_fun itself (peak = sa/0.45*Kt)."""
    f = lambda sa, u, r: float(ft.g_fun(np.array([sa / 0.45 * KF]), u, r, crit)[0])
    return np.array([brentq(f, 1e-3, 5000, args=(u, r), xtol=1e-10) for u, r in zip(uts, ratio)])


def v1_end(crit, uts, ratio):
    A, b = se.TABLE2[crit]
    return A * 1e6**b * ratio * uts / (0.5 * 900.0)


def table_rows(crit, var, scat, n):
    q = np.percentile(n, [5, 50, 95]); f = np.mean(n < 1e6)
    return {"criterion": crit, "variant": var, "scatter_included": scat, "N_p05": q[0], "N_median": q[1], "N_p95": q[2],
            "P_N_below_1e6": f, "mc_se": np.sqrt(f * (1 - f) / len(n))}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    gp, q2, rmse = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
    print(f"GP on {len(df)} runs: Q2 {q2:.4f}, CV RMSE {rmse:.2f} MPa")

    print("\n[1] formulas (R = 0.1 load line, Kf = 1.66)")
    for c, f in FORMULAS.items():
        print(f"  {c:8s} {f}")
    rng = np.random.default_rng(0)
    uts_t, ratio_t = ft.UTS.rvs(300, random_state=1), ft.RATIO.rvs(300, random_state=2)
    chk = {c: float(np.max(np.abs(sa_star(c, uts_t, ratio_t) - sa_star_numeric(c, uts_t, ratio_t)))) for c in FORMULAS}
    print("  max |closed form - numerical root| over 300 random (UTS, ratio):", {k: f"{v:.2e} MPa" for k, v in chk.items()})

    print("\n[2] nominal sa* (UTS 900, ratio 0.5)")
    nom = {c: float(sa_star(c, 900.0, 0.5)) for c in FORMULAS}
    print("  ", {k: round(v, 3) for k, v in nom.items()}, "| Table 2: goodman 198.14, gerber 241.84")
    pd.DataFrame([{"criterion": c, "sa_star_nominal": v, "table2_1e6": {"goodman": 198.14, "gerber": 241.84}.get(c, np.nan),
                   "formula": FORMULAS[c], "closed_form_vs_numeric_max_abs_MPa": chk[c]} for c, v in nom.items()]).to_csv(OUT / "nominal_sa_star.csv", index=False)

    rng = np.random.default_rng(3)
    x = ft.to_physical(rng.random((NMC, len(ft.NAMES))))
    peak = ft.peak_stress(gp, x)
    uts, ratio = x[:, -2], x[:, -1]
    sa = (1 - ft.R) / 2 * peak / KF
    eps = rng.standard_normal(NMC)

    # ---------- 3. V2' ----------
    rows, tab, dis = [], [], []
    for c in ("goodman", "gerber", "swt"):
        fail_ls = ft.g_fun(peak, uts, ratio, c) < 0
        ss = sa_star(c, uts, ratio)
        dis.append({"criterion": c, "method": "V2prime: sa > sa*", "disagreements_with_g<0": int(((sa > ss) != fail_ls).sum())})
        if c == "swt":
            continue
        A, b = se.TABLE2[c]
        n2 = 1e6 * (sa / ss) ** (1 / b)
        ss0 = sa_star(c, np.full(NMC, 900.0), np.full(NMC, 0.5))
        n2_0 = 1e6 * (sa / ss0) ** (1 / b)
        n1 = se.basquin(sa, c, ratio * uts / 450.0)
        n1_0 = se.basquin(sa, c, 1.0)
        dis.append({"criterion": c, "method": "V2prime: N<1e6", "disagreements_with_g<0": int(((n2 < 1e6) != fail_ls).sum())})
        dis.append({"criterion": c, "method": "V1 fixed line: N<1e6", "disagreements_with_g<0": int(((n1 < 1e6) != fail_ls).sum()),
                    "g<0_but_N>=1e6": int((fail_ls & ~(n1 < 1e6)).sum()), "g>=0_but_N<1e6": int((~fail_ls & (n1 < 1e6)).sum())})
        tab += [table_rows(c, "V2prime", "inputs+force only", n2_0), table_rows(c, "V2prime", "+ material scatter", n2),
                table_rows(c, "V2prime", "+ life scatter, log10 SD 0.2 (ASSUMED)", n2 * 10 ** (0.2 * eps)),
                table_rows(c, "V1_fixed_line", "inputs+force only", n1_0), table_rows(c, "V1_fixed_line", "+ material scatter", n1)]
        rows.append({"criterion": c, "Pf_limit_state": fail_ls.mean(), "Pf_V2prime": (n2 < 1e6).mean(), "Pf_V1": (n1 < 1e6).mean(),
                     "mc_se_V2prime": np.sqrt((n2 < 1e6).mean() * (1 - (n2 < 1e6).mean()) / NMC)})
    pd.DataFrame(rows).to_csv(OUT / "pf_summary.csv", index=False)
    d = pd.DataFrame(dis); d.to_csv(OUT / "classification_disagreement.csv", index=False)
    t = pd.DataFrame(tab); t.to_csv(OUT / "table12_v2prime.csv", index=False)
    print("\n[3] Pf"); print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.5g}"))
    print(d.to_string(index=False)); print("\n[5] Table 12 format"); print(t.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    # ---------- 4. endpoint vs ratio ----------
    r_grid = np.array([0.4, 0.5, 0.6, 0.65])
    er, ex = [], []
    for c in ("goodman", "gerber"):
        v1 = v1_end(c, 900.0, r_grid); ss = sa_star(c, 900.0, r_grid)
        for r_, a_, b_ in zip(r_grid, v1, ss):
            er.append({"criterion": c, "ratio": r_, "V1_endpoint_MPa": a_, "sa_star_MPa": b_, "V1_over_sa_star": a_ / b_})
        v1s, sss = v1_end(c, uts, ratio), sa_star(c, uts, ratio)
        er.append({"criterion": c, "ratio": "mean over 400k samples", "V1_endpoint_MPa": v1s.mean(), "sa_star_MPa": sss.mean(),
                   "V1_over_sa_star": (v1s / sss).mean()})
        er.append({"criterion": c, "ratio": "sd over samples", "V1_over_sa_star": (v1s / sss).std()})
        er.append({"criterion": c, "ratio": "UTS-only variation at ratio 0.5 (min..max of V1/sa* over UTS 5-95%)",
                   "V1_over_sa_star": float(np.ptp(v1_end(c, np.percentile(uts, [5, 95]), 0.5) / sa_star(c, np.percentile(uts, [5, 95]), 0.5)))})
        # where do the disagreements sit?
        fail_ls = ft.g_fun(peak, uts, ratio, c) < 0
        n1 = se.basquin(sa, c, ratio * uts / 450.0)
        dg = fail_ls != (n1 < 1e6)
        for lab, m in (("ratio > 0.55", ratio > 0.55), ("ratio > 0.50", ratio > 0.50), ("ratio <= 0.50", ratio <= 0.50)):
            ex.append({"criterion": c, "region": lab, "share_of_samples": m.mean(), "disagreements": int((dg & m).sum()),
                       "share_of_disagreements": (dg & m).sum() / dg.sum(), "mean_V1_over_sa_star": (v1_end(c, uts, ratio)[m] / sa_star(c, uts, ratio)[m]).mean()})
    e = pd.DataFrame(er); e.to_csv(OUT / "endpoint_ratio.csv", index=False)
    x_ = pd.DataFrame(ex); x_.to_csv(OUT / "disagreement_by_ratio.csv", index=False)
    print("\n[4] V1 endpoint / sa*"); print(e.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(x_.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    # linearity: V1/sa* is linear in ratio -> fit
    for c in ("goodman", "gerber"):
        rg = np.linspace(0.3, 0.7, 9)
        q = np.polyfit(rg, v1_end(c, 900.0, rg) / sa_star(c, 900.0, rg), 1)
        print(f"  {c}: V1/sa* ~ {q[1]:.4f} + {q[0]:.4f} * ratio  (linear fit over ratio 0.3-0.7; residual max "
              f"{np.max(np.abs(np.polyval(q, rg) - v1_end(c, 900.0, rg) / sa_star(c, 900.0, rg))):.2e})")

    # ---------- figure ----------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    rg = np.linspace(0.3, 0.7, 200)
    for a, c in zip(ax, ("goodman", "gerber")):
        a.plot(rg, v1_end(c, 900.0, rg), color="#4c78a8", label="V1 endpoint (Table 2 A, b; linear in ratio)")
        a.plot(rg, sa_star(c, 900.0, rg), color="#f58518", label="sa*: load-line runout amplitude (V2')")
        a.axvspan(0.55, 0.7, color="grey", alpha=.15, label="ratio > 0.55")
        a.axvline(0.5, color="k", lw=.8, ls="--")
        a.set_xlabel("endurance ratio (UTS = 900 MPa)"); a.set_ylabel("amplitude at 1e6 cycles (MPa)")
        a.set_title(c.capitalize()); a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "endpoints_vs_ratio.png", dpi=600)


if __name__ == "__main__":
    main()
