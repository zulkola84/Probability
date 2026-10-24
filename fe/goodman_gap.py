"""Phase 2, step 6: why does the Goodman P_f from the Basquin life (Stage E, 0.937) differ from the limit state (0.980)?

Imports stage_e (constants TABLE2, basquin), fatigue and sensitivity unchanged. Same MC samples as Stage E part (1):
GP on all 416 runs, rng seed 3, 400,000 samples, eps drawn after the inputs. Variants of the life model:
  V1  fixed line (Stage E): N = (sa / (A * scale))^(1/b), A, b = Table 2 constants, scale = ratio*UTS/450
  V2  per-sample line: log-log line through (1e3, anchor_i) and (1e6, allowable_i), allowable_i from the criterion with
      that sample's UTS, ratio and mean stress; anchor_i = 0.9 * 0.45 * UTS_i (amplitude at R = 0.1; this reproduces the
      Table 2 value 364.5 MPa at UTS 900, so it is the Table 2 reading, ASSUMED; 0.9*UTS_i as an amplitude would not)
  V3  DIAGNOSTIC: fixed line like V1 but A recalibrated once so that A*1e6^b equals the limit-state allowable at the
      nominal inputs (isolates the endpoint offset from the sample-to-sample mean-stress dependence)
Life scatter log10 SD 0.2 is ASSUMED (as Stage E). All outputs go to OUT (outside results/).
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fatigue as ft
import params as p
import sensitivity as sens
import stage_e as se
import uncertainty as unc

OUT = Path(__file__).resolve().parents[1] / "results" / "phase2" / "step6_goodman_gap"
NMC = 400_000
CRITS = ("goodman", "gerber")
LIFE_SD = se.LIFE_SD_LOG10


def allowable(peak, uts, ratio, crit, kf=ft.KT_B):
    """Allowable amplitude at 1e6 cycles of the limit state (g = allowable - sa), as in ft.g_fun."""
    return ft.g_fun(peak, uts, ratio, crit, kf) + (1 - ft.R) / 2 * peak / ft.KT_B


def life_v2(sa, allow, uts):
    anchor = 0.9 * 0.45 * uts                                  # ASSUMED reading, see docstring
    with np.errstate(all="ignore"):
        b = (np.log(allow) - np.log(anchor)) / (np.log(1e6) - np.log(1e3))
        n = 1e6 * (sa / allow) ** (1 / b)
    n = np.where(allow <= 0, 0.0, n)
    return n, b


def stats_row(crit, var, scat, n):
    q = np.percentile(n, [5, 50, 95])
    return {"criterion": crit, "variant": var, "scatter_included": scat, "N_p05": q[0], "N_median": q[1], "N_p95": q[2],
            "P_N_below_1e6": np.mean(n < 1e6), "mc_se": np.sqrt(np.mean(n < 1e6) * (1 - np.mean(n < 1e6)) / len(n))}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    gp, q2, rmse = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
    print(f"GP on {len(df)} runs: Q2 {q2:.4f}, CV RMSE {rmse:.2f} MPa")

    # ---------- 1. constants and formulas ----------
    print("\n[1] Basquin constants (Bologna 2023 Table 2): FIXED constants, not recomputed per sample")
    for c, (A, b) in se.TABLE2.items():
        print(f"  {c:8s} A = {A} MPa, b = {b}; A*1e3^b = {A * 1e3 ** b:.2f}, A*1e6^b = {A * 1e6 ** b:.2f} MPa")
    print("  Per-sample material scatter enters only as A_i = A * ratio_i*UTS_i/450 (stage_e.py: scale), b stays fixed.")
    print("  Nominal amplitude, life (stage_e.py):   sa = (1-R)/2 * peak / KT_B = 0.45 * peak / 1.66  (R = 0.1, KT_B = 1.66)")
    print("  Nominal amplitude, limit state (fatigue.py g_fun): sn = peak/KT_B; sa = (1-R)/2*sn = 0.45*peak/1.66 -> identical")
    print("  Resistance, life:        sa_fail(1e6) = A*scale*1e6^b = (A*1e6^b/450) * ratio*UTS   [no Kf, no per-sample mean stress]")
    print("  Resistance, limit state: allow = ratio*UTS/Kf * (1 - sm/UTS) (Goodman) or *(1-(sm/UTS)^2) (Gerber),"
          " sm = (1+R)/2*sn = 0.55*peak/1.66, Kf = Kt = 1.66")
    print("  Differences: (i) Kf: 1.66 in the limit state, not explicit in the life line; (ii) the mean-stress correction is"
          " recomputed per sample (via sm and UTS_i) in the limit state, and fixed at the Table 2 endpoint in V1.")

    # ---------- samples ----------
    rng = np.random.default_rng(3)
    x = ft.to_physical(rng.random((NMC, len(ft.NAMES))))
    peak = ft.peak_stress(gp, x)
    uts, ratio = x[:, -2], x[:, -1]
    sa = (1 - ft.R) / 2 * peak / ft.KT_B
    scale = ratio * uts / (0.5 * 900.0)
    eps = rng.standard_normal(NMC)

    # ---------- 4. line-through-mean check ----------
    print("\n[4] fixed-line check at nominal inputs")
    nom_x = np.array([[d.mean() for d in ft.DISTS]])
    pk_nom = ft.peak_stress(gp, nom_x)
    rows4 = []
    for c, (A, b) in se.TABLE2.items():
        end = A * 1e6 ** b
        paper_end = {"goodman": 198.14, "gerber": 241.84}[c]
        ls_end = float(allowable(pk_nom, nom_x[:, -2], nom_x[:, -1], c)[0])
        ls_end_paper = float(allowable(np.array([p.PAPER_PEAK_FE]), np.array([900.0]), np.array([0.5]), c)[0])
        ls_end_an = float(allowable(np.array([p.PAPER_PEAK_ANALYTICAL]), np.array([900.0]), np.array([0.5]), c)[0])
        n_pe = {k: se.basquin((1 - ft.R) / 2 * v / ft.KT_B, c) for k, v in
                (("paper_analytical_898.1", p.PAPER_PEAK_ANALYTICAL), ("paper_FE_857.8", p.PAPER_PEAK_FE))}
        n_gp = se.basquin((1 - ft.R) / 2 * float(pk_nom[0]) / ft.KT_B, c)
        rows4.append({"criterion": c, "A_b_endpoint_1e6_MPa": end, "table2_1e6_MPa": paper_end,
                      "endpoint_vs_table2_pct": 100 * (end / paper_end - 1),
                      "limit_state_allowable_nominal_GP_MPa": ls_end, "limit_state_allowable_paper_FE_peak_MPa": ls_end_paper,
                      "limit_state_allowable_paper_analytical_peak_MPa": ls_end_an,
                      "limit_state_vs_fixed_endpoint_pct": 100 * (ls_end / end - 1),
                      "nominal_GP_peak_MPa": float(pk_nom[0]), "N_at_nominal_GP_peak": n_gp,
                      "N_paper_analytical_peak": n_pe["paper_analytical_898.1"], "N_paper_FE_peak": n_pe["paper_FE_857.8"]})
    c4 = pd.DataFrame(rows4)
    mean_an = np.mean([r["N_paper_analytical_peak"] for r in rows4]); mean_fe = np.mean([r["N_paper_FE_peak"] for r in rows4])
    c4["mean_of_criteria_N_analytical"] = mean_an; c4["mean_of_criteria_N_FE"] = mean_fe
    c4["published_N_analytical"] = 575_000; c4["published_N_FE"] = 1_000_000
    c4.to_csv(OUT / "line_through_mean_check.csv", index=False)
    print(c4.T.to_string())

    # ---------- 2/3. classification and variants ----------
    nominal_allow = {c: float(allowable(pk_nom, nom_x[:, -2], nom_x[:, -1], c)[0]) for c in CRITS}
    pf_rows, dis_rows, tab_rows, lives = [], [], [], {}
    for c in CRITS:
        A, b = se.TABLE2[c]
        allow = allowable(peak, uts, ratio, c)
        g = allow - sa
        fail_ls = g < 0
        n1 = se.basquin(sa, c, scale)                                   # V1
        a3 = nominal_allow[c] / 1e6 ** b                                # V3 recalibrated A
        n3 = (sa / (a3 * scale)) ** (1 / b)
        n2, b2 = life_v2(sa, allow, uts)                                # V2
        n1_mean = se.basquin(sa, c, 1.0)
        lives[c] = {"V1": n1, "V2": n2, "V3": n3}
        for lab, n in (("limit_state_P(g<0)", None), ("V1_fixed_line", n1), ("V2_per_sample_line", n2), ("V3_recalibrated_fixed_line_DIAGNOSTIC", n3)):
            f = fail_ls if n is None else n < 1e6
            pf_rows.append({"criterion": c, "method": lab, "Pf": f.mean(), "mc_se": np.sqrt(f.mean() * (1 - f.mean()) / NMC)})
            if n is not None:
                dis_rows.append({"criterion": c, "method": lab, "disagree_total": int((f != fail_ls).sum()),
                                 "g_lt0_but_N_ge_1e6": int((fail_ls & ~f).sum()), "g_ge0_but_N_lt_1e6": int((~fail_ls & f).sum())})
        # numbers supporting the cause
        end1 = A * scale * 1e6 ** b
        pf_rows.append({"criterion": c, "method": "info: mean fixed-line endpoint (MPa)", "Pf": end1.mean(), "mc_se": np.nan})
        pf_rows.append({"criterion": c, "method": "info: mean limit-state allowable (MPa)", "Pf": allow.mean(), "mc_se": np.nan})
        pf_rows.append({"criterion": c, "method": "info: mean(endpoint_fixed/allowable)", "Pf": (end1 / allow).mean(), "mc_se": np.nan})
        pf_rows.append({"criterion": c, "method": "info: sd(endpoint_fixed/allowable)", "Pf": (end1 / allow).std(), "mc_se": np.nan})
        pf_rows.append({"criterion": c, "method": "info: corr(log endpoint_fixed, log allowable)",
                        "Pf": np.corrcoef(np.log(end1), np.log(np.maximum(allow, 1e-9)))[0, 1], "mc_se": np.nan})
        pf_rows.append({"criterion": c, "method": "info: V2 non-positive allowable (count)", "Pf": int((allow <= 0).sum()), "mc_se": np.nan})
        pf_rows.append({"criterion": c, "method": "info: V2 slope b_i mean", "Pf": np.nanmean(b2), "mc_se": np.nan})
        pf_rows.append({"criterion": c, "method": "info: V2 slope b_i sd", "Pf": np.nanstd(b2), "mc_se": np.nan})
        pf_rows.append({"criterion": c, "method": "info: V3 recalibrated A (MPa)", "Pf": a3, "mc_se": np.nan})
        pf_rows.append({"criterion": c, "method": "info: nominal limit-state allowable (MPa)", "Pf": nominal_allow[c], "mc_se": np.nan})

        # ---------- 6. Table 12 format ----------
        one = np.ones(NMC)
        sa0 = sa
        # inputs+force only: ratio = 0.5, UTS = 900 fixed
        allow0 = allowable(peak, np.full(NMC, 900.0), np.full(NMC, 0.5), c)
        n2_0, _ = life_v2(sa0, allow0, np.full(NMC, 900.0))
        for var, nfun in (("V1_fixed_line", (se.basquin(sa, c, 1.0), n1)), ("V2_per_sample_line", (n2_0, n2)),
                          ("V3_recalibrated_fixed_line_DIAGNOSTIC", ((sa / (a3 * 1.0)) ** (1 / b), n3))):
            tab_rows.append(stats_row(c, var, "inputs+force only", nfun[0]))
            tab_rows.append(stats_row(c, var, "+ material scatter", nfun[1]))
            tab_rows.append(stats_row(c, var, "+ life scatter, log10 SD 0.2 (ASSUMED)", nfun[1] * 10 ** (LIFE_SD * eps)))
    pf = pd.DataFrame(pf_rows); pf.to_csv(OUT / "pf_limit_state_vs_life.csv", index=False)
    dis = pd.DataFrame(dis_rows); dis.to_csv(OUT / "classification_disagreement.csv", index=False)
    tab = pd.DataFrame(tab_rows); tab.to_csv(OUT / "table12_with_variants.csv", index=False)
    print("\n[2,3] Pf and supporting numbers"); print(pf.to_string(index=False, float_format=lambda v: f"{v:.5g}"))
    print("\nclassification disagreement vs g<0"); print(dis.to_string(index=False))
    print("\n[6] Table 12 with variants"); print(tab.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    # check V1 reproduces Stage E
    stage_e_life = pd.read_csv(p.RESULTS_DIR / "stageE" / "life_at_1200N.csv")
    print("\nStage E file check (V1 should match):")
    print(stage_e_life[stage_e_life.variant != "+ life scatter (ASSUMED)"][["criterion", "variant", "P_N_below_1e6"]].to_string(index=False))

    # ---------- 7. figure ----------
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    bins = np.linspace(2, 9, 141)
    for lab, n, col in (("V1 fixed line (Stage E)", lives["goodman"]["V1"], "#4c78a8"),
                        ("V2 per-sample line", lives["goodman"]["V2"], "#f58518")):
        ok = n > 0
        ax.hist(np.log10(n[ok]), bins, alpha=.55, color=col, label=f"{lab}: P(N<1e6) = {np.mean(n < 1e6):.3f}")
    ax.axvline(6, color="k", lw=1.2); ax.set_xlabel("log10 N (cycles), Goodman, 1200 N, material scatter ASSUMED")
    ax.set_ylabel("count (400,000 samples)"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "goodman_log10N_hist.png", dpi=600)


if __name__ == "__main__":
    main()
