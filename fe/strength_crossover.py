"""Phase 2, step 8: strength-scatter crossover in the Sobol ranking of g.

Reuses stage_d.physical (UTS and endurance-ratio samples, Gaussian copula for the correlation), the Stage D 5-input GP on all
416 FE runs, fatigue.g_fun (Kf = Kt = 1.66, UTS mean 900 MPa, ratio mean 0.5, width fixed), 400,000 common-random-number samples
(seed 1, same stream as Stage D) and a Saltelli base sample of 8192 (8 variables). Grid:
  main grid  : endurance-ratio SD 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08 (ASSUMED) x UTS CV 0.5 % (ASSUMED),
               1.1 % (Okazaki 2012, doi:10.3390/ma5122981), 3 % (baseline, ASSUMED), 5 % (ASSUMED),
               8.5 % (Suzuki 2018, doi:10.4012/dmj.2017-125, component-level proxy), correlation 0 (35 cases)
  correlated : the 7 ratio SDs with UTS CV 3 % and correlation 0.7 (ASSUMED, as the brief); the Sobol 'uts' index is then the
               SHARED strength-driver index and 'endurance_ratio' the independent residual (as Table 11 / stage_d.py)
  scan       : UTS CV 0.5, 1, 1.5, 2, 2.5, 3, 4, 5 % (ASSUMED) x ratio SD 0.0001, 0.005, 0.01, 0.02, 0.03, 0.05 (ASSUMED), correlation 0,
               to locate the UTS CV at which the strength group equals the geometric/load maximum
  extension  : ratio SD 0.0025, 0.001, 0.0003, 0.0001 (ASSUMED), correlation 0, every UTS CV, to find crossovers below the grid
Strength group = ST(endurance_ratio) + ST(uts). Geometric/load group member = max ST over thickness, force, lever_arm, hole_d.
Outputs go to OUT (outside results/). Usage: python strength_crossover.py
"""
import itertools
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SALib.analyze import sobol
from SALib.sample import sobol as sobol_sample
from scipy import stats

import fatigue as ft
import params as p
import sensitivity as sens
import stage_d as sd_mod
import uncertainty as unc

OUT = Path(__file__).resolve().parents[1] / "results" / "phase2" / "step8_strength_crossover"
NAMES = ft.NAMES
CRITS = ("goodman", "gerber", "swt")
RATIO_SD = (0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08)           # ASSUMED
RATIO_SD_EXT = (0.0025, 0.001, 0.0003, 0.0001)                   # ASSUMED (extension below the grid)
UTS_CV = {0.01: "ASSUMED", 0.015: "ASSUMED", 0.02: "ASSUMED", 0.025: "ASSUMED", 0.04: "ASSUMED", 0.005: "ASSUMED", 0.011: "Okazaki 2012 (doi:10.3390/ma5122981)", 0.03: "baseline, ASSUMED", 0.05: "ASSUMED",
          0.085: "Suzuki 2018 (doi:10.4012/dmj.2017-125), component proxy"}
MAIN_CV = (0.005, 0.011, 0.03, 0.05, 0.085)
GEOM = ["thickness", "force", "lever_arm", "hole_d"]
SCAN_CV = (0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05)     # ASSUMED
SCAN_SD = (0.0001, 0.005, 0.01, 0.02, 0.03, 0.05)                # ASSUMED
N_MC, N_SOBOL = 400_000, 2**13


def form(gp, crit, cv, sd, rho, iters=100):
    def g_u(z):
        x, uts, ratio = sd_mod.physical(stats.norm.cdf(z)[None, :], cv, sd, rho)
        return float(ft.g_fun(ft.peak_stress(gp, x), uts, ratio, crit)[0])
    z, done = np.zeros(len(NAMES)), False
    for _ in range(iters):
        g0 = g_u(z)
        grad = np.array([(g_u(z + 1e-4 * e) - g0) / 1e-4 for e in np.eye(len(z))])
        if grad @ grad < 1e-20:
            break
        z_new = (grad @ z - g0) / (grad @ grad) * grad
        done = np.linalg.norm(z_new - z) < 1e-6
        z = z_new
        if done:
            break
    beta = np.linalg.norm(z)
    if g_u(np.zeros(len(z))) < 0:
        return stats.norm.cdf(beta), -beta, done
    return stats.norm.cdf(-beta), beta, done


def crossing(sds, diff):
    """Largest-to-smallest scan: first SD (descending) where diff > 0; linear interpolation in log SD between it and the
    previous (larger) SD. Returns (sd, how) with how in {'interpolated', 'exceeded at the largest SD', 'none'}."""
    o = np.argsort(sds)[::-1]
    s, d = np.asarray(sds)[o], np.asarray(diff)[o]
    if d[0] > 0:
        return np.nan, "exceeded at the largest SD"
    for i in range(1, len(s)):
        if d[i] > 0:
            f = d[i - 1] / (d[i - 1] - d[i])
            return float(np.exp(np.log(s[i - 1]) + f * (np.log(s[i]) - np.log(s[i - 1])))), "interpolated"
    return np.nan, "none"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    gp, q2, rmse = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
    print(f"GP on {len(df)} runs, Q2 {q2:.4f}, CV RMSE {rmse:.2f} MPa")
    rng = np.random.default_rng(1)
    u = rng.random((N_MC, len(NAMES)))
    peak = sd_mod.peak_of(gp, u)
    prob = {"num_vars": len(NAMES), "names": NAMES, "bounds": [[0, 1]] * len(NAMES)}
    us = sobol_sample.sample(prob, N_SOBOL, calc_second_order=False, seed=1)
    peak_s = sd_mod.peak_of(gp, us)

    cases = [("main", cv, sd, 0.0) for cv, sd in itertools.product(MAIN_CV, RATIO_SD)]
    cases += [("correlated", 0.03, sd, 0.7) for sd in RATIO_SD]
    cases += [("extension", cv, sd, 0.0) for cv, sd in itertools.product(MAIN_CV, RATIO_SD_EXT)]
    cases += [("scan", cv, sd, 0.0) for cv, sd in itertools.product(SCAN_CV, SCAN_SD)]
    rows = []
    for k, (grp, cv, sd, rho) in enumerate(cases):
        _, uts, ratio = sd_mod.physical(u, cv, sd, rho)
        _, uts_s, ratio_s = sd_mod.physical(us, cv, sd, rho)
        for crit in CRITS:
            f = ft.g_fun(peak, uts, ratio, crit) < 0
            pf = f.mean()
            pff, beta, ok = form(gp, crit, cv, sd, rho)
            s = sobol.analyze(prob, ft.g_fun(peak_s, uts_s, ratio_s, crit), calc_second_order=False, seed=1)
            st = dict(zip(NAMES, s["ST"])); conf = dict(zip(NAMES, s["ST_conf"]))
            geo = {n: st[n] for n in GEOM}
            gname = max(geo, key=geo.get)
            row = {"group": grp, "ratio_sd": sd, "ratio_cv_pct": 100 * sd / 0.5, "uts_cv_pct": 100 * cv, "uts_cv_source": UTS_CV[cv],
                   "rho": rho, "criterion": crit, "Pf": pf, "mc_se": np.sqrt(pf * (1 - pf) / N_MC), "Pf_FORM": pff, "beta": beta,
                   "FORM_converged": ok}
            row.update({f"ST_{n}": st[n] for n in NAMES}); row.update({f"STconf_{n}": conf[n] for n in NAMES})
            row.update({"strength_group_ST": st["endurance_ratio"] + st["uts"], "geom_load_max_ST": geo[gname], "geom_load_max_name": gname,
                        "shared_driver_ST_uts": st["uts"] if rho > 0 else np.nan})
            rows.append(row)
        if (k + 1) % 10 == 0:
            print(f"  {k + 1}/{len(cases)} cases done")
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "crossover_grid_results.csv", index=False)
    pd.set_option("display.width", 250)
    print(res[res.group != "correlated"][["group", "ratio_sd", "uts_cv_pct", "criterion", "Pf", "mc_se", "Pf_FORM", "ST_endurance_ratio", "ST_uts",
                                          "strength_group_ST", "geom_load_max_ST", "geom_load_max_name"]].round(4).to_string(index=False))
    print(res[res.group == "correlated"][["ratio_sd", "criterion", "Pf", "Pf_FORM", "ST_endurance_ratio", "ST_uts", "strength_group_ST",
                                          "geom_load_max_ST", "geom_load_max_name"]].round(4).to_string(index=False))

    # ---------- crossovers (rho = 0; grid + extension)
    r0 = res[(res.rho == 0) & res.group.isin(["main", "extension"]) & res.uts_cv_pct.isin([100 * c for c in (0.005, 0.011, 0.03, 0.05, 0.085)])]
    xr = []
    for crit in CRITS:
        for cvp in sorted(r0.uts_cv_pct.unique()):
            d = r0[(r0.criterion == crit) & (r0.uts_cv_pct == cvp)].sort_values("ratio_sd")
            for what, a in (("strength group", d.strength_group_ST), ("endurance ratio alone", d.ST_endurance_ratio)):
                sdx, how = crossing(d.ratio_sd.values, (d.geom_load_max_ST - a).values)
                in_grid = bool(np.isfinite(sdx) and sdx >= min(RATIO_SD) - 1e-12)
                xr.append({"criterion": crit, "uts_cv_pct": cvp, "compared_with": what, "crossover_ratio_sd": sdx,
                           "crossover_ratio_cv_pct": 100 * sdx / 0.5 if np.isfinite(sdx) else np.nan, "how": how,
                           "inside_main_grid": in_grid if how == "interpolated" else False,
                           "smallest_sd_tested": float(d.ratio_sd.min()),
                           "geom_load_ST_at_smallest_sd": float(d.sort_values("ratio_sd").geom_load_max_ST.iloc[0]),
                           "strength_group_ST_at_smallest_sd": float(d.sort_values("ratio_sd").strength_group_ST.iloc[0])})
            # strength group crosses 0.5 ("most of the variance")
            sdh, howh = crossing(d.ratio_sd.values, (0.5 - d.strength_group_ST).values)
            xr.append({"criterion": crit, "uts_cv_pct": cvp, "compared_with": "strength group < 0.5 (fails 'most of the variance')",
                       "crossover_ratio_sd": sdh, "crossover_ratio_cv_pct": 100 * sdh / 0.5 if np.isfinite(sdh) else np.nan, "how": howh})
    xr = pd.DataFrame(xr); xr.to_csv(OUT / "crossover_table.csv", index=False)
    print(xr.round(5).to_string(index=False))

    # UTS CV at which the strength group equals the geometric/load maximum (rho = 0), at fixed ratio SD
    sc = res[(res.rho == 0) & res.group.isin(["main", "extension", "scan"])]
    br = []
    for crit in CRITS:
        for sdv in SCAN_SD:
            d = sc[(sc.criterion == crit) & np.isclose(sc.ratio_sd, sdv)].drop_duplicates("uts_cv_pct").sort_values("uts_cv_pct")
            diff = (d.strength_group_ST - d.geom_load_max_ST).values; cvp = d.uts_cv_pct.values
            val, how = np.nan, "none (strength group larger at all tested UTS CV)" if (diff > 0).all() else "none (geometric/load larger at all tested UTS CV)"
            for i in range(1, len(cvp)):
                if diff[i - 1] < 0 <= diff[i]:
                    val, how = float(cvp[i - 1] + diff[i - 1] / (diff[i - 1] - diff[i]) * (cvp[i] - cvp[i - 1])), "interpolated (linear in CV)"
                    break
            br.append({"criterion": crit, "ratio_sd": sdv, "ratio_cv_pct": 100 * sdv / 0.5, "uts_cv_pct_at_crossover": val, "how": how,
                       "n_uts_cv_points": len(d)})
    br = pd.DataFrame(br); br.to_csv(OUT / "uts_cv_crossover_by_ratio_sd.csv", index=False); print(br.round(3).to_string(index=False))

    # does the claim hold in every case of the main grid and correlated rows?
    mg = res[res.group.isin(["main", "correlated"])]
    claim = mg.groupby(["criterion"]).apply(lambda d: pd.Series({
        "cases": len(d), "cases_strength_group_ST_ge_0.5": int((d.strength_group_ST >= 0.5).sum()),
        "cases_strength_group_gt_geom_load": int((d.strength_group_ST > d.geom_load_max_ST).sum()),
        "min_strength_group_ST": d.strength_group_ST.min()}))
    claim.to_csv(OUT / "claim_check.csv"); print(claim.to_string())
    failing = mg[(mg.strength_group_ST < 0.5)][["group", "criterion", "ratio_sd", "ratio_cv_pct", "uts_cv_pct", "rho", "strength_group_ST", "geom_load_max_ST"]]
    failing.to_csv(OUT / "cases_where_claim_fails.csv", index=False); print(failing.round(4).to_string(index=False))
    # Goodman/Gerber order of Pf
    pv = res.pivot_table(index=["group", "ratio_sd", "uts_cv_pct", "rho"], columns="criterion", values="Pf").reset_index()
    pv["goodman_gt_gerber"] = pv.goodman > pv.gerber
    pv["swt_ge_goodman"] = pv.swt >= pv.goodman
    pv.to_csv(OUT / "pf_order_check.csv", index=False)
    print("Goodman > Gerber in", int(pv.goodman_gt_gerber.sum()), "of", len(pv), "cases; SWT >= Goodman in", int(pv.swt_ge_goodman.sum()))
    print(pv[~pv.goodman_gt_gerber].round(4).to_string(index=False))

    # ---------- figures
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    sel = res[(res.rho == 0) & (res.uts_cv_pct == 3.0)]
    colors = {"endurance_ratio": "#d62728", "uts": "#ff7f0e", "thickness": "#1f77b4", "force": "#2ca02c", "lever_arm": "#9467bd", "hole_d": "#8c564b"}
    for a, crit in zip(axs, CRITS):
        d = sel[sel.criterion == crit].sort_values("ratio_sd")
        for n, c in colors.items():
            a.plot(d.ratio_sd, d[f"ST_{n}"], "-o", ms=3, color=c, label=n)
        a.plot(d.ratio_sd, d.strength_group_ST, "k--", lw=1, label="strength group")
        a.axvspan(0.03, 0.08, color="0.85", zorder=0)
        a.set_xscale("log"); a.set_xlabel("endurance-ratio SD (UTS CV 3 %)"); a.set_title(crit.capitalize() if crit != "swt" else "SWT")
    axs[0].set_ylabel("total-order Sobol index of g"); axs[0].legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(OUT / "fig_sobol_vs_ratio_sd.png", dpi=600); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    cvs = sorted(res[res.group == "main"].uts_cv_pct.unique())
    for a, crit in zip(axs, CRITS):
        d = res[(res.group == "main") & (res.criterion == crit)]
        S = d.pivot(index="uts_cv_pct", columns="ratio_sd", values="strength_group_ST").loc[cvs]
        G = d.pivot(index="uts_cv_pct", columns="ratio_sd", values="geom_load_max_ST").loc[cvs]
        im = a.imshow(S.values, origin="lower", cmap="viridis", vmin=0, vmax=1, aspect="auto")
        for i in range(S.shape[0]):
            for j in range(S.shape[1]):
                a.text(j, i, f"{S.values[i, j]:.2f}", ha="center", va="center", fontsize=6, color="w" if S.values[i, j] < 0.6 else "k")
        diff = (S - G).values
        if diff.min() < 0 < diff.max():
            a.contour(diff, levels=[0], colors="r", linewidths=1.5)
        a.set_xticks(range(S.shape[1])); a.set_xticklabels([f"{c:g}" for c in S.columns], fontsize=7)
        a.set_yticks(range(S.shape[0])); a.set_yticklabels([f"{c:g}" for c in S.index], fontsize=7)
        a.set_xlabel("endurance-ratio SD"); a.set_title(crit.capitalize() if crit != "swt" else "SWT", fontsize=9)
    axs[0].set_ylabel("UTS CV (%)")
    fig.colorbar(im, ax=axs, shrink=0.85, label="strength-group index ST(ratio) + ST(UTS); red line: = largest geometric/load index")
    fig.savefig(OUT / "fig_strength_group_heatmap.png", dpi=600, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
