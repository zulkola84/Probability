"""Stage E: (1) probabilistic Basquin life at 1200 N, (2) design chart: thickness for beta = 2, 3, 4 at 800-1400 N.

(1) Basquin, Bologna et al. 2023 Table 2 (values checked against the paper's HTML text; the equation itself was not
    readable): N = (sigma_a / A)^(1/b), A = 670.53 / 549.36 MPa, b = -0.088 / -0.059 (Goodman / Gerber).
    sigma_a = 0.45 * peak / KT_B (nominal alternating stress at R = 0.1, peak from the FE surrogate): this reading is
    ASSUMED and is tested against the paper's ~575,000 (analytical) and ~1e6 (FE) cycles in `basquin_check.csv`.
    Scatter: FE inputs + force (uncertainty.py), material scatter as in Stage D "assumed" (UTS CV 3 %, ratio SD 0.05;
    scales A by ratio*UTS/450, ASSUMED), optional life scatter about the line, log10(N) SD 0.2 (ASSUMED; no source).
(2) Limit state of fatigue.py (Goodman, Gerber, Kf = Kt). Peak stress at design thickness t and load F:
    peak = GP(x) * R(t) * F/1200 with R(t) = FE(t)/FE(5.12) from a nominal-input FE sweep (`fe_thickness.csv`);
    tolerance scatter of the thickness is kept in the GP input. Width and hole size stay at the calibrated values.
Writes results/stageE/*. Usage: python stage_e.py
"""
import tempfile  # noqa: F401
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

import fatigue as ft
import params as p
import run_samples
import sensitivity as sens
import uncertainty as unc

OUT = p.RESULTS_DIR / "stageE"
TABLE2 = {"goodman": (670.53, -0.088), "gerber": (549.36, -0.059)}
N_FE = len(unc.FE_NAMES)
T_GRID = [4.0, 4.5, 5.12, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0]
LOADS = list(range(800, 1401, 100))
BETAS = (2, 3, 4)
LIFE_SD_LOG10 = 0.2       # ASSUMED


def basquin(sa, crit, scale=1.0):
    A, b = TABLE2[crit]
    return (sa / (A * scale)) ** (1 / b)


def fe_nominal(t):
    x = [d.mean() for d in unc.FE_INPUTS.values()]
    x[0] = t
    return run_samples.run_one(x)[0]


def fe_point(x):
    return run_samples.run_one(x)[0]


def thickness_sweep():
    path = OUT / "fe_thickness.csv"
    if path.exists():
        return pd.read_csv(path)
    with Pool(8) as pool:
        pk = pool.map(fe_nominal, T_GRID)
    df = pd.DataFrame({"thickness": T_GRID, "peak_vm_1200N": pk})
    df.to_csv(path, index=False)
    return df


def validate_scaling(gp, rfun, ts=(6.0, 7.5), n=24):
    """FE at random inputs with thickness shifted to ~t, vs GP(x at 5.12 + dev) * R(t)."""
    rng = np.random.default_rng(5)
    rows = []
    x = unc.to_physical(rng.random((n, N_FE)))
    base = ft.peak_stress(gp, np.column_stack([x, np.full(n, p.RUNOUT_FORCE)]))
    for t in ts:
        x2 = x.copy(); x2[:, 0] += t - p.THICKNESS
        with Pool(8) as pool:
            fe = np.array(pool.map(fe_point, list(x2)))
        pred = base * rfun(t)
        err = 100 * (pred - fe) / fe
        rows.append({"thickness": t, "n": n, "bias_pct": err.mean(), "rmse_pct": np.sqrt(np.mean(err**2)),
                     "max_abs_pct": np.abs(err).max()})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(exist_ok=True)
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    gp, q2, rmse = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
    print(f"surrogate: {len(df)} runs, Q2 {q2:.4f}, RMSE {rmse:.2f}")

    # ---------- (1) life ----------
    rows = []
    for lab, pk in (("analytical", p.PAPER_PEAK_ANALYTICAL), ("FE", p.PAPER_PEAK_FE)):
        sa = (1 - ft.R) / 2 * pk / ft.KT_B
        n = {c: basquin(sa, c) for c in TABLE2}
        rows.append({"peak_source": lab, "peak_MPa": pk, "sigma_a": sa, "N_goodman": n["goodman"], "N_gerber": n["gerber"],
                     "N_mean": np.mean(list(n.values())), "N_geomean": np.sqrt(n["goodman"] * n["gerber"]),
                     "paper_N": 575_000 if lab == "analytical" else 1_000_000})
    chk = pd.DataFrame(rows); chk["mean_vs_paper_pct"] = 100 * (chk.N_mean / chk.paper_N - 1)
    chk.to_csv(OUT / "basquin_check.csv", index=False); print(chk.round(3).to_string(index=False))

    rng = np.random.default_rng(3)
    nmc = 400_000
    x = ft.to_physical(rng.random((nmc, len(ft.NAMES))))
    peak = ft.peak_stress(gp, x)
    sa = (1 - ft.R) / 2 * peak / ft.KT_B
    scale = x[:, -1] * x[:, -2] / (0.5 * 900.0)
    eps = rng.standard_normal(nmc)
    rows = []
    for crit in TABLE2:
        for var, sc, ls in (("inputs+force only", 1.0, 0.0), ("+ material scatter", scale, 0.0),
                            ("+ life scatter (ASSUMED)", scale, LIFE_SD_LOG10)):
            n = basquin(sa, crit, sc) * 10 ** (ls * eps)
            q = np.percentile(n, [5, 50, 95])
            rows.append({"criterion": crit, "variant": var, "N_p05": q[0], "N_median": q[1], "N_p95": q[2],
                         "P_N_below_1e6": np.mean(n < 1e6)})
    life = pd.DataFrame(rows); life.to_csv(OUT / "life_at_1200N.csv", index=False)
    print(life.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    # ---------- (2) design chart ----------
    sweep = thickness_sweep()
    ref = float(sweep.loc[sweep.thickness == 5.12, "peak_vm_1200N"].iloc[0])
    rfun = PchipInterpolator(sweep.thickness, np.log(sweep.peak_vm_1200N / ref))
    R = lambda t: float(np.exp(rfun(t)))
    val = validate_scaling(gp, R); val.to_csv(OUT / "thickness_scaling_check.csv", index=False); print(val.round(3).to_string(index=False))
    print("FE nominal peak vs thickness:", dict(zip(sweep.thickness, sweep.peak_vm_1200N.round(1))))

    n_d = 2_000_000
    rng = np.random.default_rng(9)
    xd = ft.to_physical(rng.random((n_d, len(ft.NAMES))))
    peak1200 = ft.peak_stress(gp, np.column_stack([xd[:, :N_FE], np.full(n_d, p.RUNOUT_FORCE) * (xd[:, N_FE] / p.RUNOUT_FORCE)]))
    zf = xd[:, N_FE] / p.RUNOUT_FORCE                    # force sample / 1200 (relative scatter about 1)
    peak_unit = peak1200 / zf                            # peak at 1200 N, scatter in force removed
    uts, ratio = xd[:, -2], xd[:, -1]

    def pf(t, load, crit):
        pk = peak_unit * R(t) * load / p.RUNOUT_FORCE * zf     # force scatter (CV 2 %) kept, mean = design load
        return np.mean(ft.g_fun(pk, uts, ratio, crit) < 0)

    rows = []
    for crit in ("goodman", "gerber"):
        for load in LOADS:
            for beta in BETAS:
                target = stats.norm.cdf(-beta)
                f = lambda t: pf(t, load, crit) - target
                lo, hi = 4.0, 12.0
                if f(lo) < 0:
                    t_req, note = lo, "<= 4.0 mm (grid lower end)"
                elif f(hi) > 0:
                    t_req, note = np.nan, "> 12 mm"
                else:
                    t_req, note = brentq(f, lo, hi, xtol=1e-3), ""
                rows.append({"criterion": crit, "load_N": load, "beta": beta, "t_required_mm": t_req,
                             "Pf_check": pf(t_req, load, crit) if np.isfinite(t_req) else np.nan, "note": note,
                             "Pf_at_5.12": pf(5.12, load, crit)})
        print(crit, "done")
    chart = pd.DataFrame(rows); chart.to_csv(OUT / "design_chart.csv", index=False)
    print(chart.pivot_table(index=["criterion", "load_N"], columns="beta", values="t_required_mm").round(2).to_string())

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for a, crit in zip(ax, ("goodman", "gerber")):
        for beta, col in zip(BETAS, ("#4c78a8", "#f58518", "#54a24b")):
            d = chart[(chart.criterion == crit) & (chart.beta == beta)]
            a.plot(d.load_N, d.t_required_mm, "-o", color=col, label=f"beta = {beta}")
        a.axhline(5.12, color="k", ls="--", lw=1, label="paper plate 5.12 mm")
        a.set_title(f"{crit.capitalize()} limit state (runout at 1e6)"); a.set_xlabel("applied load (N)"); a.legend()
    ax[0].set_ylabel("required thickness (mm)")
    fig.tight_layout(); fig.savefig(OUT / "design_chart.png", dpi=150)


if __name__ == "__main__":
    main()
