"""Phase 3b: does the model-form share survive matched factor ranges?

Same variance decomposition as model_form_variance.py, with levels chosen so that the notch factor and the mean tensile
strength move the endurance limit by comparable amounts:
  matched A: Kf = Kt x (0.90, 1.00, 1.10); mean UTS = 900 x (0.90, 1.00, 1.10)
  matched B: Kf = Kt x (0.90, 1.00, 1.10); mean UTS = 900 x (0.85, 1.00, 1.15)
Criterion (Goodman, Gerber, SWT) and the scatter set (low, base, high) are as before.
Writes results/phase3/model_form_ranges.json
"""
import itertools
import json

import numpy as np
import pandas as pd

import fatigue as ft
import model_form_variance as mv
import params as p
import sensitivity as sens
import stage_d as sd
import uncertainty as unc

OUT = p.RESULTS_DIR / "phase3"


def run(gp, kf_scale, uts_scale):
    rng = np.random.default_rng(7)
    u = rng.random((mv.N_MC, len(sd.NAMES)))
    peak = sd.peak_of(gp, u)
    rows = []
    for sname, (rsd, cv) in mv.SCATTER.items():
        _, uts0, ratio = sd.physical(u, cv, rsd, 0.0)
        for us in uts_scale:
            uts = uts0 * us
            for crit, ks in itertools.product(mv.CRIT, kf_scale):
                g = ft.g_fun(peak, uts, ratio, crit, ft.KT_B * ks)
                rows.append({"scatter": sname, "uts_mean": round(900 * us, 1), "criterion": crit,
                             "kf": f"{ks:.2f}Kt", "Pf": float(np.mean(g < 0))})
    return pd.DataFrame(rows)


def main():
    df = pd.read_csv(p.RESULTS_DIR / "stageA" / "fe_samples_al.csv")
    gp, _, _ = sens.fit_surrogate(df[unc.FE_NAMES].to_numpy(), df["peak_vm"].to_numpy())
    fac = ["criterion", "kf", "uts_mean", "scatter"]
    out = {}
    for label, uts_scale in (("matched A (Kf +-10%, UTS +-10%)", (0.9, 1.0, 1.1)),
                             ("matched B (Kf +-10%, UTS +-15%)", (0.85, 1.0, 1.15))):
        d = run(gp, (0.9, 1.0, 1.1), uts_scale)
        d.to_csv(OUT / f"Pf_{label.split()[1]}.csv", index=False)
        out[label] = mv.summarise(d, label, fac)
    json.dump(out, open(OUT / "model_form_ranges.json", "w"), indent=2)
    for k, s in out.items():
        print(k, {x: round(v, 3) for x, v in s.items() if isinstance(v, float)})
        for f, v in s["anova"].items():
            print("   ", f, round(v["S1"], 3), round(v["ST"], 3))


if __name__ == "__main__":
    main()
