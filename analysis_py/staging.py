"""Neuropathological staging axes against the CPS grid.

Writes data/pathology_axes.csv, data/pathology_peraxis.csv and
data/pathology_engine.csv. Nothing else: the marker-panel tables belong to
final.py and conditions.py, and no two scripts in this directory write the
same file.

Two conventions appear here and the difference is deliberate.

  joint run (pathology_engine.csv)
      The four neuropathological axes are placed on exactly the grid the
      44-marker panel receives: nucleus-weighted bin means, nucleus-weighted
      bin-centre CPS. This is what makes the two analyses comparable, which
      is the whole point of running it.

  per-axis run (pathology_peraxis.csv)
      Staging is a donor-level quantity, so each axis is summarised by the
      unweighted donor average of its bin and the axis is the donor-mean CPS
      of that bin. Reported separately in Results and labelled there.

pathology_axes.csv carries the Spearman correlations against CPS and the
nucleus-weighted bin means, i.e. the joint-run inputs.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import sys, pathlib, numpy as np, pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from engine import detect_boundaries, piecewise_fit, dp_changepoint

D, E = DATA, INPUT
CUTS = np.arange(0.2, 0.91, 0.1)      # nine bins: <0.2, 0.2-0.3, ..., >=0.9
BP, TOL = 0.206829687177018, 0.05     # the panel's Ward first boundary
WIN = (0.40, 0.60)                    # the window the atlas authors report
SEED, NPERM = 42, 1000                # as everywhere else in this pipeline

ORD = {
    "braak": {"Braak 0": 0, "Braak I": 1, "Braak II": 2, "Braak III": 3,
              "Braak IV": 4, "Braak V": 5, "Braak VI": 6},
    "thal": {f"Thal {i}": i for i in range(6)},
    "cerad": {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3},
    "adnc": {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3},
    "cogstat": {"No dementia": 0, "Dementia": 1},
    "lewy": {"Not Identified (olfactory bulb not assessed)": 0,
             "Not Identified (olfactory bulb assessed)": 0, "Olfactory bulb only": 1,
             "Amygdala-predominant": 2, "Brainstem-predominant": 2,
             "Limbic (Transitional)": 3, "Neocortical (Diffuse)": 4},
    "late": {"Not Identified": 0, "LATE Stage 1": 1, "LATE Stage 2": 2,
             "LATE Stage 3": 3},
}
JOINT = ["braak", "thal", "cerad", "adnc"]

if __name__ == "__main__":
    mn = pd.read_csv(E + "donor_mean_Astrocyte.csv")
    mn = mn[~mn.cps.isna()].reset_index(drop=True)
    DONOR = "donor" if "donor" in mn.columns else "donor_id"
    NC, CPS = mn.n_cells.to_numpy(float), mn.cps.to_numpy(float)
    LAB = np.clip(np.digitize(CPS, CUTS, right=False), 0, 8)
    BIDX = [np.where(LAB == b)[0] for b in range(9)]
    meta = pd.read_csv(E + "donor_meta.csv").set_index("donor_id")

    # ---- 1. correlations and the nucleus-weighted bin means -----------------
    rows, mat = [], {}
    for a in ORD:
        if a not in meta.columns:
            continue
        y = meta[a].reindex(mn[DONOR]).map(ORD[a]).to_numpy(float)
        if np.isnan(y).all():
            continue
        ok = ~np.isnan(y)
        rho, p = spearmanr(CPS[ok], y[ok])
        row, bmean = {"axis": a, "rho_cps": rho, "p_cps": p, "n": int(ok.sum())}, []
        for b in range(9):
            k = BIDX[b][ok[BIDX[b]]]
            v = (y[k] * NC[k]).sum() / NC[k].sum() if len(k) else np.nan
            row[f"bin{b + 1}"] = v
            bmean.append(v)
        rows.append(row)
        if a in JOINT:
            mat[a] = bmean
    PA = pd.DataFrame(rows)
    PA.to_csv(D + "pathology_axes.csv", index=False)
    print("pathology ~ CPS")
    print(PA[["axis", "rho_cps", "p_cps", "n"]].round(4).to_string(index=False))

    # ---- 2. per-axis breakpoints, donor-level convention --------------------
    ax_d = np.array([CPS[b].mean() for b in BIDX])
    per = []
    for a in ORD:
        if a not in meta.columns:
            continue
        y = meta[a].reindex(mn[DONOR]).map(ORD[a]).to_numpy(float)
        if np.isnan(y).all():
            continue
        ok = ~np.isnan(y)
        v = np.array([y[b[ok[b]]].mean() if len(b[ok[b]]) else np.nan for b in BIDX])
        if np.isnan(v).any():
            continue
        per.append({"axis": a, "piecewise": float(piecewise_fit(ax_d, v)),
                    "dp": float(dp_changepoint(ax_d, v))})
    PX = pd.DataFrame(per)
    PX.to_csv(D + "pathology_peraxis.csv", index=False)
    print("\nper-axis breakpoints (donor-mean axis, unweighted donor averaging)")
    print(PX.round(4).to_string(index=False))

    # ---- 3. joint run on the panel's own grid -------------------------------
    Mp = np.array([mat[a] for a in JOINT], float).T
    axp = np.array([(CPS[b] * NC[b]).sum() / NC[b].sum() for b in BIDX])
    bd = {k: float(v) for k, v in detect_boundaries(Mp, axp).items()
          if not k.startswith("Segmented") and v == v}
    pd.DataFrame([{"component": k, "breakpoint": v} for k, v in bd.items()]
                 ).sort_values("breakpoint").to_csv(D + "pathology_engine.csv", index=False)
    # ---- 4. null for the staging boundary ----------------------------------
    # The question Figure 5B asks is whether 0.388 is unusual for ANY four-column
    # matrix read on this grid, so the null is 1,000 draws of four background
    # genes put through the same bin averaging and the same engine, with no
    # z-scoring, exactly as the staging matrix above is built. Its median is the
    # grid's forced position, which is the point of the comparison.
    GEN = [g for g in mn.columns
           if g not in (DONOR, "donor", "donor_id", "cps", "n_cells")
           and mn[g].std() > 0]
    V = mn[GEN].to_numpy(float)
    rng = np.random.default_rng(SEED)
    pn = []
    for _ in range(NPERM):
        X = V[:, rng.choice(len(GEN), Mp.shape[1], replace=False)]
        Mn = np.array([[(X[b, j] * NC[b]).sum() / NC[b].sum()
                        for j in range(X.shape[1])] for b in BIDX])
        pn.append(detect_boundaries(Mn, axp).get("Hier_k3_bp1", np.nan))
    pn = np.array(pn, float)
    pd.DataFrame({"bp": pn}).to_csv(D + "pathology_null.csv", index=False)
    p_stage = float(np.nansum(np.abs(pn - bd["Hier_k3_bp1"]) <= TOL) / NPERM)
    print(f"staging null: median {np.nanmedian(pn):.4f} | "
          f"p({bd['Hier_k3_bp1']:.3f}) = {p_stage:.3f}")

    vv = np.array(list(bd.values()))
    print(f"\npathology engine: Ward bp1 = {bd['Hier_k3_bp1']:.4f}; "
          f"{int(((vv >= WIN[0]) & (vv <= WIN[1])).sum())} of {len(vv)} in "
          f"{WIN[0]:.2f}-{WIN[1]:.2f}; "
          f"{int((np.abs(vv - BP) <= TOL).sum())} of {len(vv)} within {TOL} of {BP:.3f}")
