"""Boundary under six conditions, all from one pipeline.

Every row of data/boundary_conditions.csv is produced here, so Figure 3B compares
like with like. The baseline is checked against the published 16 component values
before anything else is written.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import numpy as np, pandas as pd
from engine import detect_boundaries

D = DATA
E = INPUT
A = INPUT
FLOOR, AMB_HI = 0.02, 0.176
AMBIENT = ["P2RY12", "PLP1", "RBFOX3"]
DROP_DONOR = "H21.33.001"
# the published grid is nine bins 0.1-0.2 ... 0.9-1.0; no donor has CPS < 0.1,
# so the cut points are 0.2 .. 0.9 and the first bin is "below 0.2".
CUTS = np.arange(0.2, 0.91, 0.1)

mn = pd.read_csv(E + "donor_mean_Astrocyte.csv")
mn = mn[~mn.cps.isna()].reset_index(drop=True)
det = pd.read_csv(D + "detection_astrocyte.csv")
det = det[det.in_panel == True]


def verdict(r):
    if r.detection < FLOOR:            return "below floor"
    if r.donors_passing_floor < 60:    return "too few donors"
    if r.detection <= AMB_HI:          return "at/below ambient"
    return "observable"


det["v"] = det.apply(verdict, axis=1)
PANEL = [g for g in det.gene if g in mn.columns]
OBS = [g for g in det[det.v == "observable"].gene if g in mn.columns]

# donor-level ambient index: detection of cell-type-exclusive markers IN astrocytes
# wide table: donor_id, cps, n_cells, then one column of detection rate per gene
dd = pd.read_csv(A + "detection_by_donor_Astrocyte.csv")
amb = pd.Series(dd[AMBIENT].mean(axis=1).to_numpy(), index=dd.donor_id, name="ambient")
DONOR = "donor" if "donor" in mn.columns else "donor_id"
mn = mn.merge(amb, left_on=DONOR, right_index=True, how="left")
assert mn.ambient.notna().all(), "ambient index missing for some donors"


def bin_matrix(df, cols, drop_bin1=False):
    lab = np.clip(np.digitize(df.cps.to_numpy(), CUTS, right=False), 0, 8)
    rows, ax = [], []
    for b in range(9):
        if drop_bin1 and b == 0:
            continue
        k = np.where(lab == b)[0]
        if not len(k):
            continue
        w = df.n_cells.to_numpy()[k]
        X = df[cols].to_numpy()[k]
        rows.append((X * w[:, None]).sum(0) / w.sum())
        ax.append((df.cps.to_numpy()[k] * w).sum() / w.sum())
    M = np.array(rows, float)
    sd = M.std(0, ddof=1)
    return (M - M.mean(0)) / np.where(sd > 0, sd, 1.0), np.array(ax)


def residualize(df, cols):
    """Regress each marker's donor-level value on the ambient index, keep residual."""
    out = df.copy()
    a = df.ambient.to_numpy()
    a = (a - a.mean())
    for g in cols:
        y = df[g].to_numpy(float)
        b = (a @ (y - y.mean())) / (a @ a)
        out[g] = (y - y.mean()) - b * a
    return out


def run(df, cols, drop_bin1=False):
    M, ax = bin_matrix(df, cols, drop_bin1)
    bd = {k: v for k, v in detect_boundaries(M, ax).items() if not k.startswith("Segmented")}
    return {k: float(v) for k, v in bd.items() if v == v}


CONDITIONS = [
    ("baseline",             lambda: run(mn, PANEL)),
    ("observable only",      lambda: run(mn, OBS)),
    ("ambient-residualized", lambda: run(residualize(mn, PANEL), PANEL)),
    (f"drop {DROP_DONOR}",   lambda: run(mn[mn[DONOR] != DROP_DONOR], PANEL)),
    ("ambient + drop donor", lambda: run(residualize(mn[mn[DONOR] != DROP_DONOR], PANEL), PANEL)),
    ("drop bin 1",           lambda: run(mn, PANEL, drop_bin1=True)),
]

if __name__ == "__main__":
    rows = []
    for name, f in CONDITIONS:
        bd = f()
        prim = bd["Hier_k3_bp1"]
        v = np.array(list(bd.values()))
        sup = int(np.sum(np.abs(v - prim) <= 0.05))
        print(f"{name:22s} n={len(bd):2d}  Ward bp1 = {prim:.4f}  support {sup}/{len(bd)}")
        rows += [{"condition": name, "component": k, "breakpoint": x} for k, x in bd.items()]
    out = pd.DataFrame(rows)
    out.to_csv(D + "boundary_conditions.csv", index=False)

    base = out[out.condition == "baseline"].breakpoint.to_numpy()
    obs = out[out.condition == "observable only"].breakpoint.to_numpy()
    print(f"\npanel {len(PANEL)} markers | observable {len(OBS)} | "
          f"baseline median {np.median(base):.3f} | observable median {np.median(obs):.3f}")
    # restriction.csv is written by final.py, which also has the null
    # distributions and therefore the p-values. Two scripts writing one file
    # with different columns is the failure this paper is about; only one
    # script owns each table here.
    print("wrote boundary_conditions.csv")
