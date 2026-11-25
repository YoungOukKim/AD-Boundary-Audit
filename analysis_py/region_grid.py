"""Grid bias across the nine SEA-AD brain regions, from published CPS alone.

The published claim: on the real SEA-AD grid, pure noise puts Ward's first
boundary at the first midpoint 100% of the time, because the first bin is
small and its mean is therefore noisy.

Here the same simulation is run on the donor-CPS grid of every brain region
for which Allen publishes a CPS (Travaglini 2026), plus the Gabitto 2024 axis
the manuscript used.  No expression data is needed: the claim is about the
grid, not about the genes.  Weighting is by DONORS, the paper's own unit of
replication (the published run weighted by nuclei, which are unavailable
outside MTG).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.cluster import KMeans

X = INPUT + "cps.xlsx"
CUTS = np.arange(0.2, 0.91, 0.1)
NG, B, SEED = 44, 1000, 42


def bin_of(v):
    return np.clip(np.digitize(np.asarray(v, float), CUTS, right=False), 0, 8)


def first_change(cl, ax):
    for i in range(1, len(cl)):
        if cl[i] != cl[i - 1]:
            return (ax[i - 1] + ax[i]) / 2.0
    return np.nan


def ward1(M, ax):
    return first_change(fcluster(linkage(M, method="ward"), 3, criterion="maxclust"), ax)


def km1(M, ax):
    lab = KMeans(n_clusters=3, n_init=10, random_state=SEED).fit_predict(M)
    return first_change(lab, ax)


def run(cps, tag):
    """cps: 1-D donor CPS values for one region/axis."""
    lab = bin_of(cps)
    keep = sorted(set(lab.tolist()))
    n = np.array([int((lab == b).sum()) for b in keep], float)
    ax = np.array([float(np.mean(cps[lab == b])) for b in keep])
    o = np.argsort(ax); n, ax = n[o], ax[o]
    if len(ax) < 4:
        return None
    firstmid = (ax[0] + ax[1]) / 2.0
    rng = np.random.default_rng(SEED)
    hits_w = hits_k = 0
    for _ in range(B):
        M = rng.standard_normal((len(ax), NG)) / np.sqrt(n)[:, None]
        hits_w += abs(ward1(M, ax) - firstmid) < 1e-9
        hits_k += abs(km1(M, ax) - firstmid) < 1e-9
    return dict(grid=tag, n_bins=len(ax), n_donors=int(n.sum()), first_bin_n=int(n[0]),
                first_mid=round(firstmid, 4), chance=round(1.0 / (len(ax) - 1), 4),
                ward_first_mid=hits_w / B, km_first_mid=hits_k / B,
                comp="·".join(str(int(v)) for v in n))


T = pd.read_excel(X, "Travaglini 2026")
G = pd.read_excel(X, "Gabitto 2024")

rows = []
r = run(G["CPS"].to_numpy(float), "MTG — Gabitto 2024 (원고가 쓴 축)")
if r: rows.append(r)
for col, nm in [("CPS_Global", "Global"), ("CPS_Local", "Local")]:
    for reg in ["MTG", "DFC", "V1C", "HIP", "STG", "MEC", "AnG", "ITG", "FI"]:
        s = T.loc[T["Brain Region"] == reg, col].to_numpy(float)
        r = run(s, f"{reg} — Travaglini {nm}")
        if r: rows.append(r)

out = pd.DataFrame(rows)
out.to_csv(DATA + "region_grid_bias.csv", index=False)
pd.set_option("display.width", 220)
print(out.to_string(index=False))
print()
tr = out[out.grid.str.contains("Travaglini")]
print(f"Travaglini 격자 {len(tr)}개 | Ward 첫 중점 적중률 "
      f"min {tr.ward_first_mid.min():.3f} / 중앙값 {tr.ward_first_mid.median():.3f} / max {tr.ward_first_mid.max():.3f}")
print(f"우연 기대치 범위 {tr.chance.min():.3f}–{tr.chance.max():.3f}")
print(f"적중률 0.90 이상인 격자: {(tr.ward_first_mid>=0.90).sum()} / {len(tr)}")
