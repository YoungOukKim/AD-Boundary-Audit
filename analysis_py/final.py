"""Every number in the revised paper, from one script.

The observed panel is the authors' own published bin matrix, taken verbatim from
the repository (inputs/SEAAD_44panel_bin_means.csv = AD-Boundary-Audit
data/SEAAD_44panel_bin_means.csv). That removes every preprocessing question
from the observed side. Background panels for the null are rebuilt from the
donor-level extract with the repository's own convention: z-score each gene
across all astrocyte nuclei, then take the bin mean.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import numpy as np, pandas as pd
from engine import detect_boundaries
from scipy.stats import spearmanr, kstest
from scipy.cluster.hierarchy import linkage, fcluster

D, E, A = DATA, INPUT, INPUT
SEG = ("Segmented_PC1", "Segmented_PC2", "Segmented_PC3")
TOL, NPERM, NIND, SEED = 0.05, 1000, 300, 42
FLOOR, AMB_HI = 0.02, 0.176
AMBIENT = ["P2RY12", "PLP1", "RBFOX3"]

# ---- observed: the authors' published matrix --------------------------------
pub = pd.read_csv(INPUT + "SEAAD_44panel_bin_means.csv")
PANEL = [c for c in pub.columns if c not in ("cps_bin", "cps_center", "n_cells")]
AX = pub.cps_center.to_numpy(float)
MOBS = pub[PANEL].to_numpy(float)
NCELL = pub.n_cells.to_numpy(int)

# ---- donor-level extract, on the repository's z-scoring convention -----------
mn = pd.read_csv(E + "donor_mean_Astrocyte.csv")
sq = pd.read_csv(E + "donor_sumsq_Astrocyte.csv")
GEN = [c for c in mn.columns if c not in ("donor_id", "n_cells", "cps")]
Nall = mn.n_cells.to_numpy(float); Ntot = Nall.sum()
MU = (mn[GEN].to_numpy(float) * Nall[:, None]).sum(0) / Ntot
VAR = (sq[GEN].to_numpy(float).sum(0) - Ntot * MU ** 2) / (Ntot - 1)
SD = np.sqrt(np.maximum(VAR, 1e-300)); SD[SD <= 0] = 1.0
GI = {g: i for i, g in enumerate(GEN)}

m = mn[~mn.cps.isna()].reset_index(drop=True)
DONOR = "donor_id"
NK = m.n_cells.to_numpy(float); CPS = m.cps.to_numpy(float)
LAB = np.clip(np.digitize(CPS, np.arange(0.2, 0.91, 0.1), right=False), 0, 8)
BIDX = [np.where(LAB == b)[0] for b in range(9)]
ALL = np.arange(len(m))

dd = pd.read_csv(A + "detection_by_donor_Astrocyte.csv")
m["ambient"] = m[DONOR].map(pd.Series(dd[AMBIENT].mean(axis=1).to_numpy(), index=dd.donor_id))
assert m.ambient.notna().all()

det = pd.read_csv(D + "detection_astrocyte.csv"); det = det[det.in_panel == True].copy()
def verdict(r):
    if r.detection < FLOOR:         return "below floor"
    if r.donors_passing_floor < 60: return "too few donors"
    if r.detection <= AMB_HI:       return "at/below ambient"
    return "observable"
det["v"] = det.apply(verdict, axis=1)
OBSV = [g for g in det[det.v == "observable"].gene if g in GEN]

SKIP = {"donor_id", "n_cells", "cps", "ambient"}
BG = [g for g in GEN if g not in SKIP and g not in set(PANEL) and mn[g].to_numpy().std() > 0]
LOGMU = {g: np.log10(mn[g].to_numpy().mean() + 1e-6) for g in BG + PANEL}
BGS = sorted(BG, key=lambda g: LOGMU[g]); BGV = np.array([LOGMU[g] for g in BGS])


def zmat(cols, rows=None, X=None):
    """bin matrix on the repository's convention, for an arbitrary gene set."""
    idx = [GI[g] for g in cols]
    Z = (m[cols].to_numpy(float) - MU[idx]) / SD[idx] if X is None else X
    keep = np.zeros(len(m), bool); keep[ALL if rows is None else rows] = True
    out, ax = [], []
    for b in range(9):
        k = BIDX[b][keep[BIDX[b]]]
        if len(k) == 0: continue
        w = NK[k]
        out.append((Z[k] * w[:, None]).sum(0) / w.sum())
        ax.append((CPS[k] * w).sum() / w.sum())
    return np.array(out, float), np.array(ax)


def comps(M, ax, drop_seg=True):
    bd = detect_boundaries(M, ax)
    return {k: float(v) for k, v in bd.items()
            if v == v and not (drop_seg and k in SEG)}


def ward1(M, ax):
    if len(ax) < 4: return np.nan
    cl = fcluster(linkage(M, method="ward"), 3, criterion="maxclust")
    for i in range(1, len(cl)):
        if cl[i] != cl[i - 1]: return (ax[i - 1] + ax[i]) / 2
    return np.nan


def matched(rng, target):
    out, used = [], set()
    for g in target:
        j = int(np.searchsorted(BGV, LOGMU[g]))
        lo, hi = max(0, j - 25), min(len(BGS), j + 25)
        cand = [x for x in BGS[lo:hi] if x not in used] or [x for x in BGS if x not in used]
        pick = cand[int(rng.integers(len(cand)))]
        used.add(pick); out.append(pick)
    return out


if __name__ == "__main__":
    rng = np.random.default_rng(SEED)
    OBS = comps(MOBS, AX)
    BP = OBS["Hier_k3_bp1"]
    sup = int(np.sum(np.abs(np.array(list(OBS.values())) - BP) <= TOL))
    print(f"observed (published matrix): Ward bp1 = {BP:.5f}  support {sup}/{len(OBS)}")
    pd.DataFrame([{"component": k, "breakpoint": v} for k, v in sorted(OBS.items(), key=lambda t: t[1])]
                 ).to_csv(D + "observed_components.csv", index=False)

    # ---- 1. the same boundary from our own reconstruction -------------------
    Mrec, axr = zmat(PANEL)
    REC = comps(Mrec, axr)
    diff = np.abs(Mrec - MOBS)
    frag = [(k, OBS[k], REC.get(k, np.nan)) for k in OBS]
    nmove = sum(1 for _, a, b in frag if abs(a - b) > TOL)
    print(f"reconstruction: max cell |diff| {diff.max():.4f}, mean {diff.mean():.4f}; "
          f"{nmove} of {len(OBS)} components move > {TOL}")
    pd.DataFrame(frag, columns=["component", "published_matrix", "reconstructed"]
                 ).sort_values("component").to_csv(
                     D + "reconstruction_fragility.csv", index=False)
    for k, a, b in frag:
        if abs(a - b) > TOL: print(f"   {k:24s} {a:.4f} -> {b:.4f}")

    # ---- 2. specificity, like-for-like and as published ---------------------
    wnull, pnull = [], []
    for _ in range(NPERM):
        cols = matched(rng, PANEL)
        M2, a2 = zmat(cols)
        bd = comps(M2, a2)
        wnull.append(bd.get("Hier_k3_bp1", np.nan))
        pnull.append(bd.get("PCA_PC1", np.nan))
    wnull = np.array(wnull, float); pnull = np.array(pnull, float)
    f = lambda n, x: int(np.nansum(np.abs(n - x) <= TOL))
    pc1_obs = OBS["PCA_PC1"]
    print(f"\nspecificity, 1000 expression-matched panels")
    print(f"  Ward {BP:.3f} vs PC1 null  (as published): {f(pnull,BP):4d}/1000  p={f(pnull,BP)/1000:.3f}")
    print(f"  Ward {BP:.3f} vs Ward null (like for like): {f(wnull,BP):4d}/1000  p={f(wnull,BP)/1000:.3f}")
    print(f"  PC1  {pc1_obs:.3f} vs PC1 null             : {f(pnull,pc1_obs):4d}/1000  p={f(pnull,pc1_obs)/1000:.3f}")
    print(f"  null medians: Ward {np.nanmedian(wnull):.3f} | PC1 {np.nanmedian(pnull):.3f}")
    pd.DataFrame({"ward_null": wnull, "pc1_null": pnull}).to_csv(D + "nulls.csv", index=False)

    # ---- 3. restriction to the observable markers ---------------------------
    Mo, axo = zmat(OBSV)
    OBSO = comps(Mo, axo); BPo = OBSO["Hier_k3_bp1"]
    onull = []
    for _ in range(NPERM):
        M2, a2 = zmat(matched(rng, OBSV))
        onull.append(comps(M2, a2).get("Hier_k3_bp1", np.nan))
    onull = np.array(onull, float)
    print(f"\nrestriction: {len(OBSV)} observable markers -> Ward bp1 {BPo:.4f}  "
          f"(all 44: {BP:.4f}); Ward-null {f(onull,BPo)}/1000  p={f(onull,BPo)/1000:.3f}")
    pd.DataFrame([{"n_panel": len(PANEL), "n_observable": len(OBSV),
                   "n_below_floor": int((det.v == "below floor").sum()),
                   "n_too_few": int((det.v == "too few donors").sum()),
                   "n_ambient": int((det.v == "at/below ambient").sum()),
                   "bp_all": BP, "bp_observable": BPo,
                   "p_ward_all": f(wnull, BP) / 1000, "p_ward_obs": f(onull, BPo) / 1000}]
                 ).to_csv(D + "restriction.csv", index=False)

    # ---- 4. component independence ------------------------------------------
    panels = []
    for _ in range(NIND):
        M2, a2 = zmat(matched(rng, PANEL))
        panels.append(comps(M2, a2))
    keys = sorted({k for b in panels for k in b})
    Mx = pd.DataFrame([{k: b.get(k, np.nan) for k in keys} for b in panels])
    pairs = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            ok = Mx[[a, b]].dropna()
            if len(ok) < 30: continue
            pairs.append({"a": a, "b": b, "identical": float((ok[a] - ok[b]).abs().lt(1e-9).mean())})
    P = pd.DataFrame(pairs).sort_values("identical", ascending=False)
    nd = Mx.apply(lambda c: c.dropna().round(6).nunique()).rename("n_distinct")
    eff = Mx.apply(lambda r: r.dropna().round(6).nunique(), axis=1)
    P.to_csv(D + "independence.csv", index=False); nd.to_csv(D + "independence_distinct.csv")
    print(f"\nindependence over {NIND} matched panels")
    print(P.head(4).to_string(index=False))
    print(f"  near-constant components: {nd[nd<6].to_dict()}")
    print(f"  median distinct values per panel: {eff.median():.0f} of {len(keys)}")

    # ---- 5. donor dimension --------------------------------------------------
    lodo = [{"donor": m[DONOR].iloc[i], "cps": CPS[i], "n_astro": int(NK[i]),
             "bp": ward1(*zmat(PANEL, rows=np.delete(ALL, i)))} for i in range(len(m))]
    L = pd.DataFrame(lodo); L.to_csv(D + "leave_one_donor.csv", index=False)
    moved = L[(L.bp - BP).abs() > TOL]
    bs = [ward1(*zmat(PANEL, rows=rng.choice(ALL, len(ALL), replace=True))) for _ in range(NPERM)]
    B = pd.DataFrame({"bp": bs}); B.to_csv(D + "bootstrap.csv", index=False)
    bv = B.bp.dropna().to_numpy(); lo, hi = np.percentile(bv, [2.5, 97.5])
    print(f"\nleave-one-donor: {len(moved)} of {len(L)} move > {TOL}")
    print(f"bootstrap: median {np.median(bv):.3f}  95% CI {lo:.3f}-{hi:.3f}  "
          f"within {TOL}: {np.mean(np.abs(bv-BP)<=TOL)*100:.1f}%")
    ks = kstest(CPS, "uniform")
    dpb = [len(k) for k in BIDX]
    print(f"CPS uniform? KS D={ks.statistic:.3f} p={ks.pvalue:.2e} | donors/bin {dpb} | nuclei/bin {NCELL.tolist()}")
    pd.DataFrame({"bin": range(1, 10), "cps_centre": AX, "n_donors": dpb,
                  "n_cells": NCELL}).to_csv(D + "bin_composition.csv", index=False)
    print("WROTE observed_components, nulls, restriction, independence, leave_one_donor, bootstrap, bin_composition")
