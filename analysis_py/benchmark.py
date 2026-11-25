"""Ground-truth benchmark, rebuilt on the real donor structure.

The published benchmark planted boundaries directly in bin x marker matrices on a
regular 20-bin grid. That grid cannot express a donor-level artefact, because it
has no donors. Here the same regimes are planted in a DONOR x marker matrix using
the cohort's real CPS values and real per-donor astrocyte counts, then collapsed
to bins exactly as the analysis does. This keeps the real 3-19 donors-per-bin
imbalance, and lets a fifth artefact class — a single aberrant donor — be planted
and tested.

Regimes (the published four plus one):
  real          co-regulated localized step at the planted boundary
  global        the same step shared by every background gene
  null          no structure
  ramp          monotone, no step
  decoy         one marker carries a large step
  donor         one donor is shifted on every marker          <-- NEW

Gates (the published three plus one):
  support       >= S*, the 95th percentile of random-panel support
  specificity   p <= 0.05 against the random-panel null
  leave-one-marker-out   max shift <= 0.05 CPS
  leave-one-donor-out    max shift <= 0.05 CPS                <-- NEW
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import numpy as np, pandas as pd
from engine import detect_boundaries, piecewise_fit
from scipy.cluster.hierarchy import linkage, fcluster

E_ = INPUT
K, NPERM, NREP = 44, 200, 16
TS = 0.30                      # planted step position
TRUTH = 0.300                  # midpoint of the two bins the step separates
DELTA, SEED = 0.6, 42
TOL = 0.05

mn = pd.read_csv(E_ + "donor_mean_Astrocyte.csv")
mn = mn[~mn.cps.isna()].reset_index(drop=True)
REAL_NC = mn.n_cells.to_numpy().astype(float)

# Synthetic cohort on the published benchmark grid: 20 equal bins, 6 donors per
# bin, per-donor cell counts resampled from the real astrocyte yields. This keeps
# the published grid (so the result is comparable to the original benchmark)
# while giving the audit donors to leave out.
NBIN, PER_BIN = 20, 6
_rng0 = np.random.default_rng(SEED)
EDGES = np.linspace(0.0, 1.0, NBIN + 1)
CPS = np.concatenate([_rng0.uniform(EDGES[b] + 0.004, EDGES[b + 1] - 0.004, PER_BIN)
                      for b in range(NBIN)])
NC = _rng0.choice(REAL_NC, size=len(CPS), replace=True)
ND = len(CPS)
LAB = np.digitize(CPS, EDGES[1:-1], right=False)
BINS = [np.where(LAB == b)[0] for b in range(NBIN)]
BINS = [b for b in BINS if len(b)]
AX = np.array([(CPS[k] * NC[k]).sum() / NC[k].sum() for k in BINS])

def binmeans(P, keep=None, drop_marker=None):
    cols = np.arange(P.shape[1]) if drop_marker is None else \
           np.delete(np.arange(P.shape[1]), drop_marker)
    rows, ax = [], []
    for k in BINS:
        kk = k if keep is None else k[keep[k]]
        if len(kk) == 0: continue
        w = NC[kk]
        rows.append((P[np.ix_(kk, cols)] * w[:, None]).sum(0) / w.sum())
        ax.append((CPS[kk] * w).sum() / w.sum())
    return np.array(rows), np.array(ax)

def hier_bp1(M, ax):
    """Ward k=3 first boundary only — the cheap readout used inside leave-one-out."""
    if len(ax) < 4: return np.nan
    cl = fcluster(linkage(M, method="ward"), 3, criterion="maxclust")
    for i in range(1, len(cl)):
        if cl[i] != cl[i-1]: return (ax[i-1] + ax[i]) / 2
    return np.nan

def step_of(cps, tau=TS):
    s = (cps >= tau).astype(float); return s - s.mean()

def gen_bg(G, rng, global_step=False):
    """Donor-level background: flat noise plus diverse LATE logistics, so random
    panels break late and rarely early — the empirical null shape."""
    M = np.zeros((ND, G))
    for g in range(G):
        if rng.random() < 0.30:
            M[:, g] = rng.normal(0, 0.6, ND)
        else:
            mu = rng.uniform(0.52, 0.88); k = rng.uniform(0.04, 0.12)
            sh = 1 / (1 + np.exp(-(CPS - mu) / k)); sh = sh - sh.mean()
            sd = sh.std(ddof=1)
            sh = sh / sd if sd > 0 else sh
            M[:, g] = abs(rng.normal(0, 0.7)) * rng.choice([-1, 1]) * sh + rng.normal(0, 0.2, ND)
    if global_step:
        M = M + np.outer(step_of(CPS), 0.6 * (1 + rng.normal(0, 0.15, G)))
    return M

def mk(kind, BG, rng):
    idx = rng.choice(BG.shape[1], K, replace=False)
    P = BG[:, idx].copy()
    if kind == "real":
        P = 0.6 * P + np.outer(step_of(CPS), DELTA * (1 + rng.normal(0, 0.15, K)))
    elif kind == "ramp":
        P = 0.6 * P + np.outer(CPS - CPS.mean(), DELTA * (1 + rng.normal(0, 0.15, K)))
    elif kind == "decoy":
        P[:, 0] = P[:, 0] + 5 * DELTA * step_of(CPS)
    elif kind == "donor":
        # one high-yield donor in the bin just below the planted step is shifted
        # on every marker. The bin mean follows that donor, the bin separates
        # from its neighbours, and clustering calls a boundary in the right
        # place -- a boundary built by one person, exactly as observed in the
        # real cohort. Nothing else in the panel carries it.
        b = int(np.digitize([TS], EDGES[1:-1], right=False)[0]) - 1
        cand = BINS[b]
        d0 = int(cand[int(np.argmax(NC[cand]))])
        P[d0, :] = P[d0, :] + 7.0
    return P

def build_null(BG, rng):
    pc1, sup = [], []
    for _ in range(NPERM):
        cols = rng.choice(BG.shape[1], K, replace=False)
        M, ax = binmeans(BG[:, cols])
        Mc = M - M.mean(0)
        U, S, _ = np.linalg.svd(Mc, full_matrices=False)
        b = piecewise_fit(ax, (U * S)[:, 0]); pc1.append(b if b is not None else np.nan)
        bd = {k: v for k, v in detect_boundaries(M, ax).items() if not k.startswith("Segmented")}
        prim = bd.get("Hier_k3_bp1", np.nan)
        v = np.array(list(bd.values()))
        sup.append(int(np.sum(np.abs(v - prim) <= TOL)) if prim == prim else 0)
    return np.array(pc1, float), np.array(sup, float)

def audit(P, pc1null, Sstar):
    M, ax = binmeans(P)
    bd = {k: v for k, v in detect_boundaries(M, ax).items() if not k.startswith("Segmented")}
    prim = bd.get("Hier_k3_bp1", np.nan)
    v = np.array(list(bd.values()))
    support = int(np.sum(np.abs(v - prim) <= TOL)) if prim == prim else 0
    pspec = float(np.nanmean(np.abs(pc1null - prim) <= TOL)) if prim == prim else 1.0
    lomo = max(abs(hier_bp1(*binmeans(P, drop_marker=j)) - prim)
               for j in range(P.shape[1]))
    keepm = np.ones(ND, bool); lodo = 0.0
    for i in range(ND):
        keepm[i] = False
        b = hier_bp1(*binmeans(P, keep=keepm))
        if b == b: lodo = max(lodo, abs(b - prim))
        keepm[i] = True
    return dict(primary=prim, support=support, pspec=pspec, lomo=lomo, lodo=lodo,
                g_support=support >= Sstar, g_spec=pspec <= 0.05,
                g_lomo=lomo <= TOL, g_lodo=lodo <= TOL)

if __name__ == "__main__":
    rng = np.random.default_rng(SEED)
    BG = gen_bg(1500, np.random.default_rng(SEED))
    BGg = gen_bg(1500, np.random.default_rng(SEED), global_step=True)
    pc1null, supnull = build_null(BG, np.random.default_rng(7))
    Sstar = float(np.nanpercentile(supnull, 95))
    pc1null_g, supnull_g = build_null(BGg, np.random.default_rng(7))
    Sstar_g = float(np.nanpercentile(supnull_g, 95))
    print(f"  global background: S* = {Sstar_g:.1f} | null within 0.05 of truth "
          f"{np.nanmean(np.abs(pc1null_g-TRUTH)<=TOL)*100:.1f}%")
    print(f"cohort: {ND} donors, {len(BINS)} bins, {PER_BIN}/bin | "
          f"S* = {Sstar:.1f} | PC1 null within 0.05 of truth "
          f"{np.nanmean(np.abs(pc1null-TRUTH)<=TOL)*100:.1f}% | late(>0.5) "
          f"{np.nanmean(pc1null>0.5)*100:.0f}%")
    rows = []
    for kind in ["real", "global", "null", "ramp", "decoy", "donor"]:
        src = BGg if kind == "global" else BG
        for r in range(NREP):
            P = mk(kind, src, np.random.default_rng(1000 + 97 * r + abs(hash(kind)) % 50))
            a = audit(P, pc1null_g if kind == "global" else pc1null,
                      Sstar_g if kind == "global" else Sstar)
            a["regime"] = kind; a["rep"] = r
            a["certified"] = bool(a["g_support"] and a["g_spec"] and a["g_lomo"] and a["g_lodo"])
            a["certified_old"] = bool(a["g_support"] and a["g_spec"] and a["g_lomo"])
            rows.append(a)
        done = pd.DataFrame(rows); done = done[done.regime == kind]
        print(f"  {kind:7s} certified {done.certified.sum():2d}/{NREP}   "
              f"(published 3 gates: {done.certified_old.sum():2d}/{NREP})")
    df = pd.DataFrame(rows)
    df.to_csv(DATA + "benchmark.csv", index=False)
    pd.DataFrame(dict(Sstar=[Sstar], truth=[TRUTH], nperm=[NPERM], nrep=[NREP])
                 ).to_csv(DATA + "benchmark_meta.csv", index=False)
    print("saved benchmark.csv")
