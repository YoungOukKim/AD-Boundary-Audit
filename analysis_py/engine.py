#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Faithful Python port of the 9-algorithm / 19-component boundary-detection engine
embedded verbatim in R/analysis/celltype_generalization_audit.R
(= R/analysis/boundary_detection.R).

Validated against the committed reference output
  output/tables/celltype_audit/multimethod_{Astrocyte,Microglia}.csv
to 15 significant digits before use.

Port notes (why this is exact):
  * every PC-based detector is SSE- or abs()-based, so the arbitrary sign of a
    principal component does not affect any breakpoint.
  * prcomp(center=TRUE, scale.=FALSE) == column-mean-centering + SVD.
  * scipy linkage(method='ward') == R hclust(method='ward.D2').
  * .cluster_boundaries only reads label CHANGES between adjacent rows, so
    cluster label permutation is irrelevant.
  * the reference output used the .HAVE_SEG == FALSE fallback grid
    (proved by Segmented_PC3 landing exactly on axis_vals[2]).
"""
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist


def _lm_sse(x, y):
    """SSE of a simple OLS fit y ~ x (intercept + slope)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    X = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    return float(r @ r)


def piecewise_fit(x, y):
    """R: piecewise_fit() -- 50-point grid from x[2]+0.02 to x[n-1]-0.02 (1-indexed)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    if n < 4 or np.nanstd(y, ddof=1) < 1e-10:
        return np.nan
    best_bp, best_sse = np.nan, np.inf
    # R seq(x[2]+0.02, x[length(x)-1]-0.02, length.out=50); R 1-indexed
    lo = x[1] + 0.02
    hi = x[n - 2] - 0.02
    for bp in np.linspace(lo, hi, 50):
        left = x <= bp
        right = x > bp
        if left.sum() < 2 or right.sum() < 2:
            continue
        sse = _lm_sse(x[left], y[left]) + _lm_sse(x[right], y[right])
        if sse < best_sse:
            best_sse, best_bp = sse, bp
    return best_bp


def cusum_bp(v, x):
    """R: x[which.max(abs(cumsum(v-mean(v))))]"""
    v = np.asarray(v, float)
    c = np.cumsum(v - v.mean())
    return float(np.asarray(x, float)[int(np.argmax(np.abs(c)))])


def dp_changepoint(x, y, min_seg=2):
    """R: dp_changepoint(); k over (min_seg+1):(n-min_seg), 1-indexed."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    bi, bs = None, np.inf
    for k in range(min_seg + 1, n - min_seg + 1):          # R k = 3..n-2 inclusive
        s = _lm_sse(x[:k], y[:k]) + _lm_sse(x[k:], y[k:])  # R 1:k and (k+1):n
        if s < bs:
            bs, bi = s, k
    if bi is None:
        return np.nan
    return float((x[bi - 1] + x[bi]) / 2.0)                # R (x[bi]+x[bi+1])/2


def _cluster_boundaries(cl, a):
    """R: .cluster_boundaries() -- midpoints wherever the adjacent label changes."""
    out = []
    for i in range(1, len(cl)):
        if cl[i] != cl[i - 1]:
            out.append((a[i - 1] + a[i]) / 2.0)
    return out


def _segmented_fallback(x, y):
    """R .HAVE_SEG==FALSE branch: 200-point grid, hinge basis y ~ x1 + x2."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    bs, ps = np.inf, np.nan
    for kk in np.linspace(x[1], x[n - 2], 200):            # R axis_vals[2] .. axis_vals[n-1]
        x1 = np.minimum(x, kk)
        x2 = np.maximum(x - kk, 0.0)
        X = np.column_stack([np.ones_like(x), x1, x2])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ beta
        s = float(r @ r)
        if s < bs:
            bs, ps = s, kk
    return ps


def detect_boundaries(expr_bin, axis_vals, kmeans_labels=None):
    """
    Returns dict {algorithm: breakpoint}. 19 components.
    kmeans_labels: optional precomputed k=3 labels (to control the RNG path).
    """
    M = np.asarray(expr_bin, float)
    a = np.asarray(axis_vals, float)
    assert M.shape[0] == len(a)
    res = {}

    # --- PCA (prcomp center=TRUE scale=FALSE) ---
    Mc = M - M.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
    scores = U * S                                          # == pca$x
    npc = min(3, scores.shape[1])

    for k in range(npc):
        bp = piecewise_fit(a, scores[:, k])
        if not np.isnan(bp):
            res[f"PCA_PC{k+1}"] = bp
    for k in range(npc):
        res[f"CUSUM_PC{k+1}"] = cusum_bp(scores[:, k], a)

    # --- max consecutive euclidean gap ---
    gaps = np.array([np.sqrt(((M[i] - M[i - 1]) ** 2).sum()) for i in range(1, M.shape[0])])
    gp = (a[:-1] + a[1:]) / 2.0
    res["Mahal_max_gap"] = float(gp[int(np.argmax(gaps))])

    # --- Ward.D2 hierarchical, k=3  (PRIMARY) ---
    Z = linkage(pdist(M, metric="euclidean"), method="ward")
    cl = fcluster(Z, t=3, criterion="maxclust")
    b3 = _cluster_boundaries(cl, a)
    if len(b3) >= 1:
        res["Hier_k3_bp1"] = b3[0]
    if len(b3) >= 2:
        res["Hier_k3_bp2"] = b3[1]

    # --- k-means k=3 ---
    if kmeans_labels is None:
        kmeans_labels = _kmeans_global(M, 3)
    bk = _cluster_boundaries(kmeans_labels, a)
    if len(bk) >= 1:
        res["KMeans_k3_bp1"] = bk[0]
    if len(bk) >= 2:
        res["KMeans_k3_bp2"] = bk[1]

    # --- dynamic-programming changepoint ---
    for k in range(npc):
        bp = dp_changepoint(a, scores[:, k])
        if not np.isnan(bp):
            res[f"DP_PC{k+1}"] = bp

    # --- variance jump (row variances, ddof=1) ---
    rv = M.var(axis=1, ddof=1)
    vd = np.diff(rv)
    if len(vd) > 0:
        j = int(np.argmax(np.abs(vd)))
        res["Variance_jump"] = float((a[j] + a[j + 1]) / 2.0)

    # --- trajectory inflection (2nd difference of distance from bin 1) ---
    if len(a) >= 4:
        d0 = np.array([np.sqrt(((M[i] - M[0]) ** 2).sum()) for i in range(M.shape[0])])
        res["Trajectory_inflection"] = float(a[int(np.argmax(np.abs(np.diff(np.diff(d0))))) + 1])

    # --- segmented (fallback grid, as used for the reference output) ---
    for k in range(npc):
        res[f"Segmented_PC{k+1}"] = _segmented_fallback(a, scores[:, k])

    return res


def _kmeans_global(M, k):
    """
    Exhaustive global optimum for contiguous-agnostic k-means on <=12 points.
    R used kmeans(nstart=10, seed 42); with 9 points and k=3 Hartigan-Wong
    with 10 restarts reaches the global optimum, so the global optimum is the
    faithful target and is reproducible without R's RNG.
    """
    from itertools import product
    n = M.shape[0]
    best, bl = np.inf, None
    for lab in product(range(k), repeat=n):
        lab = np.array(lab)
        if len(set(lab.tolist())) != k:
            continue
        w = 0.0
        for c in range(k):
            P = M[lab == c]
            w += ((P - P.mean(axis=0)) ** 2).sum()
        if w < best:
            best, bl = w, lab
    return bl
