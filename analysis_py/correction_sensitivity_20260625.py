#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sensitivity of the boundary consensus to the SEA-AD specimen-metadata
correction published 2026-06-25 on the Sage Bionetworks AD Knowledge Portal.

  Run:  python3 analysis_py/correction_sensitivity_20260625.py
  Out:  output/tables/SuppTable_S7_correction_sensitivity.csv
        output/tables/SuppTable_S7_correction_sensitivity_readme.txt

WHY THIS IS PYTHON AND NOT R
----------------------------
The 19-component engine is re-implemented here and VERIFIED against the two
reference outputs already committed to this repository:
  * output/tables/celltype_audit/multimethod_Astrocyte.csv  -> 19/19 exact
  * output/tables/celltype_audit/multimethod_Microglia.csv   -> 19/19 exact
  * output/tables/SuppTable_S2_consensus.csv                 -> 16/19 exact
The three components that differ in the last comparison are Segmented_PC1..3
and the cause is documented below (.HAVE_SEG). An R re-implementation should be
added and cross-checked once an R installation is available.

SEGMENTED COMPONENTS ARE EXCLUDED, FOLLOWING THE MANUSCRIPT
----------------------------------------------------------
The manuscript (STAR Methods, 'Boundary-detection engine') states: "Components
(9) depend on the segmented-regression implementation and are excluded from
every analysis reported here; the remaining 16 reproduce the deposited values
exactly." This script therefore reports 16-component consensus counts only.

The reason the manuscript excludes them is reproduced here for the record:
boundary_detection.R branches on `.HAVE_SEG <- requireNamespace("segmented")`.
With the package present, Segmented_PC1/2/3 = 0.665898 / 0.654355 / 0.416294;
without it, the 200-point fallback grid gives 0.835621 / 0.280384 / 0.243573.
Grid-point tests on each panel's own axis confirm which branch produced each
committed table: SuppTable_S2_consensus.csv was generated WITH the package,
multimethod_Astrocyte.csv and multimethod_Microglia.csv WITHOUT it.

OUTSTANDING REPOSITORY ACTION (not a manuscript defect)
-------------------------------------------------------
The deposited tables still carry 19 components and their region counts, while
the manuscript reports 16. Reviewers who regenerate from this repository will
find consensus counts that the manuscript does not state:
    SuppTable_S2_consensus.csv            19 components, Segmented included
    celltype_audit/generalization_summary 4/19 (astrocyte), 6/19 (microglia)
    celltype_audit/headtohead_summary     3/19 (astrocyte), on a different
                                          cell set: 70,009 astrocytes and
                                          40,000 microglia, against 67,419 in
                                          the deposited bin matrix
These tables should be regenerated on 16 components before deposition.

WHAT THE CORRECTION ACTUALLY SAYS
---------------------------------
Columns in the change log are "correct" then "org." (the original, wrong value).
Three MTG rows exist; only two can touch an MTG RNA-seq object:

  M2AT_220121_206_D01  donor        correct H21.33.047  org. H18.30.002  M/M  MTG/MTG
      -> M2AT = ATAC modality. Cannot affect an RNA-seq analysis. Its original
         (wrong) donor H18.30.002 is a Reference-series donor and is absent
         from the 84-donor cohort in any case.
  M2TX_210727_207_E01  donor;tissue correct H20.33.016/DFC  org. H21.33.037/MTG  F/F
  M1TX_210723_159_G01  donor;tissue correct H21.33.037/MTG  org. H20.33.016/DFC  F/F
      -> a clean two-way swap between two TRANSCRIPTOMIC specimens.

Consequence for a pre-correction MTG RNA object: the nuclei labelled
H21.33.037 are in fact the DFC tissue of donor H20.33.016, whose own CPS is
0.616531 (bin [0.6,0.7)); and H21.33.037's true MTG specimen was filed under
DFC and therefore excluded. Both specimens are female, so there is no sex
confound. CPS is a donor-level property, so the substituted specimen re-enters
bin 1 at the same axis position, CPS 0.180319.

STRUCTURE OF THE TEST
---------------------
Bins are FIXED CPS intervals, not donor quantiles, so removing a donor cannot
move a bin edge or move any other donor between bins; exactly one row of the
9 x 44 matrix changes. That row changes in two separable ways:

 (A) LOCATION -- exact, needs no expression data.
     centre' = (centre - w*CPS_d)/(1-w); every detector reports a bin centre or
     a midpoint of two adjacent centres, so the shift in the reported
     breakpoint is exact.
 (B) FRAGILITY -- the minimum perturbation radius r* over ALL directions in
     R^44 that moves the Ward primary call, inverted into the donor deviation
     ||mu_bin1 - mu_contaminant|| it would require.
 (C) MECHANISM -- the contaminant is a mid/late-trajectory sample sitting in
     the lowest bin, so removing it must move bin 1 AWAY from the bulk, whereas
     every breaking direction points TOWARD the bulk. Verified empirically by
     using each of the 8 other bins in turn as the contaminant's profile proxy.
"""
import csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from engine import detect_boundaries                        # noqa: E402
from scipy.cluster.hierarchy import linkage, fcluster       # noqa: E402
from scipy.spatial.distance import pdist                    # noqa: E402

PANELS = [
    ('44panel',             'data/SEAAD_44panel_bin_means.csv'),
    ('paper1panel',         'data/SEAAD_paper1panel_bin_means.csv'),
    ('celltype_Astrocyte',  'data/celltype_audit/binmeans_Astrocyte.csv'),
    ('celltype_Microglia',  'data/celltype_audit/binmeans_Microglia.csv'),
]
PRIMARY = 'Hier_k3_bp1'
SEGCOMP = {'Segmented_PC1', 'Segmented_PC2', 'Segmented_PC3'}
CPS_037, CPS_016 = 0.180318627999973, 0.616531

# H21.33.037's share of the first bin.
#   W_MEAS is the measured astrocyte share, 415 of the bin's 1,511 astrocytes
#   (p3fig/data/leave_one_donor.csv; the bin's four donors contribute
#   428 + 314 + 415 + 354 = 1,511, matching the deposited bin matrix exactly).
#   W_LP is the linear-programming upper bound derived from n_cells and
#   cps_center alone, kept for the microglia panel, where no per-donor cell
#   count is available.
#   The measured value also validates the closed-form location analysis used
#   here: at W_MEAS this script returns 0.204893, and the exhaustive
#   leave-one-donor-out run returns 0.204893 for the same donor.
W_MEAS, W_LP = 415 / 1511, 0.6667


def load(path):
    rows = list(csv.reader(open(os.path.join(REPO, path))))
    hdr = rows[0]; drop = {'cps_bin', 'cps_center', 'n_cells'}
    gi = [i for i, h in enumerate(hdr) if h not in drop]
    lab = [r[hdr.index('cps_bin')] for r in rows[1:]]
    ctr = np.array([float(r[hdr.index('cps_center')]) for r in rows[1:]])
    M = np.array([[float(r[i]) for i in gi] for r in rows[1:]])
    return M, ctr, lab


def primary_of(M, ctr):
    cl = fcluster(linkage(pdist(M, metric='euclidean'), method='ward'), t=3, criterion='maxclust')
    for i in range(1, len(cl)):
        if cl[i] != cl[i - 1]:
            return (ctr[i - 1] + ctr[i]) / 2.0
    return float('nan')


def consensus(d, p, branch_independent=False):
    it = {k: v for k, v in d.items() if k not in SEGCOMP} if branch_independent else d
    return sum(abs(v - round(p, 3)) <= 0.05 for v in it.values()), len(it)


def min_radius(M, ctr, row=0, ntrial=700, seed=11, rmax=16.0, nr=240):
    rng = np.random.default_rng(seed)
    base = primary_of(M, ctr)
    dirs = []
    for j in range(M.shape[0]):
        if j != row:
            v = M[j] - M[row]; n = np.linalg.norm(v)
            if n > 0:
                dirs += [v / n, -v / n]
    for _ in range(ntrial):
        v = rng.normal(size=M.shape[1]); dirs.append(v / np.linalg.norm(v))
    radii = np.concatenate([np.linspace(0.005, 2.0, nr // 2), np.linspace(2.0, rmax, nr - nr // 2)])
    best, bn = np.inf, float('nan')
    T = M.copy()
    for u in dirs:
        for r in radii:
            if r >= best:
                break
            T[row] = M[row] + r * u
            p = primary_of(T, ctr)
            if abs(p - base) >= 1e-12:
                best, bn = float(r), float(p)
                break
        T[row] = M[row]
    return best, bn


out = []
for name, path in PANELS:
    M, ctr, lab = load(path)
    base_all = detect_boundaries(M, ctr)
    base = base_all[PRIMARY]
    span = float(np.linalg.norm(M[-1] - M[0]))
    nn = float(sorted(np.sqrt(((M - M[0]) ** 2).sum(1))[1:])[0])
    r_star, breaks_to = min_radius(M, ctr)
    proxy = [i for i, L in enumerate(lab) if L in ('CPS0.65', '[0.6,0.7)')][0]

    for w, wtag in ((W_MEAS, 'measured_astrocyte_share_415_of_1511'), (W_LP, 'LP_provable_upper_bound')):
        c2 = ctr.copy(); c2[0] = (ctr[0] - w * CPS_037) / (1 - w)
        # (C) mechanism: contaminant profile = H20.33.016's own CPS bin
        d = w / (1 - w) * (M[0] - M[proxy])
        T = M.copy(); T[0] = M[0] + d
        got = detect_boundaries(T, c2)
        p = got[PRIMARY]
        c_bi, _ = consensus(got, p, branch_independent=True)
        b_bi, _ = consensus(base_all, base, branch_independent=True)
        # worst proxy over all 8 bins
        worst = 0.0
        for j in range(1, M.shape[0]):
            dj = w / (1 - w) * (M[0] - M[j])
            Tj = M.copy(); Tj[0] = M[0] + dj
            worst = max(worst, abs(primary_of(Tj, c2) - base))
        out.append(dict(
            panel=name, weight_scenario=wtag, w_bin1_cell_share=round(w, 6),
            published_primary_Ward=round(base, 9),
            corrected_primary_Ward=round(p, 9),
            shift=round(p - base, 9),
            published_rounded=round(base, 2), corrected_rounded=round(p, 2),
            consensus_published_16=f'{b_bi}/16', consensus_corrected_16=f'{c_bi}/16',
            bin1_centre_published=round(float(ctr[0]), 9),
            bin1_centre_corrected=round(float(c2[0]), 9),
            min_breaking_radius_r_star=round(r_star, 6),
            r_star_pct_of_nearest_row_gap=round(r_star / nn * 100, 2),
            required_donor_deviation_to_break=round(r_star * (1 - w) / w, 6),
            required_deviation_SD_per_marker=round(r_star * (1 - w) / w / np.sqrt(M.shape[1]), 4),
            trajectory_span=round(span, 6),
            max_abs_shift_over_all_8_profile_proxies=round(worst, 9),
            verdict='conclusion unchanged' if worst <= 0.02 else 'REVIEW',
        ))

OUT = os.path.join(REPO, 'output/tables')
os.makedirs(OUT, exist_ok=True)
f = os.path.join(OUT, 'SuppTable_S7_correction_sensitivity.csv')
with open(f, 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
    wr.writeheader(); wr.writerows(out)
print(f'wrote {f}  ({len(out)} rows)')
for r in out:
    print(f"  {r['panel']:20s} {r['weight_scenario']:30s} "
          f"{r['published_primary_Ward']:.6f} -> {r['corrected_primary_Ward']:.6f} "
          f"({r['shift']:+.6f})  worst-proxy |shift| {r['max_abs_shift_over_all_8_profile_proxies']:.6f}"
          f"  {r['verdict']}")

readme = os.path.join(OUT, 'SuppTable_S7_correction_sensitivity_readme.txt')
open(readme, 'w').write(__doc__.strip() + '\n')
print(f'wrote {readme}')
