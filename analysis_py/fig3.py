"""Figure 3 - what the boundary does survive: observability and the donors."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import numpy as np, pandas as pd
from style import setup, OI, INK, INK2, MUTED, GRID, panel_label, ygrid, xgrid
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

setup()
D = DATA
bc = pd.read_csv(D + "bin_composition.csv")
det = pd.read_csv(D + "detection_astrocyte.csv"); det = det[det.in_panel == True].copy()
res = pd.read_csv(D + "restriction.csv").iloc[0]
L = pd.read_csv(D + "leave_one_donor.csv")
B = pd.read_csv(D + "bootstrap.csv")
BP, TOL = float(res.bp_all), 0.05
FLOOR, AMB_LO, AMB_HI = 0.02, 0.134, 0.176

fig = plt.figure(figsize=(7.3, 5.3))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.46,
                      width_ratios=[1.0, 1.35, 1.0])

# ---- (a) donors and nuclei per bin ------------------------------------------
ax = fig.add_subplot(gs[0, 0]); panel_label(ax, "a", dx=-0.30)
x = np.arange(9)
ax.bar(x, bc.n_donors, 0.66, color=OI["blue"], zorder=3)
for i, v in enumerate(bc.n_donors):
    ax.text(i, v + 0.5, str(v), ha="center", fontsize=6.0, color=INK)
ax.set_xticks(x); ax.set_xticklabels([f"{c:.2f}" for c in bc.cps_centre],
                                     fontsize=5.6, rotation=60)
ax.set_xlabel("CPS bin centre"); ax.set_ylabel("donors in bin")
ax.set_ylim(0, 22); ygrid(ax)
ax.axvline(0.5, color=OI["vermillion"], lw=1.0, ls=(0, (4, 2)), zorder=4)
ax.text(0.62, 20.5, "boundary", fontsize=6.0, color=OI["vermillion"])
ax.set_title("equal-width bins are\nnot equal-sized", fontsize=7.2, pad=5)

# ---- (b) detection ----------------------------------------------------------
ax = fig.add_subplot(gs[0, 1:]); panel_label(ax, "b", dx=-0.13)
def verdict(r):
    if r.detection < FLOOR:         return "below floor"
    if r.donors_passing_floor < 60: return "too few donors"
    if r.detection <= AMB_HI:       return "at or below ambient"
    return "observable"
det["v"] = det.apply(verdict, axis=1)
COL = {"below floor": OI["vermillion"], "too few donors": OI["orange"],
       "at or below ambient": OI["purple"], "observable": OI["blue"]}
d = det.sort_values("detection").reset_index(drop=True)
ax.barh(np.arange(len(d)), d.detection, 0.74,
        color=[COL[v] for v in d.v], zorder=3)
ax.axvspan(AMB_LO, AMB_HI, color=MUTED, alpha=0.25, zorder=1)
ax.axvline(FLOOR, color=INK2, lw=0.9, ls=(0, (3, 2)), zorder=4)
ax.set_xscale("log"); ax.set_xlim(2e-4, 1.3)
ax.set_yticks([]); ax.set_ylim(-1, len(d) + 1.6)   # 맨 위에 라벨용 빈 줄 한 칸
ax.set_xlabel("detection rate in astrocytes (log scale)")
ax.text(FLOOR * 0.92, 1.0, "floor 0.02", fontsize=6.1, color=INK2, ha="right")
# 라벨은 맨 위 막대 위에 비워 둔 줄에 둔다. 아래쪽은 범례가 차지하고,
# 띠 안쪽은 막대가 지나가므로 둘 다 겹친다.
ax.text(np.sqrt(AMB_LO * AMB_HI), len(d) + 0.75, "ambient", fontsize=6.1,
        color=INK2, ha="center", va="center")
n = det.v.value_counts()
ax.legend(handles=[Patch(facecolor=COL[k], label=f"{k} ({n.get(k,0)})")
                   for k in ["below floor", "too few donors", "at or below ambient",
                             "observable"]],
          # 막대를 덮지 않도록 우하단 빈 영역으로 옮긴다
          loc="lower right", bbox_to_anchor=(0.995, 0.035), fontsize=6.0,
          handlelength=1.1, framealpha=0.96, borderpad=0.35, labelspacing=0.28)
ax.set_title(f"{int(res.n_panel - res.n_observable)} of {int(res.n_panel)} panel markers "
             f"cannot be adjudicated", fontsize=7.2, pad=5)
xgrid(ax)

# ---- (c) restriction --------------------------------------------------------
ax = fig.add_subplot(gs[1, 0]); panel_label(ax, "c", dx=-0.30)
vals = [(f"all {int(res.n_panel)}\nmarkers", res.bp_all, res.p_ward_all),
        (f"{int(res.n_observable)} observable\nonly", res.bp_observable, res.p_ward_obs)]
for i, (lab, bp, pv) in enumerate(vals):
    ax.scatter([bp], [i], s=70, marker="D", color=OI["blue"], zorder=4)
    ax.text(bp, i + 0.30, f"{bp:.3f}", ha="center", fontsize=6.6,
            color=INK, fontweight="bold")
    ax.text(0.93, i - 0.30, f"p = {pv:.2f}", ha="right", fontsize=6.2, color=INK2)
ax.set_yticks([0, 1]); ax.set_yticklabels([v[0] for v in vals], fontsize=6.4,
                                          linespacing=1.25)
ax.set_ylim(-0.6, 1.6); ax.set_xlim(0.12, 0.95)
ax.set_xlabel("Ward first boundary (CPS)"); xgrid(ax)
ax.set_title("dropping the 23 changes\nnothing", fontsize=7.2, pad=5)

# ---- (d) leave-one-donor ----------------------------------------------------
ax = fig.add_subplot(gs[1, 1]); panel_label(ax, "d", dx=-0.22)
rng = np.random.default_rng(2)
yy = rng.uniform(0.2, 0.8, len(L))
inside = (L.bp - BP).abs() <= TOL
ax.scatter(L.bp[inside], yy[inside.to_numpy()], s=13, color=OI["blue"],
           alpha=0.75, zorder=3, linewidth=0)
ax.scatter(L.bp[~inside], yy[(~inside).to_numpy()], s=28, color=OI["vermillion"],
           zorder=4, linewidth=0)
for _, r in L[~inside].iterrows():
    ax.annotate(r.donor, (r.bp, yy[L.index[L.donor == r.donor][0]]),
                xytext=(0, 9), textcoords="offset points", fontsize=5.8,
                color=OI["vermillion"], ha="center")
ax.axvspan(BP - TOL, BP + TOL, color=OI["blue"], alpha=0.10, zorder=1)
ax.axvline(BP, color=INK2, lw=0.9, ls=(0, (4, 2)), zorder=2)
ax.set_yticks([]); ax.set_ylim(0, 1.05); ax.set_xlim(0.15, 0.40)
ax.set_xlabel("boundary with one donor removed (CPS)")
ax.set_title(f"{int((~inside).sum())} of {len(L)} donors move it", fontsize=7.2, pad=5)
xgrid(ax)
for sp in ("left",): ax.spines[sp].set_visible(False)

# ---- (e) bootstrap ----------------------------------------------------------
ax = fig.add_subplot(gs[1, 2]); panel_label(ax, "e", dx=-0.30)
bv = B.bp.dropna().to_numpy()
ax.hist(bv, bins=np.arange(0.15, 0.95, 0.03), color=OI["blue"], zorder=3)
lo, hi = np.percentile(bv, [2.5, 97.5])
ax.axvline(BP, color=OI["vermillion"], lw=1.2, zorder=5)
ax.set_xlabel("donor bootstrap (CPS)"); ax.set_ylabel("resamples")
ax.set_xlim(0.15, 0.80); ygrid(ax)
ax.text(0.97, 0.95, f"median {np.median(bv):.3f}\n95% CI {lo:.3f}–{hi:.3f}\n"
                    f"{np.mean(np.abs(bv-BP)<=TOL)*100:.0f}% within 0.05",
        transform=ax.transAxes, fontsize=6.0, ha="right", va="top",
        color=INK2, linespacing=1.4)
ax.set_title("but the interval is wide", fontsize=7.2, pad=5)

fig.savefig(FIG + "Figure3.png")
fig.savefig(FIG + "Figure3.pdf")
print(f"saved Figure3 | {int((~inside).sum())}/{len(L)} donors | CI {lo:.3f}-{hi:.3f} | "
      f"{int(res.n_panel-res.n_observable)} of {int(res.n_panel)} unadjudicable")
