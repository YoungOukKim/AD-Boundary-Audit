"""Figure 2 — ground-truth benchmark, and the artefact the published gates miss."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import numpy as np, pandas as pd
from style import setup, OI, INK, INK2, MUTED, GRID, panel_label, ygrid, xgrid
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

setup()
D = DATA
# NOTE: the regime literal "null" is a pandas NA token by default, so the whole
# null-panel arm reads back as NaN unless default NA handling is switched off.
d = pd.read_csv(D + "benchmark.csv", keep_default_na=False, na_values=[""])
for c in ["g_support", "g_spec", "g_lomo", "g_lodo", "certified", "certified_old"]:
    d[c] = d[c].astype(str).str.strip().str.lower().isin(["true", "1"])
for c in ["primary", "support", "pspec", "lomo", "lodo"]:
    d[c] = pd.to_numeric(d[c], errors="coerce")
meta = pd.read_csv(D + "benchmark_meta.csv")
NREP = int(meta.nrep[0])

NAME = {"real": "real boundary", "global": "global offset", "null": "null panel",
        "ramp": "monotone ramp", "decoy": "single-marker decoy",
        "donor": "single-donor artefact"}
ORD = ["real", "global", "null", "ramp", "decoy", "donor"]

fig = plt.figure(figsize=(7.3, 3.35))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.1, 0.85], wspace=0.52)

# ---- (a) certification rate, published gates vs published + donor gate ------
ax = fig.add_subplot(gs[0, 0]); panel_label(ax, "a", dx=-0.34)
old = [d[d.regime == k].certified_old.sum() for k in ORD]
new = [d[d.regime == k].certified.sum() for k in ORD]
y = np.arange(len(ORD))[::-1]; h = 0.36
ax.barh(y + h/2, old, h, color=MUTED, zorder=3, label="published gates (3)")
ax.barh(y - h/2, new, h, color=OI["blue"], zorder=3, label="+ leave-one-donor-out")
ax.barh(y[-1] + h/2, old[-1], h, color=OI["vermillion"], zorder=4)
ax.set_yticks(y); ax.set_yticklabels([NAME[k] for k in ORD], fontsize=6.8)
ax.set_xlim(0, NREP + 3.4); ax.set_xlabel(f"panels certified (of {NREP})")
xgrid(ax)
ax.text(old[-1] + 0.5, y[-1] + h/2, "16/16", fontsize=6.6, va="center",
        color=OI["vermillion"], fontweight="bold")
ax.text(0.4, y[-1] - h/2, "0/16", fontsize=6.6, va="center", ha="left",
        color=OI["blue"], fontweight="bold")
ax.legend(loc="upper center", bbox_to_anchor=(0.50, -0.26), fontsize=6.4,
          ncol=2, columnspacing=1.1)
ax.set_title("what the audit certifies", fontsize=7.3, pad=5)

# ---- (b) which gate does the work -------------------------------------------
ax = fig.add_subplot(gs[0, 1]); panel_label(ax, "b", dx=-0.26)
G = ["g_support", "g_spec", "g_lomo", "g_lodo"]
GL = ["support", "speci-\nficity", "marker\nLOO", "donor\nLOO"]
M = np.array([[d[d.regime == k][g].mean() for g in G] for k in ORD])
im = ax.imshow(M, cmap="RdYlBu", vmin=0, vmax=1, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, f"{M[i,j]*100:.0f}", ha="center", va="center", fontsize=6.4,
                color=INK if 0.25 < M[i, j] < 0.8 else "white", fontweight="bold")
ax.set_xticks(range(4)); ax.set_xticklabels(GL, fontsize=6.3, linespacing=1.2)
ax.set_yticks(range(len(ORD))); ax.set_yticklabels([NAME[k] for k in ORD], fontsize=6.6)
ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(ORD), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.4)
ax.tick_params(which="minor", length=0)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_title("% of panels passing each gate", fontsize=7.3, pad=5)

# ---- (c) the quantity the new gate measures ---------------------------------
ax = fig.add_subplot(gs[0, 2]); panel_label(ax, "c", dx=-0.42)
for i, k in enumerate(["real", "donor"]):
    v = d[d.regime == k].lodo.to_numpy()
    xx = np.full_like(v, i, dtype=float) + (np.random.default_rng(i).random(len(v)) - .5) * 0.22
    ax.scatter(xx, v, s=15, color=(OI["blue"] if k == "real" else OI["vermillion"]),
               zorder=3, linewidth=0, alpha=0.85)
ax.axhline(0.05, color=INK2, lw=0.8, ls=(0, (3, 2)), zorder=2)
ax.text(-0.44, 0.072, "gate 0.05", fontsize=6.3, color=INK2, ha="left",
        va="bottom")
ax.set_xticks([0, 1]); ax.set_xticklabels(["real\nboundary", "single-donor\nartefact"],
                                          fontsize=6.5, linespacing=1.25)
ax.set_xlim(-0.5, 1.5); ax.set_ylim(-0.03, 0.60)
ax.set_ylabel("largest shift when one\ndonor is removed (CPS)", fontsize=7)
ygrid(ax)
ax.set_title("what leave-one-donor-out sees", fontsize=7.3, pad=5)

fig.savefig(FIG + "Figure2.png")
fig.savefig(FIG + "Figure2.pdf")
sens = d[d.regime == "real"].certified.mean()
spec_new = 1 - d[d.regime != "real"].certified.mean()
spec_old = 1 - d[d.regime != "real"].certified_old.mean()
print(f"saved Figure2 | sensitivity {sens*100:.0f}% | specificity {spec_old*100:.0f}% -> {spec_new*100:.0f}%")
