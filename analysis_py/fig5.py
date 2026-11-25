"""Figure 5 - the positive control: on the same grid, staging data escapes the
position the grid forces, and the panel does not."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import numpy as np, pandas as pd
from style import setup, OI, SEQ, MARKERS, INK, INK2, MUTED, GRID, panel_label, ygrid, xgrid
import matplotlib.pyplot as plt

setup()
D = DATA
pa = pd.read_csv(D + "pathology_axes.csv")
pe = pd.read_csv(D + "pathology_engine.csv")
pn = pd.read_csv(D + "pathology_null.csv").bp.to_numpy(float)
nulls = pd.read_csv(D + "nulls.csv")
bc = pd.read_csv(D + "bin_composition.csv")
ax_c = bc.cps_centre.to_numpy()
WIN, OURS, TOL = (0.40, 0.60), 0.20683, 0.05
PBP = float(pe[pe.component == "Hier_k3_bp1"].breakpoint.iloc[0])

fig = plt.figure(figsize=(7.3, 2.95))
gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.18, 0.92], wspace=0.44)

# ---- (a) staging along CPS --------------------------------------------------
ax = fig.add_subplot(gs[0, 0]); panel_label(ax, "a", dx=-0.26)
for i, (k, lab) in enumerate([("braak", "Braak"), ("thal", "Thal"),
                              ("cerad", "CERAD"), ("adnc", "ADNC")]):
    r = pa[pa.axis == k].iloc[0]
    y = np.array([r[f"bin{j+1}"] for j in range(9)], float)
    y = (y - y.min()) / (y.max() - y.min())
    ax.plot(ax_c, y, marker=MARKERS[i], color=SEQ[i], ms=3.3, lw=1.2, label=lab,
            zorder=3, markeredgewidth=0)
ax.axvspan(*WIN, color=OI["green"], alpha=0.12, zorder=0)
ax.axvline(OURS, color=OI["vermillion"], lw=1.0, ls=(0, (4, 2)), zorder=2)
ax.text(np.mean(WIN), 1.07, "reported 0.4–0.6", fontsize=6.2, color=OI["green"],
        ha="center", va="bottom")
ax.text(OURS - 0.012, 1.02, "panel boundary", fontsize=6.0, color=OI["vermillion"],
        rotation=90, va="top", ha="right")
ygrid(ax); ax.set_xlim(0.1, 0.96); ax.set_ylim(-0.06, 1.20)
ax.set_xlabel("CPS"); ax.set_ylabel("staging, scaled 0–1")
ax.legend(loc="lower right", fontsize=6.1, ncol=2, columnspacing=0.8,
          handlelength=1.2, borderpad=0.3, labelspacing=0.25)
ax.set_title("neuropathological staging", fontsize=7.3, pad=6)

# ---- (b) the same engine, the same grid, two inputs -------------------------
ax = fig.add_subplot(gs[0, 1]); panel_label(ax, "b", dx=-0.20)
rng = np.random.default_rng(3)
v = pe.breakpoint.to_numpy()
w = nulls.ward_null.dropna().to_numpy()
ax.axvspan(*WIN, color=OI["green"], alpha=0.13, zorder=0)
ax.axvspan(OURS - TOL, OURS + TOL, color=OI["vermillion"], alpha=0.11, zorder=0)
# null density strip
h, e = np.histogram(w, bins=np.arange(0.10, 1.0, 0.035), density=True)
ax.bar((e[:-1] + e[1:]) / 2, h / h.max() * 0.30, width=0.033, bottom=0.03,
       color=MUTED, alpha=0.75, zorder=2)
ax.text(0.935, 0.20, "Ward null\n(random panels)", fontsize=5.9, color=INK2,
        ha="right", va="center", linespacing=1.25)
ax.scatter(v, rng.uniform(0.52, 0.88, len(v)), s=24, facecolor="none",
           edgecolor=INK2, linewidth=0.8, zorder=4)
ax.scatter([OURS], [0.44], s=52, marker="D", color=OI["vermillion"], zorder=5)
ax.text(OURS, 0.38, "panel", fontsize=6.0, color=OI["vermillion"], ha="center",
        va="top")
nin = int(((v >= WIN[0]) & (v <= WIN[1])).sum()); nour = int((np.abs(v - OURS) <= TOL).sum())
ax.text(np.mean(WIN), 0.96, f"{nin} of 16", fontsize=7.6, color=OI["green"],
        ha="center", fontweight="bold")
ax.text(OURS, 0.96, f"{nour} of 16", fontsize=7.6, color=OI["vermillion"],
        ha="center", fontweight="bold")
ax.set_yticks([]); ax.set_ylim(0, 1.06); ax.set_xlim(0.10, 0.96)
ax.set_xlabel("breakpoint on the four staging axes (CPS)")
xgrid(ax)
for sp in ("left",): ax.spines[sp].set_visible(False)
ax.set_title(f"staging escapes the forced position", fontsize=7.3, pad=6)

# ---- (c) is CPS privileged to tau? ------------------------------------------
ax = fig.add_subplot(gs[0, 2]); panel_label(ax, "c", dx=-0.36)
ordr = ["adnc", "thal", "braak", "cerad", "cogstat", "late", "lewy"]
nm = {"adnc": "ADNC", "thal": "Thal (Aβ)", "braak": "Braak (tau)", "cerad": "CERAD",
      "cogstat": "dementia", "late": "LATE", "lewy": "Lewy"}
r = [float(pa[pa.axis == k].rho_cps.iloc[0]) for k in ordr]
p = [float(pa[pa.axis == k].p_cps.iloc[0]) for k in ordr]
col = [OI["blue"]] * 5 + [MUTED] * 2
yy = np.arange(len(ordr))[::-1]
ax.barh(yy, r, color=col, height=0.64, zorder=3)
ax.barh(yy[2], r[2], color=OI["vermillion"], height=0.64, zorder=4)
for i, (rv, pv) in enumerate(zip(r, p)):
    ax.text(rv + 0.03, yy[i], f"{rv:+.2f}" + ("" if pv < .01 else " n.s."),
            va="center", fontsize=6.0, color=INK2)
ax.set_yticks(yy); ax.set_yticklabels([nm[k] for k in ordr], fontsize=6.4)
ax.set_xlim(0, 1.05); xgrid(ax)
ax.set_xlabel("Spearman $\\rho$ with CPS")
ax.set_title("CPS is not tau-privileged", fontsize=7.3, pad=6)

fig.savefig(FIG + "Figure5.png")
fig.savefig(FIG + "Figure5.pdf")
print(f"saved Figure5 | pathology Ward {PBP:.3f} | in window {nin}/16 | at panel boundary {nour}/16 "
      f"| null median {np.median(w):.3f}")
