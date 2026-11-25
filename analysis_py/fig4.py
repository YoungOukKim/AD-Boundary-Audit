"""Figure 4 - why the boundary fails: the null it was tested against, the grid
that produces it, and the components that are not independent."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import numpy as np, pandas as pd
from style import setup, OI, INK, INK2, MUTED, GRID, panel_label, ygrid, xgrid
import matplotlib.pyplot as plt

setup()
D = DATA
nulls = pd.read_csv(D + "nulls.csv")
grid = pd.read_csv(D + "grid_bias.csv")
ind = pd.read_csv(D + "independence.csv")
frag = pd.read_csv(D + "reconstruction_fragility.csv")
res = pd.read_csv(D + "restriction.csv").iloc[0]
BP = float(res.bp_all)
TOL = 0.05

fig = plt.figure(figsize=(7.3, 5.9))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.40,
                      height_ratios=[1.0, 0.95])

# ---- (a) the two nulls ------------------------------------------------------
ax = fig.add_subplot(gs[0, 0]); panel_label(ax, "a", dx=-0.21)
bins = np.arange(0.10, 1.0, 0.035)
w = nulls.ward_null.dropna().to_numpy(); p = nulls.pc1_null.dropna().to_numpy()
ax.hist(p, bins=bins, color=MUTED, alpha=0.95, zorder=2,
        label="PC1 breakpoint  (the null as published)")
ax.hist(w, bins=bins, color=OI["blue"], alpha=0.80, zorder=3,
        label="Ward boundary  (like for like)")
ax.axvline(BP, color=OI["vermillion"], lw=1.4, zorder=5)
ax.text(BP + 0.012, ax.get_ylim()[1] * 0.55, "observed\n0.207", fontsize=6.4,
        color=OI["vermillion"], va="top", ha="left", fontweight="bold", linespacing=1.25)
ax.set_xlabel("breakpoint in 1,000 expression-matched random panels (CPS)")
ax.set_ylabel("panels"); ygrid(ax); ax.set_xlim(0.10, 0.95)
ax.legend(loc="upper right", fontsize=5.9, handlelength=1.1, framealpha=0.92,
          borderpad=0.35, labelspacing=0.3)
ax.set_title("the same boundary, two nulls", fontsize=7.4, pad=5)

# ---- (b) what each comparison returns ---------------------------------------
ax = fig.add_subplot(gs[0, 1]); panel_label(ax, "b", dx=-0.24)
f = lambda n, x: float(np.mean(np.abs(n - x) <= TOL))
rows = [("Ward 0.207\nvs PC1 null", f(p, BP), OI["vermillion"], "as published"),
        ("Ward 0.207\nvs Ward null", f(w, BP), OI["blue"], "like for like"),
        ("PC1 0.543\nvs PC1 null", f(p, 0.542979), MUTED, "like for like")]
y = np.arange(len(rows))[::-1]
for i, (lab, v, c, note) in enumerate(rows):
    ax.barh(y[i], v, 0.52, color=c, zorder=3)
    inside = v > 0.45
    # 0/1000 은 p = 0.000 이 아니라 p < 0.001 이다. p값 오용을 다루는 논문이므로
    # 자기 그림에서 0.000 으로 적지 않는다.
    ptxt = "p < 0.001" if v < 0.001 else f"p = {v:.3f}"
    ax.text(v - 0.02 if inside else v + 0.02, y[i], ptxt, va="center",
            ha="right" if inside else "left", fontsize=6.8,
            color="white" if inside else INK, fontweight="bold", zorder=5)
    ax.text(1.0, y[i] - 0.33, note, va="center", ha="right", fontsize=6.0,
            color=INK2, style="italic")
ax.axvline(0.05, color=INK2, lw=0.8, ls=(0, (3, 2)), zorder=4)
ax.text(0.058, y[0] + 0.38, "0.05", fontsize=6.2, color=INK2)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=6.6, linespacing=1.25)
ax.set_xlim(0, 1.02); ax.set_xlabel("fraction of random panels within 0.05 CPS")
xgrid(ax)
ax.set_title("the p-value is a property of the null", fontsize=7.4, pad=5)

# ---- (c) the grid produces the boundary -------------------------------------
ax = fig.add_subplot(gs[1, 0]); panel_label(ax, "c", dx=-0.21)
g = grid[grid.tag.str.startswith("9 bins")].sort_values("first_bin_n")
ax.plot(g.first_bin_n, g.ward_first_mid * 100, "-o", color=OI["blue"], ms=4,
        lw=1.4, zorder=4, label="Ward", markeredgewidth=0)
ax.plot(g.first_bin_n, g.km_first_mid * 100, "-s", color=OI["orange"], ms=3.4,
        lw=1.2, zorder=3, label="k-means", markeredgewidth=0)
ax.axhline(grid.random_label_baseline.iloc[0] * 100, color=INK2, lw=0.8,
           ls=(0, (3, 2)), zorder=2)
ax.text(5.6, grid.random_label_baseline.iloc[0] * 100 - 6.5,
        "any random 3-labelling", fontsize=6.1, color=INK2, ha="center")
ax.axhline(12.5, color=MUTED, lw=0.8, ls=(0, (1, 2)), zorder=2)
ax.text(5.6, 15.5, "chance (1 of 8 positions)", fontsize=6.1, color=MUTED, ha="center")
real = grid[grid.tag.str.startswith("real")].iloc[0]
ax.scatter([2.0], [real.ward_first_mid * 100], s=46, marker="*",
           color=OI["vermillion"], zorder=6)
ax.text(2.35, real.ward_first_mid * 100 - 1.5, "SEA-AD grid\n(100%)", fontsize=6.2,
        color=OI["vermillion"], va="top", linespacing=1.2)
ax.invert_xaxis()
ax.set_xlabel("donors in the first bin (nine-bin grid, pure noise)")
ax.set_ylabel("first cut lands at the\nfirst midpoint (%)", fontsize=7)
ax.set_ylim(0, 108); ygrid(ax)
ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.30), fontsize=6.4, handlelength=1.4)
ax.set_title("on noise, the grid decides where the cut falls", fontsize=7.4, pad=5)

# ---- (d) components are not independent -------------------------------------
ax = fig.add_subplot(gs[1, 1]); panel_label(ax, "d", dx=-0.24)
top = ind.head(4).iloc[::-1]
lab = [f"{a.replace('_k3','').replace('Hier','Ward').replace('KMeans','k-means')}\n"
       f"= {b.replace('_k3','').replace('KMeans','k-means').replace('Variance_jump','var. jump')}"
       for a, b in zip(top.a, top.b)]
yy = np.arange(len(top))
ax.barh(yy, top.identical * 100, 0.56, color=OI["purple"], zorder=3)
for i, v in enumerate(top.identical * 100):
    ax.text(v + 1.6, yy[i], f"{v:.0f}%", va="center", fontsize=6.6, color=INK,
            fontweight="bold")
ax.set_yticks(yy); ax.set_yticklabels(lab, fontsize=6.1, linespacing=1.2)
ax.set_xlim(0, 100); ax.set_xlabel("random panels returning the identical value (%)")
xgrid(ax)
moved = frag[(frag.published_matrix - frag.reconstructed).abs() > TOL]
ax.text(0.5, -0.38,
        f"and {len(moved)} of 16 components move more than 0.05 CPS when the\n"
        f"same matrix is rebuilt with a 1% change in scaling — including\n"
        f"variance jump (0.207 → 0.886), one of the four supporting the boundary",
        transform=ax.transAxes, fontsize=6.1, color=INK2, ha="center", va="top",
        style="italic", linespacing=1.4)
ax.set_title("agreement is redundancy, not support", fontsize=7.4, pad=5)

fig.savefig(FIG + "Figure4.png")
fig.savefig(FIG + "Figure4.pdf")
print(f"saved Figure4 | ward-null p={f(w,BP):.3f} pc1-null p={f(p,BP):.3f} | "
      f"grid real {real.ward_first_mid*100:.0f}%")
