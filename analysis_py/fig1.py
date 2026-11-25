"""Figure 1 - the audit, its exits, and where this paper's own claims stop."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import numpy as np
from style import setup, OI, INK, INK2, MUTED
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

setup()
fig = plt.figure(figsize=(7.3, 5.3))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, text, fc="#ffffff", ec=INK2, fs=7.2, lw=0.8, weight=None, tc=INK, r=1.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, zorder=4, linespacing=1.35, fontweight=weight)

def arrow(x1, y1, x2, y2, color=INK2, lw=0.9):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=8, linewidth=lw, color=color, zorder=2))

ax.text(2, 97.2, "a", fontsize=10, fontweight="bold", color=INK)
ax.text(6.5, 97.2, "A reported boundary enters the audit", fontsize=8.2, color=INK, va="center")
box(3, 85.2, 17, 7.6, "reported\nboundary", fc="#f3f3f5", fs=7.4, weight="bold")

# The order is logical: a position that matched noise also returns is not worth
# asking further questions about, so specificity comes before redundancy.
gates = [
    ("G1  observability", "which markers clear the detection\nfloor and the ambient band?", OI["purple"]),
    ("G2  specificity", "is it absent from matched random\npanels, scored on the same statistic?", OI["green"]),
    ("G3  independence", "do the agreeing components carry\nindependent evidence?", OI["orange"]),
    ("G4  donor stability", "does it survive leave-one-donor-out\nand a donor bootstrap?", OI["vermillion"]),
    ("G5  region / cohort", "is it recovered where composition\nallows the same question?", OI["blue"]),
]
EXIT = [("not observable", OI["vermillion"]),
        ("not specific\n— a property of the grid", OI["green"]),
        ("unevaluable\n— agreement is redundancy", OI["orange"]),
        ("unevaluable\n— carried by one donor", OI["orange"]),
        ("dataset-specific", OI["purple"])]

GY, x0, w, h, dy = 85.2, 25, 39, 7.8, 9.9
ys = [GY - i * dy for i in range(len(gates))]
for i, ((nm, sub, col), y) in enumerate(zip(gates, ys)):
    ax.add_patch(Rectangle((x0, y), 1.5, h, color=col, zorder=4))
    box(x0, y, w, h, "", fc="#ffffff", ec=INK2, lw=0.7)
    ax.text(x0 + 3.4, y + h * 0.68, nm, fontsize=7.2, fontweight="bold", color=INK, va="center", zorder=5)
    ax.text(x0 + 3.4, y + h * 0.25, sub, fontsize=6.0, color=INK2, va="center", zorder=5, linespacing=1.3)
    if i == 0: arrow(20, 89.0, x0 - 0.6, y + h / 2)
    else:      arrow(x0 + w / 2, ys[i - 1], x0 + w / 2, y + h + 0.4)
    arrow(x0 + w + 0.4, y + h / 2, 70.5, y + h / 2, color=EXIT[i][1], lw=0.9)
    ax.text(71.5, y + h / 2, EXIT[i][0], fontsize=6.4, color=EXIT[i][1],
            va="center", ha="left", fontweight="bold", linespacing=1.25)

arrow(x0 + w / 2, ys[-1], x0 + w / 2, ys[-1] - 4.0)
box(x0 + 4.5, ys[-1] - 11.4, 30, 7.0, "generalizable", fc="#eef6f2", ec=OI["green"],
    fs=7.6, weight="bold", tc=INK, lw=1.0)
ax.text(49, ys[-1] - 14.4, "the earliest failing gate sets the verdict; every gate outcome is reported",
        fontsize=6.2, ha="center", color=INK2, style="italic")

# ---------------- the three cases --------------------------------------------
ax.text(2, 25.0, "b", fontsize=10, fontweight="bold", color=INK)
ax.text(6.5, 25.0, "The three claims audited here, and the gate at which each stops",
        fontsize=8.2, color=INK, va="center")

cases = [
    # Verdict wording follows Table 1 and Results: p = 0.079 does not certify the
    # window, but it is not the same exit as the panel's p = 0.753 — the window
    # departs from the position the grid forces, and the panel does not. Calling
    # both "not specific" would erase the contrast the paper rests on.
    ("CPS 0.4–0.6\n(atlas authors)", "G1 pass\nG2  p = 0.079", "not certified",
     OI["green"], "#eef6f2"),
    ("CPS 0.207\n(this panel)", "G1 pass on restriction\nG2  p = 0.753", "not specific",
     OI["orange"], "#fdf4e8"),
    ("HMOX1\n(largest effect)", "G1 fail\n25 of 84 donors", "not observable",
     OI["vermillion"], "#fdeee8"),
]
bw, bx0, gap = 28.5, 5.5, 3.2
for i, (who, where, verdict, col, fc) in enumerate(cases):
    x = bx0 + i * (bw + gap)
    box(x, 5.2, bw, 15.3, "", fc=fc, ec=col, lw=1.0)
    ax.text(x + bw / 2, 17.3, who, fontsize=7.4, ha="center", va="center", color=INK,
            fontweight="bold", linespacing=1.3, zorder=5)
    ax.text(x + bw / 2, 12.2, where, fontsize=6.2, ha="center", va="center", color=INK2,
            zorder=5, linespacing=1.35)
    ax.add_patch(Rectangle((x + 3.0, 6.5), bw - 6, 3.6, facecolor=col, alpha=0.18,
                           edgecolor="none", zorder=4))
    ax.text(x + bw / 2, 8.3, verdict, fontsize=6.9, ha="center", va="center", color=col,
            fontweight="bold", zorder=5, linespacing=1.2)

ax.text(50, 2.0, "No claim examined here reaches “generalizable”, and only the atlas authors’ window "
        "departs from the position the grid forces (Figure 5).\nThat the audit certifies planted "
        "boundaries in 16 of 16 synthetic panels (Figure 2) is what separates this result from a "
        "procedure that rejects\neverything. A null at G1 is not evidence of no change — only that "
        "the assay cannot see it.",
        fontsize=6.2, ha="center", va="center", color=INK2, style="italic", linespacing=1.45)

fig.savefig(FIG + "Figure1.png")
fig.savefig(FIG + "Figure1.pdf")
print("saved Figure1")
