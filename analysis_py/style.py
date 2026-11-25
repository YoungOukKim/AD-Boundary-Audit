"""Shared figure style for the P3 revision. Journal-figure conventions:
static, print-safe, colourblind-safe (Okabe-Ito), identity never carried by
colour alone (direct labels + distinct markers)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import DATA, INPUT, FIG
import matplotlib as mpl
import matplotlib.pyplot as plt

# Okabe-Ito, ordered so the deutan-close pair (#CC79A7 / #009E73) is never adjacent
OI = dict(blue="#0072B2", vermillion="#D55E00", green="#009E73", orange="#E69F00",
          purple="#CC79A7", sky="#56B4E9", yellow="#F0E442", black="#000000")
SEQ = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
MARKERS = ["o", "s", "^", "D", "v", "P"]

INK = "#1a1a1a"; INK2 = "#55555f"; MUTED = "#8a8a93"; GRID = "#e3e3e6"
SURF = "#ffffff"

def setup():
    mpl.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 400,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": 7.5, "axes.labelsize": 8, "axes.titlesize": 8.5,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.edgecolor": INK2, "axes.linewidth": 0.6,
        "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK2, "ytick.color": INK2,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.5,
        "legend.frameon": False, "legend.handlelength": 1.4,
        "legend.handletextpad": 0.5, "legend.labelspacing": 0.35,
        "lines.linewidth": 1.3, "lines.markersize": 4,
        "figure.facecolor": SURF, "axes.facecolor": SURF,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

def panel_label(ax, s, dx=-0.16, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="top", ha="left", color=INK)

def ygrid(ax):
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.5)

def xgrid(ax):
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color=GRID, linewidth=0.5)
