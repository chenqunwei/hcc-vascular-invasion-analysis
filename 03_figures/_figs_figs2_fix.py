# -*- coding: utf-8 -*-
"""FigS2 轨迹图重绘（修正 pattern 0/1/2 语义映射，从已保存 z 矩阵读取，无需重读 h5ad）。"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 9,
    "axes.linewidth": 0.6, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 6,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})
SITES = ["Normal", "Tumor", "PVTT", "Lymph"]
UP_PAT = {0: ("tumor-onset", "#5A3E9B"), 1: ("PVTT-onset", "#2A8C99"), 2: ("metastasis-onset", "#9ACC28")}
DN_PAT = {0: ("PVTT deep loss with lymph-node re-expression", "#5A3E9B"),
          1: ("progressive loss", "#2A8C99"), 2: ("PVTT-preserved", "#9ACC28")}

def tag(ax, letter, x=-0.15):
    ax.text(x, 1.08, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

fig = plt.figure(figsize=(11.5, 3.6))
gs = fig.add_gridspec(2, 4, height_ratios=[1.5, 1], hspace=0.6, wspace=0.6)
cmap = LinearSegmentedColormap.from_list("z", ["#FFFFFF", "#F9C0AA", "#D85A30", "#7A1E0D"])

for row, (name, PAT) in enumerate([("up", UP_PAT), ("dn", DN_PAT)]):
    z = pd.read_csv(f"results/trajectory_zscore_{name}.csv", index_col=0).values
    pat = pd.read_csv(f"results/trajectory_{name}_patterns.csv", index_col=0)["pattern"].values
    o = np.lexsort((z[:, 1], pat))
    zm, pm = z[o], pat[o]

    ax = fig.add_subplot(gs[row, 0:2])
    im = ax.imshow(zm, aspect="auto", cmap=cmap, vmin=-2, vmax=2, interpolation="nearest")
    ax.set_xticks(range(4)); ax.set_xticklabels(SITES, fontsize=7)
    ax.set_yticks([])
    ax.set_title(f"{name.capitalize()}-regulated program ({len(zm)} genes)", fontweight="bold", fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02, label="z-score")
    tag(ax, "A" if row == 0 else "C")

    ax = fig.add_subplot(gs[row, 2:4])
    for c in [0, 1, 2]:
        lab, color = PAT[c]
        ax.plot(range(4), zm[pm == c].mean(axis=0), "o-", color=color, lw=1.6,
                label=f"{lab} (n={int((pm == c).sum())})")
    ax.axhline(0, color="gray", ls="--", lw=0.6)
    ax.set_xticks(range(4)); ax.set_xticklabels(SITES)
    ax.set_ylabel("Mean z-score")
    ax.set_title(f"{name.capitalize()}-regulated mean trajectories", fontweight="bold", fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc="best")
    tag(ax, "B" if row == 0 else "D")

fig.savefig("figures/Fig7_trajectory_patterns.png")
print("FigS2 重绘完成（pattern 0/1/2 语义修正）")