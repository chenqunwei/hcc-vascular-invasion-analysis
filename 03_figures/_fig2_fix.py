# -*- coding: utf-8 -*-
"""Fig2 独立重画：图内标签 "Δ area" → "ΔCDF area"（与图注/正文统一）。
从 _figs_core.py 的 Fig2 部分复制，单独输出到 results/ + figures/，避免覆盖 Fig1。
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 9,
    "axes.linewidth": 0.6, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})
SUB_COLORS = {"Low VI": "#1D9E75", "Intermediate VI": "#EF9F27", "High VI": "#D85A30"}
SUB_ORDER = ["Low VI", "Intermediate VI", "High VI"]

def tag(ax, letter, x=-0.12):
    ax.text(x, 1.08, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

clusters = pd.read_csv("results/tcga_lihc_clusters.csv")
vi_order = {1: "Low VI", 3: "Intermediate VI", 2: "High VI"}
clusters["subtype"] = clusters["cluster"].map(vi_order)
clusters = clusters[clusters["mvi"].isin([0, 1])]

cons_mats = np.load("results/consensus_matrices.npy", allow_pickle=True)

fig = plt.figure(figsize=(10, 3.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1], wspace=0.4)

# A 共识热图 k=3
ax = fig.add_subplot(gs[0, 0])
cons = cons_mats[1]
dist = 1 - cons
np.fill_diagonal(dist, 0)
Z = linkage(squareform(dist), method="average")
order = leaves_list(Z)
cs = cons[np.ix_(order, order)]
im = ax.imshow(cs, cmap="Reds", vmin=0, vmax=1, aspect="auto")
ax.set_title("Consensus matrix (k=3)", fontweight="bold")
ax.set_xticks([]); ax.set_yticks([])
fig.colorbar(im, ax=ax, shrink=0.7, pad=0.04, label="Consensus")
tag(ax, "A")

# B 共识 CDF 曲线（ΔCDF area）
ax = fig.add_subplot(gs[0, 1])
cons_all = np.load("results/consensus_matrices.npy", allow_pickle=True)
vals = np.linspace(0, 1, 101)
areas = {}
for i, k in enumerate(range(2, 6)):
    triu = cons_all[i][np.triu_indices_from(cons_all[i], k=1)]
    cdf = np.array([np.mean(triu <= v) for v in vals])
    areas[k] = np.trapezoid(cdf, vals)
    ax.plot(vals, cdf, lw=1.4, label=f"k = {k}")
delta3 = areas[3] - areas[2]
ax.set_xlabel("Consensus index")
ax.set_ylabel("CDF")
ax.set_title("Consensus CDF", fontweight="bold")
ax.legend(frameon=False, fontsize=7, loc="lower right",
          title=f"ΔCDF area (k=3 vs 2) = {delta3:.3f}", title_fontsize=7)
tag(ax, "B", x=-0.15)

# C 三型 VI score 箱线
ax = fig.add_subplot(gs[0, 2])
data = [clusters[clusters["subtype"] == s]["ssgsea_vi"].values for s in SUB_ORDER]
bp = ax.boxplot(data, tick_labels=SUB_ORDER, patch_artist=True, widths=0.6)
for patch, s in zip(bp["boxes"], SUB_ORDER):
    patch.set_facecolor(SUB_COLORS[s]); patch.set_alpha(0.6)
for w in bp["whiskers"]: w.set_color("gray")
for c in bp["caps"]: c.set_color("gray")
for m in bp["medians"]: m.set_color("black")
ax.set_ylabel("VI score")
ax.set_title("VI score by subtype", fontweight="bold")
tag(ax, "C", x=-0.15)

fig.savefig("results/Fig2_molecular_subtyping.png")
fig.savefig("figures/Fig2_molecular_subtyping.png")
print("Fig2 重画完成（Δ area → ΔCDF area；results/ + figures/）")
