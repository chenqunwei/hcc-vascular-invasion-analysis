# -*- coding: utf-8 -*-
"""Fig6 单细胞定位重排：panel 顺序与正文论述逻辑一致。
A=VI by site（先讲梯度）→ B=UMAP VI score → C=UMAP cluster2 高亮 → D=top clusters。
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 9,
    "axes.linewidth": 0.6, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})
C_HI = "#D85A30"
SITE_COLORS = {"Tumor": "#85B7EB", "PVTT": "#D85A30", "Lymph": "#7F77DD"}

def fmt_p(p):
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}"

def tag(ax, letter, x=-0.12):
    ax.text(x, 1.08, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

print("读取 h5ad UMAP 坐标（backed 模式）...")
adata = sc.read_h5ad("results/gse149614_processed_qc.h5ad", backed="r")
umap_all = adata.obsm["X_umap"]
malig_mask = ((adata.obs["celltype"] == "Hepatocyte") & (adata.obs["site"] != "Normal")).values
# Unify to A (full-matrix VI scoring) per Methods; retain B's leiden cluster labels.
A = pd.read_csv("results/malignant_hepatocyte_vi_score.csv", index_col=0)
B = pd.read_csv("results/malig_cluster_vi_score.csv", index_col=0)
mc = B.loc[adata.obs_names[malig_mask], ["site", "leiden"]].copy()
mc["VI_score"] = A.loc[mc.index, "VI_score"].values
umap_m = umap_all[malig_mask]
print(f"    恶性肝细胞: {len(mc)}")

fig = plt.figure(figsize=(11, 3.2))
gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.1], wspace=0.45)

# A VI score by site
ax = fig.add_subplot(gs[0, 0])
groups = [mc[mc["site"] == s]["VI_score"].values for s in ["Tumor", "PVTT", "Lymph"]]
H, p = stats.kruskal(*groups)
vp = ax.violinplot(groups, positions=range(3), widths=0.8, showextrema=False)
for i, b in enumerate(vp["bodies"]):
    b.set_facecolor(list(SITE_COLORS.values())[i]); b.set_alpha(0.7)
for i, g in enumerate(groups):
    ax.hlines(np.median(g), i - 0.18, i + 0.18, color="black", lw=1.2)
ax.set_xticks(range(3)); ax.set_xticklabels(["Tumor", "PVTT", "Lymph"])
ax.set_ylabel("VI score")
ax.set_title("VI score by site", fontweight="bold")
ax.text(0.5, 0.95, fmt_p(p), transform=ax.transAxes, ha="center", fontsize=7)
tag(ax, "A")

# B UMAP colored by VI score
ax = fig.add_subplot(gs[0, 1])
order_s = np.argsort(mc["VI_score"].values)
sc_ = ax.scatter(umap_m[order_s, 0], umap_m[order_s, 1], s=1, c=mc["VI_score"].values[order_s],
                 cmap="Reds", vmin=-0.15, vmax=0.35, rasterized=True)
plt.colorbar(sc_, ax=ax, shrink=0.7, pad=0.02, label="VI score")
ax.set_title("Vascular invasion score", fontweight="bold")
ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
tag(ax, "B")

# C UMAP cluster 2 highlighted
ax = fig.add_subplot(gs[0, 2])
is_c2 = (mc["leiden"] == 2).values
ax.scatter(umap_m[~is_c2, 0], umap_m[~is_c2, 1], s=1, c="#C9C9C9", alpha=0.5, rasterized=True)
ax.scatter(umap_m[is_c2, 0], umap_m[is_c2, 1], s=1, c=C_HI, alpha=0.7, label="cluster 2 (high VI)", rasterized=True)
ax.set_title("Malignant hepatocyte clusters", fontweight="bold")
ax.legend(frameon=False, markerscale=4, loc="upper left", fontsize=6)
ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
tag(ax, "C")

# D Top 5 clusters
ax = fig.add_subplot(gs[0, 3])
top5 = mc.groupby("leiden")["VI_score"].median().sort_values(ascending=False).head(5)
for i, (cl, med) in enumerate(top5.items()):
    vals = mc[mc["leiden"] == cl]["VI_score"].values
    vp = ax.violinplot([vals], positions=[i], widths=0.8, showextrema=False)
    vp["bodies"][0].set_facecolor(C_HI if cl == 2 else "#B5D4F4")
    vp["bodies"][0].set_alpha(0.75)
    ax.hlines(np.median(vals), i - 0.18, i + 0.18, color="black", lw=1.2)
ax.set_xticks(range(5))
ax.set_xticklabels([f"c{c}" + (" (high VI)" if c == 2 else "") for c in top5.index],
                   rotation=20, ha="right", fontsize=7)
ax.set_ylabel("VI score")
ax.set_title("Top 5 clusters by VI score", fontweight="bold")
tag(ax, "D", x=-0.18)

fig.savefig("results/Fig6_single_cell_localization.png")
fig.savefig("figures/Fig6_single_cell_localization.png")

# --- Fig6 实际绘制值（A 口径），供图文同口径校验 ---
import json as _json
_fig6_manifest = {
    "vi_source_file": "results/malignant_hepatocyte_vi_score.csv (A, full-matrix scoring)",
    "by_site_median": {s: round(float(mc[mc["site"] == s]["VI_score"].median()), 4)
                       for s in ["Tumor", "PVTT", "Lymph"]},
    "by_site_mean":   {s: round(float(mc[mc["site"] == s]["VI_score"].mean()), 4)
                       for s in ["Tumor", "PVTT", "Lymph"]},
    "kruskal_P": float(p),
    "cluster2_mean_vi": round(float(mc[mc["leiden"] == 2]["VI_score"].mean()), 4),
    "top5_clusters_median": {int(c): round(float(v), 4) for c, v in top5.items()},
}
with open("results/fig6_manifest.json", "w") as _f:
    _json.dump(_fig6_manifest, _f, indent=2)
print("Fig6 重排完成（A=site 梯度 / B=UMAP VI / C=cluster2 / D=top5）")
print("Fig6 manifest 已保存（Fig6 guard）：", _fig6_manifest["by_site_median"], "cluster2=", _fig6_manifest["cluster2_mean_vi"])
