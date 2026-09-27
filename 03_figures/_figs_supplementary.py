# -*- coding: utf-8 -*-
"""发表级补充图重绘（统一风格，与主图一致）：
FigS1 单细胞定位（A UMAP cluster 高亮 / B UMAP VI score / C VI by site / D top clusters）
FigS2 四层次轨迹（A up 热图 / B up 均值 / C dn 热图 / D dn 均值）
Fig6 TIDE（A TIDE / B Dysfunction / C Exclusion / D IFNG）
修复：P=0.00e+00 格式、x 轴刻度错位、无 panel 标签、风格不统一。
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 9,
    "axes.linewidth": 0.6, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})
SUB_COLORS = {"Low VI": "#1D9E75", "Intermediate VI": "#EF9F27", "High VI": "#D85A30"}
C_HI = "#D85A30"
SUB_ORDER = ["Low VI", "Intermediate VI", "High VI"]
SITES = ["Normal", "Tumor", "PVTT", "Lymph"]
SITE_COLORS = {"Tumor": "#85B7EB", "PVTT": "#D85A30", "Lymph": "#7F77DD"}

def fmt_p(p):
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}"

def tag(ax, letter, x=-0.12):
    ax.text(x, 1.08, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

print("[1] 读取 h5ad（全量，约 1-2 分钟）...")
adata = sc.read_h5ad("results/gse149614_processed_qc.h5ad")
umap_all = adata.obsm["X_umap"]
print(f"    {adata.shape}")

# 恶性肝细胞（与 _step3_2 一致）
malig_mask = ((adata.obs["celltype"] == "Hepatocyte") & (adata.obs["site"] != "Normal")).values
mc = pd.read_csv("results/malig_cluster_vi_score.csv", index_col=0)
mc = mc.loc[adata.obs_names[malig_mask]]
umap_m = umap_all[malig_mask]
print(f"    恶性肝细胞: {len(mc)}")

# 肝细胞（含 Normal）用于轨迹
heps = adata[adata.obs["celltype"] == "Hepatocyte"].copy()
up_sym = pd.read_csv("results/vascular_invasion_up_genes_symbol.csv")["symbol"].dropna().tolist()
dn_sym = pd.read_csv("results/vascular_invasion_dn_genes_symbol.csv")["symbol"].dropna().tolist()
up_hit = [g for g in up_sym if g in heps.var_names]
dn_hit = [g for g in dn_sym if g in heps.var_names]
print(f"[2] 轨迹基因: up {len(up_hit)}/{len(up_sym)}, dn {len(dn_hit)}/{len(dn_sym)}")

def site_z(genes):
    idx = [heps.var_names.get_loc(g) for g in genes]
    sub = heps.X[:, idx]
    sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
    masks = {s: (heps.obs["site"] == s).values for s in SITES}
    mat = np.column_stack([sub[masks[s]].mean(axis=0) for s in SITES])
    mu, sd = mat.mean(1, keepdims=True), mat.std(1, keepdims=True)
    sd[sd == 0] = 1.0
    return (mat - mu) / sd

print("[3] 重算四 site z-score ...")
z_up, z_dn = site_z(up_hit), site_z(dn_hit)
pat_up = pd.read_csv("results/trajectory_up_patterns.csv", index_col=0)["pattern"].values
pat_dn = pd.read_csv("results/trajectory_dn_patterns.csv", index_col=0)["pattern"].values

# ============ FigS1 ============
print("[4] FigS1 单细胞定位 ...")
fig = plt.figure(figsize=(11, 3.2))
gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.1], wspace=0.45)

ax = fig.add_subplot(gs[0, 0])
is_c2 = (mc["leiden"] == 2).values
ax.scatter(umap_m[~is_c2, 0], umap_m[~is_c2, 1], s=1, c="#C9C9C9", alpha=0.5, rasterized=True)
ax.scatter(umap_m[is_c2, 0], umap_m[is_c2, 1], s=1, c=C_HI, alpha=0.7, label="cluster 2 (high VI)", rasterized=True)
ax.set_title("Malignant hepatocyte clusters", fontweight="bold")
ax.legend(frameon=False, markerscale=4, loc="upper left", fontsize=6)
ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
tag(ax, "A")

ax = fig.add_subplot(gs[0, 1])
order_s = np.argsort(mc["VI_score"].values)
sc_ = ax.scatter(umap_m[order_s, 0], umap_m[order_s, 1], s=1, c=mc["VI_score"].values[order_s],
                 cmap="Reds", vmin=-0.6, vmax=0.6, rasterized=True)
plt.colorbar(sc_, ax=ax, shrink=0.7, pad=0.02, label="VI score")
ax.set_title("Vascular invasion score", fontweight="bold")
ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
tag(ax, "B")

ax = fig.add_subplot(gs[0, 2])
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
tag(ax, "C")

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
fig.savefig("results/FigS1_single_cell_localization.png")
plt.close(fig)

# ============ FigS2 ============
print("[5] FigS2 轨迹 ...")
fig = plt.figure(figsize=(11.5, 3.4))
gs = fig.add_gridspec(2, 4, height_ratios=[1.6, 1], hspace=0.55, wspace=0.55)
PAT_META = {1: ("P1 tumor-onset", "#5A3E9B"), 2: ("P2 PVTT-onset", "#2A8C99"), 3: ("P3 metastasis-onset", "#9ACC28")}

for row, (zmat, pat, lab) in enumerate([(z_up, pat_up, "up"), (z_dn, pat_dn, "dn")]):
    o = np.lexsort((zmat[:, 1], pat))
    zm, pm = zmat[o], pat[o]
    ax = fig.add_subplot(gs[row, 0:2])
    cmap = LinearSegmentedColormap.from_list("z", ["#FFFFFF", "#F9C0AA", "#D85A30", "#7A1E0D"])
    im = ax.imshow(zm, aspect="auto", cmap=cmap, vmin=-2, vmax=2, interpolation="nearest")
    ax.set_xticks(range(4)); ax.set_xticklabels(SITES, fontsize=7)
    ax.set_yticks([])
    ax.set_title(f"{lab.capitalize()}-regulated program ({len(zm)} genes, k-means k=3)", fontweight="bold", fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02, label="z-score")
    tag(ax, "A" if row == 0 else "C")

    ax = fig.add_subplot(gs[row, 2:4])
    for c in [1, 2, 3]:
        name, color = PAT_META[c]
        ax.plot(range(4), zm[pm == c].mean(axis=0), "o-", color=color, lw=1.6,
                label=f"{name} (n={int((pm == c).sum())})")
    ax.axhline(0, color="gray", ls="--", lw=0.6)
    ax.set_xticks(range(4)); ax.set_xticklabels(SITES)
    ax.set_ylabel("Mean z-score")
    ax.set_title(f"{lab.capitalize()}-regulated mean trajectories", fontweight="bold", fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc="best")
    tag(ax, "B" if row == 0 else "D")

fig.savefig("results/FigS2_trajectory_patterns.png")
plt.close(fig)

# ============ Fig6 ============
print("[6] Fig6 TIDE 重绘 ...")
tide = pd.read_csv("results/tcga_lihc_tide.csv")
fig, axes = plt.subplots(1, 4, figsize=(12.5, 3.2))
for ax, m, ttl, lt in zip(axes, ["TIDE", "Dysfunction", "Exclusion", "IFNG"],
                          ["TIDE score", "T-cell dysfunction", "T-cell exclusion", "IFNG (inflamed)"],
                          ["A", "B", "C", "D"]):
    data = [tide[tide["subtype"] == s][m].dropna().values for s in SUB_ORDER]
    bp = ax.boxplot(data, tick_labels=SUB_ORDER, patch_artist=True, widths=0.6)
    for patch, s in zip(bp["boxes"], SUB_ORDER):
        patch.set_facecolor(SUB_COLORS[s]); patch.set_alpha(0.6)
    for w in bp["whiskers"]: w.set_color("gray")
    for c in bp["caps"]: c.set_color("gray")
    for med in bp["medians"]: med.set_color("black")
    _, p = stats.kruskal(*data)
    ax.set_title(ttl, fontweight="bold")
    ax.text(0.5, 0.95, fmt_p(p), transform=ax.transAxes, ha="center", fontsize=7)
    ax.set_xticklabels(SUB_ORDER, rotation=30, ha="right")
    tag(ax, lt, x=-0.18)
axes[0].set_ylabel("Score")
fig.tight_layout()
fig.savefig("results/Fig6_TIDE_immunotherapy.png")
plt.close(fig)
print("\n完成：FigS1/FigS2/Fig6 重绘（统一风格 + panel 标签 + P 值格式修复）。")
