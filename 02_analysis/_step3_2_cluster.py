"""Step 3.2 补充：恶性肝细胞内部无监督聚类 + 定位高 VI score 侵袭亚群。"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("=" * 60)
print("Step 3.2 补充：恶性肝细胞内部聚类 + 定位侵袭亚群")
print("=" * 60)

# ---------- 1. 提取恶性肝细胞，子集重聚类 ----------
print("[1] 提取恶性肝细胞，子集内重新降维聚类 ...")
adata = sc.read_h5ad("results/gse149614_processed_qc.h5ad")
malig = adata[(adata.obs["celltype"] == "Hepatocyte") & (adata.obs["site"] != "Normal")].copy()
print(f"  恶性肝细胞: {malig.shape[0]}")

# 子集上重新算高变基因 + PCA + neighbors + UMAP + Leiden
sc.pp.highly_variable_genes(malig, n_top_genes=2000, subset=True)
sc.tl.pca(malig, n_comps=30)
sc.pp.neighbors(malig, n_neighbors=15, n_pcs=30)
sc.tl.umap(malig)
sc.tl.leiden(malig, resolution=0.8)
n_clusters = malig.obs["leiden"].nunique()
print(f"  Leiden 聚类数（resolution=0.8）: {n_clusters}")

# ---------- 2. score_genes 打分 ----------
print("\n[2] 血管侵犯程序打分 ...")
up_sym = pd.read_csv("results/vascular_invasion_up_genes_symbol.csv")["symbol"].dropna().tolist()
dn_sym = pd.read_csv("results/vascular_invasion_dn_genes_symbol.csv")["symbol"].dropna().tolist()
up_hit = [g for g in up_sym if g in malig.var_names]
dn_hit = [g for g in dn_sym if g in malig.var_names]
sc.tl.score_genes(malig, gene_list=up_hit, score_name="up_score", ctrl_size=50, random_state=0)
sc.tl.score_genes(malig, gene_list=dn_hit, score_name="dn_score", ctrl_size=50, random_state=0)
malig.obs["VI_score"] = malig.obs["up_score"] - malig.obs["dn_score"]

# ---------- 3. 各 Leiden 亚群的 VI score ----------
print("\n[3] 各 Leiden 亚群 VI score 均值 ...")
cluster_vi = malig.obs.groupby("leiden")["VI_score"].mean().sort_values(ascending=False)
print(cluster_vi.round(3).to_string())

top_cluster = cluster_vi.index[0]
print(f"\n  VI score 最高的亚群: cluster {top_cluster}（均值 {cluster_vi.iloc[0]:.3f}）")

# 高 VI 亚群的 site 富集
top_cells = malig.obs[malig.obs["leiden"] == top_cluster]
site_comp = top_cells["site"].value_counts(normalize=True)
all_comp = malig.obs["site"].value_counts(normalize=True)
print(f"\n  高 VI 亚群 site 组成: {site_comp.round(3).to_dict()}")
print(f"  全体恶性肝细胞 site 组成: {all_comp.round(3).to_dict()}")
print(f"  PVTT 富集倍数: {site_comp.get('PVTT',0)/all_comp.get('PVTT',0):.2f}x")

# ---------- 4. 画图 ----------
print("\n[4] 画图 ...")
plt.rcParams.update({"font.family": "Arial", "font.size": 9, "savefig.dpi": 300})

# UMAP：Leiden 亚群着色，高 VI 亚群高亮
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
# 左：Leiden 亚群
umap = malig.obsm["X_umap"]
for cl in sorted(malig.obs["leiden"].unique(), key=lambda x: int(x)):
    mask = malig.obs["leiden"] == cl
    if cl == top_cluster:
        axes[0].scatter(umap[mask, 0], umap[mask, 1], s=3, c="#D85A30", label=f"cluster {cl} (high VI)", rasterized=True)
    else:
        axes[0].scatter(umap[mask, 0], umap[mask, 1], s=2, c="#C9C9C9", alpha=0.5, rasterized=True)
axes[0].set_title("Malignant hepatocyte clusters\n(high-VI cluster highlighted)", fontweight="bold")
axes[0].legend(frameon=False, fontsize=7)
axes[0].set_xticks([]); axes[0].set_yticks([])

# 右：VI score 着色
sc_ = axes[1].scatter(umap[:, 0], umap[:, 1], s=2, c=malig.obs["VI_score"], cmap="Reds", rasterized=True)
axes[1].set_title("Vascular invasion score", fontweight="bold")
axes[1].set_xticks([]); axes[1].set_yticks([])
fig.colorbar(sc_, ax=axes[1], shrink=0.7, pad=0.02)
fig.savefig("results/FigS_malig_clusters.png")
print("  已保存 results/FigS_malig_clusters.png")
plt.close(fig)

# 小提琴图：top 5 亚群 VI score
fig, ax = plt.subplots(figsize=(7, 4))
top5 = cluster_vi.index[:5].tolist()
data = [malig.obs[malig.obs["leiden"]==c]["VI_score"].values for c in top5]
vp = ax.violinplot(data, showmedians=True)
for i in range(len(top5)):
    vp["bodies"][i].set_facecolor("#D85A30" if top5[i]==top_cluster else "#9ECAE1")
    vp["bodies"][i].set_alpha(0.6)
ax.set_xticks(range(len(top5))); ax.set_xticklabels([f"c{c}" for c in top5])
ax.set_ylabel("VI score")
ax.set_title("Top 5 clusters by VI score", fontweight="bold")
fig.savefig("results/FigS_malig_cluster_vi.png")
print("  已保存 results/FigS_malig_cluster_vi.png")
plt.close(fig)

# ---------- 5. 保存 ----------
malig.obs[["site", "patient", "leiden", "up_score", "dn_score", "VI_score"]].to_csv("results/malig_cluster_vi_score.csv")
print(f"\n[5] 已保存 results/malig_cluster_vi_score.csv")
print("\nStep 3.2 补充完成。")
