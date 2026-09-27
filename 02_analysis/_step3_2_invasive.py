"""Step 3.2 定位血管侵犯侵袭亚群：恶性肝细胞血管侵犯评分 + site 间对比。"""
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
print("Step 3.2  定位血管侵犯侵袭亚群（恶性肝细胞评分）")
print("=" * 60)

# ---------- 1. 加载 h5ad + 提取恶性肝细胞 ----------
print("\n[1] 加载 h5ad，提取恶性肝细胞（排除 Normal 的正常肝细胞）...")
adata = sc.read_h5ad("results/gse149614_processed_qc.h5ad")
print(f"  总细胞: {adata.shape[0]}")

# 提取 Hepatocyte，排除 Normal site（正常肝细胞）
malig = adata[(adata.obs["celltype"] == "Hepatocyte") & (adata.obs["site"] != "Normal")].copy()
print(f"  恶性肝细胞（Tumor/PVTT/Lymph）: {malig.shape[0]}")
print(f"  site 分布: {malig.obs['site'].value_counts().to_dict()}")

# ---------- 2. 加载血管侵犯基因集 ----------
print("\n[2] 加载血管侵犯基因集（symbol）...")
up_sym = pd.read_csv("results/vascular_invasion_up_genes_symbol.csv")["symbol"].dropna().tolist()
dn_sym = pd.read_csv("results/vascular_invasion_dn_genes_symbol.csv")["symbol"].dropna().tolist()

# 匹配到单细胞矩阵
up_hit = [g for g in up_sym if g in malig.var_names]
dn_hit = [g for g in dn_sym if g in malig.var_names]
print(f"  上调基因命中 {len(up_hit)}/{len(up_sym)}，下调命中 {len(dn_hit)}/{len(dn_sym)}")

# ---------- 3. score_genes 打分 ----------
print("\n[3] score_genes 打分（上调程序分 + 下调程序分）...")
sc.tl.score_genes(malig, gene_list=up_hit, score_name="up_score", ctrl_size=50, random_state=0)
sc.tl.score_genes(malig, gene_list=dn_hit, score_name="dn_score", ctrl_size=50, random_state=0)
malig.obs["VI_score"] = malig.obs["up_score"] - malig.obs["dn_score"]
print(f"  VI_score 范围: {malig.obs['VI_score'].min():.2f} ~ {malig.obs['VI_score'].max():.2f}")

# ---------- 4. site 间 VI score 对比 ----------
print("\n[4] site 间 VI score 对比（Tumor vs PVTT vs Lymph）...")
sites = ["Tumor", "PVTT", "Lymph"]
site_means = {}
for s in sites:
    sub = malig.obs[malig.obs["site"] == s]["VI_score"]
    site_means[s] = sub.mean()
    print(f"  {s}: n={len(sub)}, VI_score 均值={sub.mean():.3f}")

groups = [malig.obs[malig.obs["site"] == s]["VI_score"].values for s in sites]
groups = [g for g in groups if len(g) >= 2]
H, p = stats.kruskal(*groups)
print(f"  Kruskal-Wallis: H={H:.1f}, p={p:.2e}")

# Tumor vs PVTT 两两
_, p_tp = stats.mannwhitneyu(malig.obs[malig.obs["site"]=="Tumor"]["VI_score"],
                             malig.obs[malig.obs["site"]=="PVTT"]["VI_score"])
_, p_tl = stats.mannwhitneyu(malig.obs[malig.obs["site"]=="Tumor"]["VI_score"],
                             malig.obs[malig.obs["site"]=="Lymph"]["VI_score"])
print(f"  Tumor vs PVTT: p={p_tp:.2e}")
print(f"  Tumor vs Lymph: p={p_tl:.2e}")

# ---------- 5. 小提琴图 ----------
print("\n[5] 画小提琴图 ...")
plt.rcParams.update({"font.family": "Arial", "font.size": 9, "savefig.dpi": 300})
fig, ax = plt.subplots(figsize=(5, 4))
colors = {"Tumor": "#378ADD", "PVTT": "#D85A30", "Lymph": "#7A5195"}
vp = ax.violinplot([malig.obs[malig.obs["site"]==s]["VI_score"].values for s in sites],
                   positions=range(3), showmedians=True)
for i, s in enumerate(sites):
    vp["bodies"][i].set_facecolor(colors[s]); vp["bodies"][i].set_alpha(0.6)
ax.set_xticks(range(3)); ax.set_xticklabels(sites)
ax.set_ylabel("Vascular invasion score")
ax.set_title(f"Malignant hepatocyte VI score by site\nKruskal-Wallis P = {p:.2e}", fontweight="bold")
fig.savefig("results/FigS_VIscore_by_site.png")
print("  已保存 results/FigS_VIscore_by_site.png")
plt.close(fig)

# ---------- 6. 保存 ----------
malig.obs[["site", "patient", "up_score", "dn_score", "VI_score"]].to_csv("results/malignant_hepatocyte_vi_score.csv")
print(f"\n[6] 已保存 results/malignant_hepatocyte_vi_score.csv（{malig.shape[0]} 恶性肝细胞）")
print("\nStep 3.2 完成。")
