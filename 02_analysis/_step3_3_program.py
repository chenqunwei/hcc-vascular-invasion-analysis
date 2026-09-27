"""Step 3.3 修复版：Wilcoxon 程序分层 + PAGA 轨迹。"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("=" * 60)
print("Step 3.3 修复版：Wilcoxon 程序分层 + PAGA 轨迹")
print("=" * 60)

# ---------- 1. 提取肝细胞 ----------
adata = sc.read_h5ad("results/gse149614_processed_qc.h5ad")
heps = adata[adata.obs["celltype"] == "Hepatocyte"].copy()
print(f"肝细胞: {heps.shape[0]}，site: {heps.obs['site'].value_counts().to_dict()}")

up_sym = pd.read_csv("results/vascular_invasion_up_genes_symbol.csv")["symbol"].dropna().tolist()
dn_sym = pd.read_csv("results/vascular_invasion_dn_genes_symbol.csv")["symbol"].dropna().tolist()
up_hit = [g for g in up_sym if g in heps.var_names]
dn_hit = [g for g in dn_sym if g in heps.var_names]
vi_genes = up_hit + dn_hit
print(f"血管侵犯基因命中: {len(vi_genes)}")

# ---------- 2. Wilcoxon 差异检验 ----------
print("\n[2] Wilcoxon 差异检验（Normal vs Tumor / Tumor vs PVTT）...")
heps_nonormal = heps[heps.obs["site"] != "Lymph"].copy()  # Lymph 另论

# Normal vs Tumor（在全部肝细胞上做）
sc.tl.rank_genes_groups(heps, groupby="site", groups=["Tumor"], reference="Normal",
                        method="wilcoxon", key_added="NT")
res_NT = sc.get.rank_genes_groups_df(heps, group="Tumor", key="NT")
res_NT = res_NT[res_NT["names"].isin(vi_genes)].set_index("names")

# Tumor vs PVTT（在非 Lymph 肝细胞上做）
sc.tl.rank_genes_groups(heps_nonormal, groupby="site", groups=["PVTT"], reference="Tumor",
                        method="wilcoxon", key_added="TP")
res_TP = sc.get.rank_genes_groups_df(heps_nonormal, group="PVTT", key="TP")
res_TP = res_TP[res_TP["names"].isin(vi_genes)].set_index("names")

# ---------- 3. 分类：核心 vs 分期特异 ----------
print("\n[3] 分类（Wilcoxon FDR<0.05）...")
FDR = 0.05
core, stage = [], []
for g in vi_genes:
    if g not in res_NT.index or g not in res_TP.index:
        continue
    nt_p, nt_lfc = res_NT.loc[g, "pvals_adj"], res_NT.loc[g, "logfoldchanges"]
    tp_p, tp_lfc = res_TP.loc[g, "pvals_adj"], res_TP.loc[g, "logfoldchanges"]
    if nt_p < FDR and abs(nt_lfc) > 0.25:
        core.append(g)  # Normal→Tumor 已显著（核心，早期启动）
    elif nt_p >= FDR and tp_p < FDR and abs(tp_lfc) > 0.25:
        stage.append(g)  # 仅 Tumor→PVTT 显著（分期特异）

print(f"  核心连续谱程序（Normal→Tumor 显著）: {len(core)}")
print(f"  分期特异程序（仅 Tumor→PVTT 显著）: {len(stage)}")
pd.DataFrame({"gene": core}).to_csv("results/core_program_genes.csv", index=False)
pd.DataFrame({"gene": stage}).to_csv("results/stage_specific_genes.csv", index=False)

# ---------- 4. 分层热图 ----------
print("\n[4] 分层热图 ...")
plot_n = 30
plot_genes = core[:plot_n] + stage[:plot_n]
groups_labels = ["Core"] * plot_n + ["Stage-specific"] * plot_n
site_means = {}
for site in ["Normal", "Tumor", "PVTT"]:
    mask = heps.obs["site"] == site
    X_sub = heps[mask, plot_genes].X
    site_means[site] = np.asarray(X_sub.toarray()).mean(axis=0) if hasattr(X_sub, "toarray") else np.asarray(X_sub.mean(axis=0)).ravel()
X_plot = np.array([site_means[s] for s in ["Normal", "Tumor", "PVTT"]]).T
X_z = (X_plot - X_plot.mean(axis=1, keepdims=True)) / (X_plot.std(axis=1, keepdims=True) + 1e-10)

plt.rcParams.update({"font.family": "Arial", "font.size": 9, "savefig.dpi": 300})
fig, ax = plt.subplots(figsize=(3.2, 5.5))
im = ax.imshow(X_z, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
ax.set_xticks(range(3)); ax.set_xticklabels(["Normal", "Tumor", "PVTT"], rotation=30, ha="right")
ax.set_yticks([])
ax.axhline(plot_n - 0.5, color="black", lw=1.2)
ax.text(-0.55, plot_n/2, "Core", va="center", fontsize=8, fontweight="bold", rotation=90)
ax.text(-0.55, plot_n + plot_n/2, "Stage", va="center", fontsize=8, fontweight="bold", rotation=90)
ax.set_title("Program stratification\n(Wilcoxon, FDR<0.05)", fontweight="bold", fontsize=9)
fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02, label="z-score")
fig.savefig("results/FigS_program_stratification.png", bbox_inches="tight")
print("  已保存 FigS_program_stratification.png")
plt.close(fig)

# ---------- 5. PAGA 轨迹 ----------
print("\n[5] PAGA 轨迹 ...")
heps_sub = heps[heps.obs["site"].isin(["Normal", "Tumor", "PVTT"])].copy()
sc.pp.highly_variable_genes(heps_sub, n_top_genes=2000, subset=True)
sc.tl.pca(heps_sub, n_comps=30)
sc.pp.neighbors(heps_sub, n_neighbors=15, n_pcs=30)
sc.tl.leiden(heps_sub, resolution=0.4)
sc.tl.paga(heps_sub, groups="site")
sc.pl.paga(heps_sub, color="site", save="_invasion.png", show=False,
           edge_width_scale=2, node_size_scale=3)

# 各 site 的平均 VI score 沿轨迹标注
vi_by_site = {}
for s in ["Normal", "Tumor", "PVTT"]:
    mask = heps_sub.obs["site"] == s
    X_g = heps_sub[mask, vi_genes].X
    up_m = np.asarray(heps_sub[mask, up_hit].X.toarray()).mean(axis=1).mean() if hasattr(heps_sub[mask, up_hit].X, "toarray") else 0
    dn_m = np.asarray(heps_sub[mask, dn_hit].X.toarray()).mean(axis=1).mean() if hasattr(heps_sub[mask, dn_hit].X, "toarray") else 0
    vi_by_site[s] = up_m - dn_m
print(f"  各 site VI score: { {k: round(v,3) for k,v in vi_by_site.items()} }")

print("\nStep 3.3 修复版完成。")
