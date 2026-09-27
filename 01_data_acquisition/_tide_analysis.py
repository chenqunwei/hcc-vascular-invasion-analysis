# -*- coding: utf-8 -*-
"""TIDE 虚拟预测：TCGA-LIHC 三型免疫治疗响应预测 + 比较。
输入：全基因组表达(results/tcga_lihc_fullexpr.csv) + 三型标签(tcga_lihc_clusters.csv)
方法：tidepy TIDE(cancer='Other')，比较三型的 TIDE/Dysfunction/Exclusion/IFNG/CD8/CTL。
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tidepy.pred import TIDE

SUB_ORDER = ["Low VI", "Intermediate VI", "High VI"]
SUB_COLORS = {"Low VI": "#1D9E75", "Intermediate VI": "#EF9F27", "High VI": "#D85A30"}

# 统一风格
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 9,
    "axes.linewidth": 0.6, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})

print("[1] 加载数据 ...")
expr = pd.read_csv("results/tcga_lihc_fullexpr.csv", index_col=0)
expr.index = [int(x) for x in expr.index]
print(f"    表达矩阵: {expr.shape}")

clusters = pd.read_csv("results/tcga_lihc_clusters.csv")
vi_order = {1: "Low VI", 3: "Intermediate VI", 2: "High VI"}
clusters["subtype"] = clusters["cluster"].map(vi_order)
print(f"    三型样本: {clusters.subtype.value_counts().to_dict()}")

# 样本对齐（表达列是 patientId，聚类是 patient_id）
common = [c for c in expr.columns if c in set(clusters["patient_id"])]
expr_sub = expr[common]
print(f"    对齐样本: {len(common)}")

print("[2] 运行 TIDE ...")
tide = TIDE(expr_sub, cancer="Other", vthres=0.0)
print(f"    TIDE 结果: {tide.shape}")
print(tide[["TIDE", "Dysfunction", "Exclusion", "IFNG", "CD8", "Responder"]].describe().round(3).to_string())

# 合并三型标签
tide["patient_id"] = tide.index
tide = tide.merge(clusters[["patient_id", "subtype"]], on="patient_id", how="left")

print("\n[3] 三型比较 ...")
metrics = ["TIDE", "Dysfunction", "Exclusion", "IFNG", "CD8", "CTL", "MDSC", "CAF", "TAM M2", "MSI Score"]
res_rows = []
for m in metrics:
    if m not in tide.columns:
        continue
    groups = [tide[tide["subtype"] == s][m].dropna() for s in SUB_ORDER]
    if any(len(g) < 2 for g in groups):
        continue
    H, p = stats.kruskal(*groups)
    means = {s: tide[tide["subtype"] == s][m].mean() for s in SUB_ORDER}
    res_rows.append({"metric": m, "Kruskal_P": p, **{f"mean_{s}": means[s] for s in SUB_ORDER}})
    sig = "**" if p < 0.01 else ("*" if p < 0.05 else "")
    print(f"    {m:12s} p={p:.4f} {sig} | " + "  ".join(f"{s}={means[s]:.3f}" for s in SUB_ORDER))

res_df = pd.DataFrame(res_rows)
res_df.to_csv("results/tide_subtype_comparison.csv", index=False)

# Responder 比例
print("\n    Responder 比例（预测响应者）:")
resp_rate = {}
for s in SUB_ORDER:
    sub = tide[tide["subtype"] == s]
    rate = sub["Responder"].mean() * 100 if "Responder" in sub.columns else np.nan
    resp_rate[s] = rate
    print(f"      {s}: {rate:.1f}%")
# 卡方检验 Responder 分布
ct = pd.crosstab(tide["subtype"], tide["Responder"]) if "Responder" in tide.columns else None
if ct is not None and ct.shape[0] == 3:
    chi2, p_chi = stats.chi2_contingency(ct)[:2]
    print(f"      Responder 三型卡方 p = {p_chi:.4f}")

tide.to_csv("results/tcga_lihc_tide.csv", index=False)
print("\n[4] 已保存 results/tcga_lihc_tide.csv + results/tide_subtype_comparison.csv")

# 画图：TIDE / Dysfunction / Exclusion / IFNG 三型箱线
fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.4))
for ax, m, ttl in zip(axes, ["TIDE", "Dysfunction", "Exclusion", "IFNG"],
                      ["TIDE score", "T-cell dysfunction", "T-cell exclusion", "IFNG (inflamed)"]):
    data = [tide[tide["subtype"] == s][m].dropna().values for s in SUB_ORDER]
    bp = ax.boxplot(data, tick_labels=SUB_ORDER, patch_artist=True, widths=0.6)
    for patch, s in zip(bp["boxes"], SUB_ORDER):
        patch.set_facecolor(SUB_COLORS[s]); patch.set_alpha(0.6)
    for w in bp["whiskers"]: w.set_color("gray")
    for c in bp["caps"]: c.set_color("gray")
    for med in bp["medians"]: med.set_color("black")
    H, p = stats.kruskal(*data)
    ax.set_title(ttl, fontweight="bold")
    ax.text(0.5, 0.95, f"P = {p:.3f}", transform=ax.transAxes, ha="center", fontsize=7)
    ax.set_xticklabels(SUB_ORDER, rotation=30, ha="right")
axes[0].set_ylabel("Score")
fig.tight_layout()
fig.savefig("results/Fig6_TIDE_immunotherapy.png")
print("    已保存 results/Fig6_TIDE_immunotherapy.png")
print("\n完成 TIDE 虚拟预测。")
