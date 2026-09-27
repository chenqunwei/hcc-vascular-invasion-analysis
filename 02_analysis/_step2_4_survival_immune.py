import json
import time
import urllib.request
import numpy as np
import pandas as pd
from scipy import stats
import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

print("=" * 60)
print("Step 2.4  预后（KM + log-rank）+ 免疫微环境")
print("=" * 60)

# ---------- 1. 加载数据 ----------
clusters = pd.read_csv("results/tcga_lihc_clusters.csv")
expr = pd.read_csv("results/tcga_lihc_vi_genes_expr.csv", index_col=0)

# cluster 重映射：按血管侵犯程度排序（1=低, 3=中, 2=高 → Low/Int/High）
vi_order = {1: "Low VI", 3: "Intermediate VI", 2: "High VI"}
clusters["subtype"] = clusters["cluster"].map(vi_order)
order = ["Low VI", "Intermediate VI", "High VI"]
print(f"亚型分布:\n{clusters['subtype'].value_counts().to_string()}")

# ---------- 2. KM 预后 ----------
print("\n[2] KM 生存分析（OS + DFS）...")

def km_analysis(df, time_col, event_col, title, out_png):
    df = df.dropna(subset=[time_col, event_col, "subtype"]).copy()
    df = df[df[time_col] > 0]
    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(7, 6))
    colors = {"Low VI": "#1D9E75", "Intermediate VI": "#EF9F27", "High VI": "#D85A30"}
    for st in order:
        sub = df[df["subtype"] == st]
        if len(sub) < 5:
            continue
        kmf.fit(sub[time_col], sub[event_col], label=f"{st} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, color=colors[st], linewidth=2)
    # log-rank（多组）
    try:
        res = multivariate_logrank_test(df[time_col], df["subtype"], df[event_col])
        p = res.p_value
    except Exception as e:
        p = np.nan
        print(f"  log-rank 失败: {e}")
    ax.set_title(f"{title}\nlog-rank p = {p:.4f}", fontsize=12, fontweight="bold")
    ax.set_xlabel("Months", fontsize=11)
    ax.set_ylabel("Survival probability", fontsize=11)
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=10)
    ax.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  {title}: log-rank p = {p:.4f}，已保存 {out_png}")
    return p

os_p = km_analysis(clusters, "os_months", "death", "Overall survival by VI subtype", "results/km_os_subtype.png")
dfs_p = km_analysis(clusters, "dfs_months", "recurrence", "Recurrence-free survival by VI subtype", "results/km_dfs_subtype.png")

# 补充：高 vs 低 两两对比（更有功效）
print("\n  高 vs 低 两两对比:")
for time_col, event_col, name in [("os_months", "death", "OS"), ("dfs_months", "recurrence", "DFS")]:
    sub = clusters.dropna(subset=[time_col, event_col, "subtype"])
    sub = sub[sub[time_col] > 0]
    hl = sub[sub["subtype"].isin(["Low VI", "High VI"])]
    if len(hl) > 0:
        res = multivariate_logrank_test(hl[time_col], hl["subtype"], hl[event_col])
        print(f"  {name} 高 vs 低: n={len(hl)}, p = {res.p_value:.4f}")

# ---------- 3. 免疫微环境 ----------
print("\n[3] 免疫微环境分析 ...")

# 免疫检查点基因（直接比较三型表达）
checkpoints = ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1", "CD8A"]

# 免疫细胞 marker 基因集（标准 marker）
immune_sets = {
    "CD8_T": ["CD8A", "CD8B", "GZMB", "PRF1", "GZMA", "GNLY"],
    "CD4_T": ["CD4", "IL7R", "CD3D", "CD3E"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"],
    "NK": ["NKG7", "KLRD1", "KLRK1", "NCR1", "FCGR3A"],
    "B_cell": ["CD19", "MS4A1", "CD79A", "CD79B"],
    "Macrophage": ["CD68", "CD163", "CSF1R", "ITGAM"],
    "DC": ["ITGAX", "CD1C", "CLEC9A", "FLT3"],
    "Neutrophil": ["FCGR3B", "S100A8", "S100A9", "CXCR2"],
}

# 表达矩阵 entrez -> symbol
def entrez_to_symbol(ids, batch=1000):
    mapping = {}
    for i in range(0, len(ids), batch):
        b = ids[i:i+batch]
        body = json.dumps({"q": ",".join(b), "scopes": "entrezgene",
                           "fields": "symbol", "species": "human"}).encode()
        req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            res = json.load(r)
        for it in res:
            q = it.get("query"); sym = it.get("symbol")
            if q and sym:
                mapping[q] = sym
        time.sleep(0.3)
    return mapping

print("  entrez -> symbol ...")
entrez_ids = [str(i) for i in expr.index]
ent2sym = entrez_to_symbol(entrez_ids)
expr_sym = expr.copy()
expr_sym.index = [ent2sym.get(g, g) for g in entrez_ids]
expr_sym = expr_sym[~expr_sym.index.duplicated(keep="first")]

# ssGSEA 免疫细胞分数
print("  ssGSEA 免疫细胞分数 ...")
gene_sets = {k: [g for g in v if g in expr_sym.index] for k, v in immune_sets.items()}
gene_sets = {k: v for k, v in gene_sets.items() if len(v) >= 3}
ss = gp.ssgsea(data=expr_sym, gene_sets=gene_sets, outdir=None,
               min_size=3, max_size=1000, sample_norm_method="rank",
               no_plot=True, threads=1)
if "Term" in ss.res2d.columns:
    imm_mat = ss.res2d.pivot(index="Name", columns="Term", values="ES")
else:
    imm_mat = ss.res2d.set_index("Name")
imm_mat.index = [str(i)[:12] for i in imm_mat.index]
imm_mat = imm_mat.groupby(level=0).mean()

# 合并到分型
imm_df = clusters.merge(imm_mat, left_on="patient_id", right_index=True, how="left")

# 免疫检查点基因（从表达矩阵取）
for ck in checkpoints:
    if ck in expr_sym.index:
        imm_df[ck] = imm_df["patient_id"].map(expr_sym.loc[ck])

# 比较三型免疫差异（Kruskal-Wallis）
print("\n  免疫细胞分数三型比较（Kruskal-Wallis）:")
immune_cols = list(gene_sets.keys())
results = []
for col in immune_cols + [c for c in checkpoints if c in imm_df.columns]:
    groups = [imm_df[imm_df["subtype"] == st][col].dropna() for st in order]
    groups = [g for g in groups if len(g) >= 3]
    if len(groups) < 2:
        continue
    h, p = stats.kruskal(*groups)
    means = {st: imm_df[imm_df["subtype"] == st][col].mean() for st in order}
    results.append({"feature": col, "p_value": p, **{f"mean_{st}": means.get(st) for st in order}})
    sig = "**" if p < 0.01 else ("*" if p < 0.05 else "")
    print(f"    {col}: p={p:.4f} {sig} | Low={means.get('Low VI',0):.2f} Int={means.get('Intermediate VI',0):.2f} High={means.get('High VI',0):.2f}")

res_df = pd.DataFrame(results)
res_df.to_csv("results/immune_subtype_comparison.csv", index=False)

# 保存免疫分数表
imm_df[["patient_id", "subtype", "ssgsea_vi"] + immune_cols + [c for c in checkpoints if c in imm_df.columns]].to_csv(
    "results/tcga_lihc_immune.csv", index=False)

# ---------- 4. 保存预后汇总 ----------
summary = {
    "OS_logrank_p_3group": os_p,
    "DFS_logrank_p_3group": dfs_p,
}
pd.DataFrame([summary]).to_csv("results/survival_summary.csv", index=False)
print(f"\n[4] 已保存: results/km_os_subtype.png, km_dfs_subtype.png, immune_subtype_comparison.csv, tcga_lihc_immune.csv, survival_summary.csv")
print("\nStep 2.4 完成。")
