import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import chi2_contingency, mannwhitneyu

print("=" * 60)
print("Step 2.3  血管侵犯分子分型（共识聚类）")
print("=" * 60)

# ---------- 1. 加载表达矩阵，取 top 变异基因 ----------
print("\n[1] 加载血管侵犯基因表达矩阵 ...")
expr = pd.read_csv("results/tcga_lihc_vi_genes_expr.csv", index_col=0)
print(f"    矩阵: {expr.shape[0]} 基因 x {expr.shape[1]} 样本")

logexpr = np.log2(expr + 1)

# 取 top 变异基因（方差排序）
var = logexpr.var(axis=1).sort_values(ascending=False)
top_genes = var.index[:1000]
X = logexpr.loc[top_genes].T  # 样本 x 基因
print(f"    取 top {len(top_genes)} 高变异基因，X: {X.shape}")

# ---------- 2. 共识聚类 ----------
def consensus_matrix(X, k, n_iter=50, frac_sample=0.8, frac_feature=0.8, seed=0):
    n = X.shape[0]
    M = np.zeros((n, n))
    count = np.zeros((n, n))
    rng = np.random.default_rng(seed)
    for _ in range(n_iter):
        idx_s = rng.choice(n, max(int(n * frac_sample), k + 1), replace=False)
        idx_f = rng.choice(X.shape[1], max(int(X.shape[1] * frac_feature), 10), replace=False)
        sub = X.iloc[idx_s, idx_f].values
        km = KMeans(n_clusters=k, n_init=10, random_state=42)
        lab = km.fit_predict(sub)
        for a in range(len(idx_s)):
            for b in range(a + 1, len(idx_s)):
                if lab[a] == lab[b]:
                    M[idx_s[a], idx_s[b]] += 1
                    M[idx_s[b], idx_s[a]] += 1
                count[idx_s[a], idx_s[b]] += 1
                count[idx_s[b], idx_s[a]] += 1
    count[count == 0] = 1
    return M / count

def pac(cons):
    triu = cons[np.triu_indices_from(cons, k=1)]
    return np.mean((triu > 0.1) & (triu < 0.9))

def cdf_area(cons):
    triu = cons[np.triu_indices_from(cons, k=1)]
    vals = np.linspace(0, 1, 101)
    cdf = np.array([np.mean(triu <= v) for v in vals])
    return np.trapezoid(cdf, vals)

print("\n[2] 共识聚类（k = 2..5，各 50 次子采样）...")
results = {}
for k in range(2, 6):
    cons = consensus_matrix(X, k, n_iter=50)
    p = pac(cons)
    area = cdf_area(cons)
    # silhouette（用完整数据 k-means 标签）
    lab_full = KMeans(n_clusters=k, n_init=10, random_state=42).fit_predict(X.values)
    sil = silhouette_score(X.values, lab_full)
    results[k] = {"cons": cons, "pac": p, "area": area, "sil": sil, "lab": lab_full}
    print(f"    k={k}: PAC={p:.3f}  CDF面积={area:.3f}  silhouette={sil:.3f}")

# ---------- 3. 选最优 k ----------
print("\n[3] 评估最优 k ...")
ks = list(results.keys())
pac_vals = [results[k]["pac"] for k in ks]
sil_vals = [results[k]["sil"] for k in ks]
area_vals = [results[k]["area"] for k in ks]

# Δ area（相邻 k 的 CDF 面积差，越大越好）
delta_area = {k: area_vals[i] - area_vals[i-1] for i, k in enumerate(ks[1:], start=1)}

print("    ΔCDF面积（k vs k-1）:", {k: round(v, 3) for k, v in delta_area.items()})
best_k = min(ks, key=lambda k: results[k]["pac"])
print(f"    PAC 最小 → k={best_k}")

# ---------- 4. 最终分型（用共识矩阵层次聚类） ----------
print(f"\n[4] 用共识矩阵层次聚类得最终分型（k={best_k}）...")
cons_best = results[best_k]["cons"]
# 距离 = 1 - 共识
dist = 1 - cons_best
np.fill_diagonal(dist, 0)
Z = linkage(dist[np.triu_indices_from(dist, k=1)], method="average")
# 用距离矩阵的 condensed 形式
from scipy.spatial.distance import squareform
condensed = squareform(dist)
Z = linkage(condensed, method="average")
final_lab = fcluster(Z, t=best_k, criterion="maxclust")
X_final = X.copy()
X_final["cluster"] = final_lab
print(f"    分型结果分布:")
print(X_final["cluster"].value_counts().sort_index().to_string())

# ---------- 5. 分型与 MVI / VI score 关联 ----------
print("\n[5] 分型与临床关联 ...")
clin = pd.read_csv("results/tcga_lihc_clinical.csv")
clin = clin[clin["mvi"].isin([0, 1])]
ssgsea = pd.read_csv("results/tcga_lihc_ssgsea_vi_score.csv")[["patient_id", "ssgsea_vi", "meanz_vi_score"]]

cluster_df = X_final[["cluster"]].reset_index().rename(columns={"index": "patient_id"})
merged = cluster_df.merge(clin[["patient_id", "mvi", "recurrence", "death", "dfs_months", "os_months"]], on="patient_id", how="inner")
merged = merged.merge(ssgsea, on="patient_id", how="left")
# 按 ssgsea_vi 中位数升序重标号，兼容既有 VI_ORDER={1:Low,3:Int,2:High}
# （lowest VI→1, middle→3, highest→2），确保语义恒定、与随机种子无关，杜绝换种子后标签反转
_rank = merged.groupby("cluster")["ssgsea_vi"].median().sort_values().index.tolist()
_relabel = {_rank[0]: 1, _rank[1]: 3, _rank[2]: 2}
merged["cluster"] = merged["cluster"].map(_relabel)
print(f"    有分型 + MVI 标注: {len(merged)} 例")

# 各亚型的 MVI 比例
for c in sorted(merged["cluster"].unique()):
    sub = merged[merged["cluster"] == c]
    mvi_rate = sub["mvi"].mean()
    vi_mean = sub["ssgsea_vi"].mean()
    print(f"    亚型 {c}: n={len(sub)}, MVI+比例={mvi_rate:.2%}, ssGSEA VI均值={vi_mean:.1f}")

# 卡方检验（分型 vs MVI）
ct = pd.crosstab(merged["cluster"], merged["mvi"])
chi2, cp, dof, _ = chi2_contingency(ct)
print(f"\n    分型 vs MVI 卡方检验: χ²={chi2:.2f}, p={cp:.4f}")

# ---------- 6. 保存 ----------
print("\n[6] 保存结果 ...")
out = merged[["patient_id", "cluster", "mvi", "recurrence", "death", "dfs_months", "os_months", "ssgsea_vi", "meanz_vi_score"]]
out.to_csv("results/tcga_lihc_clusters.csv", index=False)
np.save("results/consensus_matrix_best.npy", cons_best)
# 同时保存 k=2..5 共识矩阵数组（供 _figs_core.py 的 Fig2A/B 使用，需与分型同源一致）
np.save("results/consensus_matrices.npy", np.array([results[k]["cons"] for k in range(2, 6)]))
print(f"    已保存 results/tcga_lihc_clusters.csv（{len(out)} 例）")
print(f"    已保存 consensus_matrix_best.npy")
print("\nStep 2.3 完成。")
