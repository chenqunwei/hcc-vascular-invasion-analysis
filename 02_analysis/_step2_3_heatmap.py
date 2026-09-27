import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("补 Step 2.3 验收：共识矩阵热图 + CDF 曲线")

# ---------- 1. 数据 ----------
expr = pd.read_csv("results/tcga_lihc_vi_genes_expr.csv", index_col=0)
logexpr = np.log2(expr + 1)
var = logexpr.var(axis=1).sort_values(ascending=False)
top_genes = var.index[:1000]
X = logexpr.loc[top_genes].T
print(f"X: {X.shape}")

# ---------- 2. 共识聚类（保存所有 k） ----------
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

ks = [2, 3, 4, 5]
cons_mats = {}
print("共识聚类 k=2..5 ...")
for k in ks:
    cons_mats[k] = consensus_matrix(X, k, n_iter=50)
    print(f"  k={k} 完成")
np.save("results/consensus_matrices.npy", np.array([cons_mats[k] for k in ks]), allow_pickle=True)

# ---------- 3. 画热图 + CDF ----------
print("绘制共识矩阵热图 + CDF ...")
fig, axes = plt.subplots(2, 5, figsize=(20, 8))

# CDF 曲线（左侧一列，占 2 行合并）
ax_cdf = plt.subplot2grid((2, 5), (0, 0), rowspan=2)
colors = {2: "#888780", 3: "#D85A30", 4: "#378ADD", 5: "#1D9E75"}
for k in ks:
    cons = cons_mats[k]
    triu = cons[np.triu_indices_from(cons, k=1)]
    vals = np.linspace(0, 1, 101)
    cdf = np.array([np.mean(triu <= v) for v in vals])
    ax_cdf.plot(vals, cdf, label=f"k={k}", color=colors[k], linewidth=2)
ax_cdf.set_xlabel("consensus index", fontsize=11)
ax_cdf.set_ylabel("CDF", fontsize=11)
ax_cdf.set_title("Consensus CDF", fontsize=12, fontweight="bold")
ax_cdf.legend()
ax_cdf.grid(alpha=0.3)

# 共识矩阵热图（右侧 2×2，k=2,3,4,5）
for i, k in enumerate(ks):
    ax = plt.subplot2grid((2, 5), (i // 2, 1 + i % 2))
    cons = cons_mats[k]
    # 用层次聚类排序样本
    dist = 1 - cons
    np.fill_diagonal(dist, 0)
    Z = linkage(squareform(dist), method="average")
    order = leaves_list(Z)
    cons_sorted = cons[np.ix_(order, order)]
    im = ax.imshow(cons_sorted, cmap="Reds", vmin=0, vmax=1, aspect="auto")
    ax.set_title(f"k = {k}", fontsize=12, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])

fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02, label="consensus index")
fig.suptitle("Consensus clustering of vascular invasion program (TCGA-LIHC, n=371)", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig("results/consensus_heatmap.png", dpi=150, bbox_inches="tight")
print("已保存 results/consensus_heatmap.png")
plt.close()

# ---------- 4. 输出 PAC 汇总表 ----------
def pac(cons):
    triu = cons[np.triu_indices_from(cons, k=1)]
    return np.mean((triu > 0.1) & (triu < 0.9))
print("\nPAC 汇总:")
for k in ks:
    print(f"  k={k}: PAC={pac(cons_mats[k]):.3f}")
print("\n验收补充完成。")
