"""发表级 Figure 1 + Figure 2 生成（统一风格）。"""
import io
import gzip
import tarfile
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- 统一发表级风格 ----------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "axes.linewidth": 0.6,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# 统一配色
C_UP = "#D85A30"      # 上调（暖）
C_DN = "#378ADD"      # 下调（冷）
C_NS = "#C9C9C9"
SUB_COLORS = {"Low VI": "#1D9E75", "Intermediate VI": "#EF9F27", "High VI": "#D85A30"}
SUB_ORDER = ["Low VI", "Intermediate VI", "High VI"]

print("=" * 60)
print("发表级 Figure 1（血管侵犯程序）+ Figure 2（分子分型）")
print("=" * 60)

# ============ 数据准备 ============
# 火山图数据（GSE77509 PVTT vs Normal）
def load_gse77509(tar_path):
    tf = tarfile.open(tar_path, "r")
    mats = {}
    for m in tf.getmembers():
        if m.name.endswith(".sf_normalized_count.txt.gz"):
            s = m.name.split("_")[1].split(".")[0]
            c = gzip.GzipFile(fileobj=tf.extractfile(m)).read().decode()
            df = pd.read_csv(io.StringIO(c), sep="\t")
            df.columns = ["gene", s]
            mats[s] = df.set_index("gene")[s]
    tf.close()
    return pd.DataFrame(mats)

print("加载 GSE77509 ...")
mat = load_gse77509("data/geo/GSE77509_GSE77509_RAW.tar")
logmat = np.log2(mat + 1)
N = [c for c in mat.columns if c[0] == "N"]
T = [c for c in mat.columns if c[0] == "T"]
P = [c for c in mat.columns if c[0] == "P"]
vol = []
for g in logmat.index:
    n_, t_, p_ = logmat.loc[g, N].values, logmat.loc[g, T].values, logmat.loc[g, P].values
    log2fc = np.median(p_) - np.median(n_)
    try:
        _, pv = stats.kruskal(n_, t_, p_)
    except Exception:
        continue
    vol.append((g, log2fc, pv))
vol = pd.DataFrame(vol, columns=["gene", "log2FC", "p"]).dropna()
vol["neglogp"] = -np.log10(vol["p"].clip(lower=1e-300))

# 富集数据
enr_up = pd.read_csv("results/enrichment_up.csv")
enr_dn = pd.read_csv("results/enrichment_dn.csv")
def top_bp(res, n=10):
    r = res[res["Library"] == "GO_Biological_Process_2023"].copy()
    r = r[r["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value").head(n)
    r["overlap_n"] = r["Overlap"].apply(lambda s: int(s.split("/")[0]))
    r["total_n"] = r["Overlap"].apply(lambda s: int(s.split("/")[1]))
    r["ratio"] = r["overlap_n"] / r["total_n"]
    r["-log10padj"] = -np.log10(r["Adjusted P-value"].clip(lower=1e-30))
    return r

# 分型数据
clusters = pd.read_csv("results/tcga_lihc_clusters.csv")
vi_order = {1: "Low VI", 3: "Intermediate VI", 2: "High VI"}
clusters["subtype"] = clusters["cluster"].map(vi_order)
clusters = clusters[clusters["mvi"].isin([0, 1])]

# ============ Figure 1 ============
print("生成 Figure 1 ...")
fig = plt.figure(figsize=(10, 3.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1], wspace=0.35)

# A 火山图
ax = fig.add_subplot(gs[0, 0])
sig_up = (vol["p"] < 0.05) & (vol["log2FC"] > 0)
sig_dn = (vol["p"] < 0.05) & (vol["log2FC"] < 0)
ax.scatter(vol.loc[~sig_up & ~sig_dn, "log2FC"], vol.loc[~sig_up & ~sig_dn, "neglogp"],
           s=2, c=C_NS, alpha=0.5, rasterized=True)
ax.scatter(vol.loc[sig_up, "log2FC"], vol.loc[sig_up, "neglogp"],
           s=3, c=C_UP, alpha=0.7, label=f"Up (n={sig_up.sum()})", rasterized=True)
ax.scatter(vol.loc[sig_dn, "log2FC"], vol.loc[sig_dn, "neglogp"],
           s=3, c=C_DN, alpha=0.7, label=f"Down (n={sig_dn.sum()})", rasterized=True)
ax.axhline(-np.log10(0.05), color="gray", ls="--", lw=0.6)
ax.axvline(0, color="gray", ls="--", lw=0.6)
ax.set_xlabel("log2 fold change (PVTT vs Normal)")
ax.set_ylabel("-log10(P)")
ax.set_title("Vascular invasion gradient", fontweight="bold")
ax.legend(frameon=False, loc="upper left", markerscale=3)
ax.set_xlim(-5, 5)
ax.text(-0.12, 1.08, "A", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

# B 上调富集
ax = fig.add_subplot(gs[0, 1])
b = top_bp(enr_up)
ax.scatter(b["ratio"], range(len(b)), s=b["overlap_n"] * 18, c=b["-log10padj"],
           cmap="Blues", edgecolors="black", linewidths=0.3)
ax.set_yticks(range(len(b)))
ax.set_yticklabels([t.split(" (")[0][:45] for t in b["Term"]], fontsize=7)
ax.set_xlabel("Gene ratio")
ax.set_title("Up-regulated: GO BP", fontweight="bold")
ax.text(-0.15, 1.08, "B", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

# C 下调富集
ax = fig.add_subplot(gs[0, 2])
b = top_bp(enr_dn)
sc = ax.scatter(b["ratio"], range(len(b)), s=b["overlap_n"] * 18, c=b["-log10padj"],
                cmap="Blues", edgecolors="black", linewidths=0.3)
ax.set_yticks(range(len(b)))
ax.set_yticklabels([t.split(" (")[0][:45] for t in b["Term"]], fontsize=7)
ax.set_xlabel("Gene ratio")
ax.set_title("Down-regulated: GO BP", fontweight="bold")
ax.text(-0.15, 1.08, "C", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
fig.colorbar(sc, ax=fig.axes[-1], shrink=0.7, pad=0.04, label="-log10(adj P)")

fig.savefig("results/Fig1_vascular_program.png")
print("  已保存 results/Fig1_vascular_program.png")
plt.close(fig)

# ============ Figure 2 ============
print("生成 Figure 2 ...")
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

cons_mats = np.load("results/consensus_matrices.npy", allow_pickle=True)
# cons_mats 是 k=2..5 的共识矩阵数组

fig = plt.figure(figsize=(10, 3.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1], wspace=0.4)

# A 共识热图 k=3（k 索引 1 = k=3）
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
ax.text(-0.12, 1.08, "A", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

# B 共识 CDF 曲线（k=2-5，Δarea 标注 k=3）
ax = fig.add_subplot(gs[0, 1])
cons_all = np.load("results/consensus_matrices.npy")  # k = 2..5
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
          title=f"Δ area (k=3 vs 2) = {delta3:.3f}", title_fontsize=7)
ax.text(-0.15, 1.08, "B", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

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
ax.text(-0.15, 1.08, "C", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

fig.savefig("results/Fig2_molecular_subtyping.png")
print("  已保存 results/Fig2_molecular_subtyping.png")
plt.close(fig)

print("\nFigure 1 + Figure 2 完成。")
