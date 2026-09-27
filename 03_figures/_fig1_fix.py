# -*- coding: utf-8 -*-
"""Fig1 重画：A 火山图高亮最终核心基因集（665 up / 1369 dn），修复与正文的三阶段梯度断层。
同时修复 B panel 气泡裁切 + A 图例与 B 标签重叠。
"""
import io, gzip, tarfile, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 9,
    "axes.linewidth": 0.6, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})
C_UP, C_DN, C_NS = "#D85A30", "#378ADD", "#C9C9C9"

def tag(ax, letter, x=-0.12):
    ax.text(x, 1.08, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

def _points_per_unit_x(ax):
    """x 轴：1 个数据单位 = 多少 points。
    关键：transData.transform 返回「像素」（基于 fig.dpi），而 scatter 的 s 与字号都是
    「points」。1 point = dpi/72 像素，必须除以 (dpi/72) 统一单位，否则半径/间距会被
    低估 (dpi/72) 倍（本例 fig.dpi=100 → 低估 1.389 倍，标签因此压入气泡）。"""
    dpi = ax.figure.dpi
    pix_per_unit = abs(ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0])
    return pix_per_unit / (dpi / 72.0)

def bubble_radius(ax, sizes):
    """scatter 的 s（points²）→ x 轴数据坐标下的气泡半径。"""
    return np.sqrt(np.asarray(sizes, dtype=float)) / 2 / _points_per_unit_x(ax)

def pts_to_data(ax, pts):
    """把 points 长度换算成 x 轴数据坐标（标签与气泡的固定字符间距）。"""
    return pts / _points_per_unit_x(ax)

# ============ 数据 ============
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
    n_, p_ = logmat.loc[g, N].values, logmat.loc[g, P].values
    log2fc = np.median(p_) - np.median(n_)
    try:
        _, pv = stats.kruskal(n_, logmat.loc[g, T].values, p_)
    except Exception:
        continue
    vol.append((g, log2fc, pv))
vol = pd.DataFrame(vol, columns=["gene", "log2FC", "p"]).dropna()
vol["neglogp"] = -np.log10(vol["p"].clip(lower=1e-300))

# 核心基因集（Ensembl）
up_core = set(pd.read_csv("results/vascular_invasion_up_genes.csv")["gene"])
dn_core = set(pd.read_csv("results/vascular_invasion_dn_genes.csv")["gene"])
up_hit = up_core & set(vol["gene"])
dn_hit = dn_core & set(vol["gene"])
print(f"核心基因映射到火山图: up {len(up_hit)}/{len(up_core)}, dn {len(dn_hit)}/{len(dn_core)}")

sig_up = (vol["p"] < 0.05) & (vol["log2FC"] > 0)
sig_dn = (vol["p"] < 0.05) & (vol["log2FC"] < 0)

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

# ============ 绘图 ============
fig = plt.figure(figsize=(7.2, 9.4))
gs = fig.add_gridspec(3, 1, height_ratios=[1.15, 1, 1], hspace=0.5)

# A 火山图：背景灰 + 核心基因高亮
ax = fig.add_subplot(gs[0])
in_core = vol["gene"].isin(up_hit | dn_hit)
sig_nc = (sig_up | sig_dn) & (~in_core)  # 灰色点 = 显著且非核心基因（与图注"significant"一致）
ax.scatter(vol.loc[sig_nc, "log2FC"], vol.loc[sig_nc, "neglogp"],
           s=2, c=C_NS, alpha=0.35, rasterized=True)
ax.scatter(vol[vol["gene"].isin(up_hit)]["log2FC"], vol[vol["gene"].isin(up_hit)]["neglogp"],
           s=4, c=C_UP, alpha=0.85, label=f"Core up (n={len(up_hit)})", rasterized=True)
ax.scatter(vol[vol["gene"].isin(dn_hit)]["log2FC"], vol[vol["gene"].isin(dn_hit)]["neglogp"],
           s=4, c=C_DN, alpha=0.85, label=f"Core down (n={len(dn_hit)})", rasterized=True)
ax.axhline(-np.log10(0.05), color="gray", ls="--", lw=0.6)
ax.axvline(0, color="gray", ls="--", lw=0.6)
ax.set_xlabel("Median log2FC (PVTT vs Normal)")
ax.set_ylabel("-log10(P)")
ax.set_title("Vascular invasion gradient", fontweight="bold")
leg = ax.legend(frameon=True, loc="upper left", markerscale=3)
leg.get_frame().set_alpha(0.85)
leg.get_frame().set_edgecolor("#C9C9C9")
ax.set_xlim(-5, 5)
# 漏斗注释（右下角，与图例上下错开）
ax.text(0.98, 0.02,
        f"Single-cohort significant (P<0.05): {int(sig_up.sum()):,} up / {int(sig_dn.sum()):,} down\n"
        f"Monotonic trend + replication in a second cohort\n"
        f"Vascular invasion core program: {len(up_hit)} up / {len(dn_hit):,} down",
        transform=ax.transAxes, fontsize=6.5, va="bottom", ha="right", color="#444441",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#C9C9C9", lw=0.5))
tag(ax, "A")

# B up 富集（通路名放气泡右侧，符合 dot plot 规范）
ax = fig.add_subplot(gs[1])
b = top_bp(enr_up)
sz = np.clip(b["overlap_n"] * 18, 20, 400)
scb = ax.scatter(b["ratio"], range(len(b)), s=sz, c=b["-log10padj"],
                cmap="Blues", edgecolors="black", linewidths=0.3)
ax.set_ylim(-0.7, len(b) - 0.3)
ax.set_xlim(0, max(b["ratio"]) * 3.0)
rad = bubble_radius(ax, sz)
gap = pts_to_data(ax, 4)  # 一个英文字符 ≈ 4pt（fontsize 7）
for yi, (r, term, rd) in enumerate(zip(b["ratio"], b["Term"], rad)):
    ax.text(r + rd + gap, yi, term.split(" (")[0], va="center", ha="left", fontsize=7)
ax.set_yticks([])
ax.set_xlabel("Gene ratio")
ax.set_title("Up-regulated: GO biological process", fontweight="bold")
fig.colorbar(scb, ax=ax, shrink=0.8, pad=0.02, label="−log10(adjusted P)")
tag(ax, "B", x=-0.08)

# C dn 富集（通路名放气泡右侧）
ax = fig.add_subplot(gs[2])
b = top_bp(enr_dn)
sz = np.clip(b["overlap_n"] * 18, 20, 400)
sc = ax.scatter(b["ratio"], range(len(b)), s=sz, c=b["-log10padj"],
                cmap="Blues", edgecolors="black", linewidths=0.3)
ax.set_ylim(-0.7, len(b) - 0.3)
ax.set_xlim(0, max(b["ratio"]) * 3.0)
rad = bubble_radius(ax, sz)
gap = pts_to_data(ax, 4)  # 一个英文字符 ≈ 4pt（fontsize 7）
for yi, (r, term, rd) in enumerate(zip(b["ratio"], b["Term"], rad)):
    ax.text(r + rd + gap, yi, term.split(" (")[0], va="center", ha="left", fontsize=7)
ax.set_yticks([])
ax.set_xlabel("Gene ratio")
ax.set_title("Down-regulated: GO biological process", fontweight="bold")
fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02, label="−log10(adjusted P)")
tag(ax, "C", x=-0.08)

fig.savefig("results/Fig1_vascular_program.png")
fig.savefig("figures/Fig1_vascular_program.png")
print("Fig1 重画完成（results/ + figures/）")
