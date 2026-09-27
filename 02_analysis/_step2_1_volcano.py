import tarfile
import gzip
import io
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("补 Step 2.1 火山图")

# ---------- 1. 加载 GSE77509 ----------
def load_gse77509(tar_path):
    tf = tarfile.open(tar_path, "r")
    mats = {}
    for member in tf.getmembers():
        if member.name.endswith(".sf_normalized_count.txt.gz"):
            sample = member.name.split("_")[1].split(".")[0]
            content = gzip.GzipFile(fileobj=tf.extractfile(member)).read().decode()
            df = pd.read_csv(io.StringIO(content), sep="\t")
            df.columns = ["gene", sample]
            mats[sample] = df.set_index("gene")[sample]
    tf.close()
    return pd.DataFrame(mats)

print("加载 GSE77509 ...")
mat = load_gse77509("data/geo/GSE77509_GSE77509_RAW.tar")
logmat = np.log2(mat + 1)
N = [c for c in mat.columns if c[0] == "N"]
T = [c for c in mat.columns if c[0] == "T"]
P = [c for c in mat.columns if c[0] == "P"]
print(f"样本: N={len(N)} T={len(T)} P={len(P)}")

# ---------- 2. 计算 log2FC (PVTT vs 正常) + p 值 ----------
print("计算 log2FC (PVTT vs Normal) + Kruskal p 值 ...")
rows = []
for g in logmat.index:
    n_grp = logmat.loc[g, N].values
    p_grp = logmat.loc[g, P].values
    t_grp = logmat.loc[g, T].values
    log2fc = np.median(p_grp) - np.median(n_grp)
    try:
        h, p = stats.kruskal(n_grp, t_grp, p_grp)
    except Exception:
        continue
    rows.append((g, log2fc, p))
df = pd.DataFrame(rows, columns=["gene", "log2FC", "pvalue"])
df["neg_log10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
print(f"基因数: {len(df)}")

# ---------- 3. 画火山图 ----------
print("绘制火山图 ...")
sig_up = (df["pvalue"] < 0.05) & (df["log2FC"] > 0)
sig_dn = (df["pvalue"] < 0.05) & (df["log2FC"] < 0)
n_up = sig_up.sum(); n_dn = sig_dn.sum()

fig, ax = plt.subplots(figsize=(8, 7))
ax.scatter(df.loc[~sig_up & ~sig_dn, "log2FC"], df.loc[~sig_up & ~sig_dn, "neg_log10p"],
           s=6, color="#D3D1C7", alpha=0.6, label="NS", rasterized=True)
ax.scatter(df.loc[sig_up, "log2FC"], df.loc[sig_up, "neg_log10p"],
           s=6, color="#E24B4A", alpha=0.7, label=f"Up in PVTT (n={n_up})", rasterized=True)
ax.scatter(df.loc[sig_dn, "log2FC"], df.loc[sig_dn, "neg_log10p"],
           s=6, color="#378ADD", alpha=0.7, label=f"Down in PVTT (n={n_dn})", rasterized=True)

ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)
ax.axvline(0, color="gray", linestyle="--", linewidth=0.8)

# 标注已知基因（若存在）
known = ["AFP", "GPC3", "MYBL2", "PLK1", "E2F2", "CD38", "CPS1", "PON1"]
# 需要 gene symbol 映射，此处跳过（Ensembl ID），仅标注若匹配
ax.set_xlabel("log2 fold change (PVTT vs Normal)", fontsize=12)
ax.set_ylabel("-log10(p value)", fontsize=12)
ax.set_title("Vascular invasion gradient: PVTT vs Normal (GSE77509)", fontsize=13, fontweight="bold")
ax.legend(frameon=False, fontsize=10, loc="upper left")
ax.set_xlim(-6, 6)
plt.tight_layout()
plt.savefig("results/volcano_gse77509.png", dpi=150, bbox_inches="tight")
print(f"已保存 results/volcano_gse77509.png")
print(f"显著上调 {n_up}，显著下调 {n_dn}（p<0.05）")
plt.close()
print("\n火山图补充完成。")
