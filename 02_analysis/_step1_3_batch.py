"""Step 1.3 批次校正：GSE77509 + GSE69164 合并 → PCA → ComBat → PCA 验收。"""
import io
import gzip
import json
import tarfile
import time
import urllib.request
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("=" * 60)
print("Step 1.3  批次校正（ComBat）+ PCA 验收")
print("=" * 60)

# ---------- 1. 加载两数据集 ----------
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

def symbol_to_ensembl(symbols):
    mapping = {}
    symbols = list(symbols)
    for i in range(0, len(symbols), 1000):
        batch = symbols[i:i+1000]
        body = json.dumps({"q": ",".join(batch), "scopes": "symbol",
                           "fields": "ensembl.gene", "species": "human"}).encode()
        req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            res = json.load(r)
        for item in res:
            q = item.get("query"); ens = item.get("ensembl", {})
            g = ens.get("gene") if isinstance(ens, dict) else (ens[0].get("gene") if isinstance(ens, list) and ens else None)
            if q and g:
                mapping[q] = g
        time.sleep(0.3)
    return mapping

print("[1] 加载 GSE77509 + GSE69164 ...")
gse77509 = load_gse77509("data/geo/GSE77509_GSE77509_RAW.tar")
with gzip.open("data/geo/GSE69164_GSE69164_Whole_33_sample_FPKM.xlsx.gz", "rb") as f:
    g69164_raw = pd.read_excel(io.BytesIO(f.read()), engine="openpyxl").set_index("Gene")
sample_cols = [c for c in g69164_raw.columns if any(k in c for k in ["_ANT", "_HCC", "_PVTT"])]
g69164_raw = g69164_raw[sample_cols].apply(pd.to_numeric, errors="coerce")

print(f"  GSE77509: {gse77509.shape}（Ensembl）")
print(f"  GSE69164: {g69164_raw.shape}（Symbol）")

# GSE69164 Symbol -> Ensembl
sym2ens = symbol_to_ensembl(g69164_raw.index)
g69164 = g69164_raw.rename(index=lambda s: sym2ens.get(s))
g69164 = g69164[g69164.index.notna()]
g69164 = g69164[~g69164.index.duplicated(keep="first")]
print(f"  GSE69164 转 Ensembl 后: {g69164.shape}")

# ---------- 2. 合并共同基因 ----------
print("\n[2] 取共同基因，合并 ...")
common = gse77509.index.intersection(g69164.index)
print(f"  共同基因: {len(common)}")

merged = pd.concat([gse77509.loc[common], g69164.loc[common]], axis=1)
logexpr = np.log2(merged + 1)
batch = np.array(["GSE77509"] * gse77509.shape[1] + ["GSE69164"] * g69164.shape[1])
print(f"  合并矩阵: {logexpr.shape}（基因 x 样本），batch: GSE77509={gse77509.shape[1]}, GSE69164={g69164.shape[1]}")

# ---------- 3. PCA（校正前） ----------
print("\n[3] PCA（校正前）...")
pca = PCA(n_components=2)
pca_before = pca.fit_transform(logexpr.T.values)

# ---------- 4. ComBat 校正 ----------
print("[4] ComBat 批次校正 ...")
adata = sc.AnnData(logexpr.T.values)
adata.obs["batch"] = batch
adata.var_names = list(logexpr.index)
sc.pp.combat(adata, key="batch")
logexpr_corrected = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
print(f"  校正完成")

# ---------- 5. PCA（校正后） ----------
print("[5] PCA（校正后）...")
pca_after = pca.fit_transform(logexpr_corrected.values)

# ---------- 6. 画图对比 ----------
print("[6] 画 PCA 对比图 ...")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
colors = {"GSE77509": "#D85A30", "GSE69164": "#378ADD"}
for ax, data, title in [(axes[0], pca_before, "Before ComBat"), (axes[1], pca_after, "After ComBat")]:
    for b in ["GSE77509", "GSE69164"]:
        idx = batch == b
        ax.scatter(data[idx, 0], data[idx, 1], s=20, c=colors[b], label=b, alpha=0.7)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("PC1", fontsize=11)
    ax.set_ylabel("PC2", fontsize=11)
    ax.legend(frameon=False, fontsize=10)
plt.tight_layout()
plt.savefig("results/batch_correction_pca.png", dpi=150, bbox_inches="tight")
plt.close()
print("  已保存 results/batch_correction_pca.png")

# ---------- 7. 定量评估批次效应 ----------
print("\n[7] 定量评估批次效应 ...")
# 用 PC1/PC2 与 batch 的关联（Kruskal 或 t 检验）
from scipy import stats
for pc_name, before, after in [("PC1", pca_before[:,0], pca_after[:,0]), ("PC2", pca_before[:,1], pca_after[:,1])]:
    g1_b = before[batch=="GSE77509"]; g2_b = before[batch=="GSE69164"]
    g1_a = after[batch=="GSE77509"]; g2_a = after[batch=="GSE69164"]
    p_before = stats.mannwhitneyu(g1_b, g2_b).pvalue
    p_after = stats.mannwhitneyu(g1_a, g2_a).pvalue
    print(f"  {pc_name}: 校正前 p={p_before:.2e}，校正后 p={p_after:.2e}")

print("\nStep 1.3 完成。")
