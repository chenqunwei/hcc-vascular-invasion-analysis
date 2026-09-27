"""Step 3.1 单细胞预处理：加载 GSE149614 → 质控 → 标准化 → 降维 → UMAP → 聚类 → 保存。

输入：data/geo/GSE149614_count.txt.gz（count 矩阵）+ metadata
输出：results/gse149614_processed.h5ad + UMAP 图
"""
import gc
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("=" * 60)
print("Step 3.1  GSE149614 单细胞预处理")
print("=" * 60)

COUNT = "data/geo/GSE149614_count.txt.gz"
META = "data/geo/GSE149614_GSE149614_HCC.metadata.updated.txt.gz"

# ---------- 1. 探测格式并加载 count ----------
print("\n[1] 探测 count 矩阵格式 ...")
import gzip as _gzip
with _gzip.open(COUNT, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    line2 = f.readline().rstrip("\n").split("\t")
print(f"  表头: {len(header)} 列（细胞数+1），第二行首列(基因): {line2[0]}")

# 手工逐行读取 + COO 稀疏构建（绕过 pandas 超宽列瓶颈，只存非零）
print("\n  逐行读取 count 矩阵，构建 COO 稀疏矩阵（只存非零）...")
from scipy.sparse import coo_matrix
rows, cols, vals = [], [], []
gene_names = []
cell_names = header  # 表头无基因列名，第一列起即 71915 个细胞
n_cells = len(cell_names)
with _gzip.open(COUNT, "rt") as f:
    f.readline()  # 跳过表头
    for i, line in enumerate(f):
        parts = line.rstrip("\n").split("\t")
        gene_names.append(parts[0])
        for j in range(1, len(parts)):
            v = parts[j]
            if v != "0" and v != "":
                rows.append(i)
                cols.append(j - 1)
                vals.append(float(v))
        if (i + 1) % 5000 == 0:
            print(f"    已读 {i+1} 个基因，非零 {len(vals)} 个", flush=True)

X = coo_matrix((vals, (rows, cols)), shape=(len(gene_names), n_cells)).tocsr().astype(np.float32)
del rows, cols, vals
gc.collect()
print(f"  count 稀疏矩阵: {X.shape[0]} 基因 x {X.shape[1]} 细胞，非零率 {X.nnz/(X.shape[0]*X.shape[1])*100:.1f}%")

# 转 AnnData（细胞 x 基因）
adata = sc.AnnData(X.T.tocsr())
adata.var_names = gene_names
adata.obs_names = cell_names
del X
gc.collect()
print(f"  AnnData: {adata.shape[0]} 细胞 x {adata.shape[1]} 基因")

# ---------- 2. 合并 metadata ----------
print("\n[2] 合并 metadata ...")
meta = pd.read_csv(META, sep="\t", compression="gzip")
meta = meta.set_index("Cell")
# 对齐 obs_names
meta = meta.reindex(adata.obs_names)
for col in ["sample", "site", "patient", "stage", "virus", "celltype", "res.3"]:
    if col in meta.columns:
        adata.obs[col] = meta[col].values
print(f"  celltype 分布:\n{adata.obs['celltype'].value_counts().to_string()}")
print(f"  site 分布:\n{adata.obs['site'].value_counts().to_string()}")

# ---------- 3. 质控 ----------
print("\n[3] 质控 ...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(f"  质控后: {adata.shape[0]} 细胞 x {adata.shape[1]} 基因")

# ---------- 4. 标准化 ----------
print("\n[4] 标准化（normalize + log1p）...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False)
print(f"  高变基因: {adata.var['highly_variable'].sum()}")

# ---------- 5. 降维 + 聚类 ----------
print("\n[5] PCA + 邻居 + UMAP + Leiden ...")
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
print(f"  Leiden 聚类数: {adata.obs['leiden'].nunique()}")

# ---------- 6. 保存 + 图 ----------
print("\n[6] 保存 + 绘图 ...")
adata.write("results/gse149614_processed.h5ad")
print("  已保存 results/gse149614_processed.h5ad")

sc.settings.figdir = "results"
for color in ["celltype", "site", "patient", "leiden"]:
    if color in adata.obs.columns:
        sc.pl.umap(adata, color=color, show=False, save=f"_gse149614_{color}.png")
        print(f"  已保存 UMAP_{color}.png")

# ---------- 7. 统计输出 ----------
print("\n[7] 汇总:")
print(adata)
print("\ncelltype x site 交叉表:")
print(pd.crosstab(adata.obs["celltype"], adata.obs["site"]).to_string())
print("\nStep 3.1 完成。")
