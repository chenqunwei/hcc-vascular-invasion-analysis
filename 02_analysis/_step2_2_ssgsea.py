import json
import time
import urllib.request
import pandas as pd
import numpy as np
from scipy import stats
import gseapy as gp

print("=" * 60)
print("补做 ssGSEA 版本血管侵犯评分，验证与 mean z-score 一致性")
print("=" * 60)

# ---------- 1. 加载基因 symbol + 表达矩阵 ----------
up_sym = pd.read_csv("results/vascular_invasion_up_genes_symbol.csv")["symbol"].dropna().tolist()
dn_sym = pd.read_csv("results/vascular_invasion_dn_genes_symbol.csv")["symbol"].dropna().tolist()
print(f"[1] 基因 symbol: up={len(up_sym)} dn={len(dn_sym)}")

print("[2] 加载表达矩阵，entrez index 转 symbol ...")
expr = pd.read_csv("results/tcga_lihc_vi_genes_expr.csv", index_col=0)
entrez_ids = [str(i) for i in expr.index]
print(f"    矩阵: {expr.shape[0]} 基因 x {expr.shape[1]} 样本")

# entrez -> symbol（mygene）
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

ent2sym = entrez_to_symbol(entrez_ids)
expr.index = [ent2sym.get(g, g) for g in entrez_ids]
# 去重 symbol（取均值）
expr = expr[~expr.index.duplicated(keep="first")]
print(f"    转 symbol 后: {expr.shape[0]} 基因（映射 {len(ent2sym)}/{len(entrez_ids)}）")

up_in = [g for g in up_sym if g in expr.index]
dn_in = [g for g in dn_sym if g in expr.index]
print(f"    上调命中 {len(up_in)}/{len(up_sym)}，下调命中 {len(dn_in)}/{len(dn_sym)}")

gene_sets = {"VI_up": up_in, "VI_dn": dn_in}

# ---------- 3. ssGSEA ----------
print("[3] 运行 ssGSEA（GSEApy）...")
ss = gp.ssgsea(data=expr,
               gene_sets=gene_sets,
               outdir=None,
               min_size=5,
               max_size=10000,
               sample_norm_method="rank",
               no_plot=True,
               threads=1)
# 结果：行=样本，列=基因集
res = ss.res2d.copy()
print(f"    ssgsea 结果: {res.shape}")
print(res.head(3).to_string())

# 解析：res2d 通常是 样本 x (Term/ES/NES...) 的长表
# 转成宽表
if "Term" in res.columns:
    es_mat = res.pivot(index="Name", columns="Term", values="ES")
else:
    es_mat = res.set_index("Name")

es_up = es_mat["VI_up"] if "VI_up" in es_mat.columns else None
es_dn = es_mat["VI_dn"] if "VI_dn" in es_mat.columns else None
print(f"    ES_up 样本数: {0 if es_up is None else len(es_up)}")
print(f"    ES_dn 样本数: {0 if es_dn is None else len(es_dn)}")

ssgsea_vi = es_up - es_dn
ssgsea_vi.name = "ssgsea_vi_score"

# ---------- 4. 与 mean z-score 版对比 ----------
print("\n[4] 与 mean z-score 版对比 ...")
meanz = pd.read_csv("results/tcga_lihc_vi_score.csv")
meanz = meanz.set_index("patient_id")["vi_score"]

# 对齐样本（ssgsea 用 sample id，前12位 = patient）
ssgsea_vi.index = [str(i)[:12] for i in ssgsea_vi.index]
ssgsea_vi = ssgsea_vi.groupby(level=0).mean()

common = ssgsea_vi.index.intersection(meanz.index)
print(f"    共有患者: {len(common)}")
r, rp = stats.spearmanr(ssgsea_vi[common], meanz[common])
print(f"    Spearman 相关 = {r:.3f}（p = {rp:.2e}）")

# ---------- 5. ssGSEA 版区分 MVI ----------
print("\n[5] ssGSEA 版 VI score 区分 MVI ...")
clin = pd.read_csv("results/tcga_lihc_clinical.csv")
clin = clin[clin["mvi"].isin([0, 1])]
merged = pd.DataFrame({"ssgsea_vi": ssgsea_vi}).join(
    clin.set_index("patient_id")[["mvi", "recurrence", "death", "dfs_months", "os_months"]],
    how="inner")
merged = merged.dropna(subset=["mvi"])
print(f"    有 ssGSEA 分数 + MVI 标注: {len(merged)} 例")

mvi_pos = merged[merged["mvi"] == 1]["ssgsea_vi"]
mvi_neg = merged[merged["mvi"] == 0]["ssgsea_vi"]
u, p = stats.mannwhitneyu(mvi_pos, mvi_neg, alternative="two-sided")
auc = u / (len(mvi_pos) * len(mvi_neg))
print(f"    MVI+ (n={len(mvi_pos)}): 均值 {mvi_pos.mean():.3f}")
print(f"    MVI- (n={len(mvi_neg)}): 均值 {mvi_neg.mean():.3f}")
print(f"    Mann-Whitney p = {p:.4e}")
print(f"    AUC = {auc:.3f}")

# ---------- 6. 保存 ----------
out = merged.reset_index().rename(columns={"index": "patient_id"})
out = out[["patient_id", "ssgsea_vi", "mvi", "recurrence", "death", "dfs_months", "os_months"]]
# 合并 mean-z 版
meanz_col = meanz.rename("meanz_vi_score")
out = out.join(meanz_col, on="patient_id")
out.to_csv("results/tcga_lihc_ssgsea_vi_score.csv", index=False)
print(f"\n[6] 已保存 results/tcga_lihc_ssgsea_vi_score.csv（{len(out)} 例）")
print("\nStep 2.2-补充（ssGSEA）完成。")
