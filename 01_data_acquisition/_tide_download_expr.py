# -*- coding: utf-8 -*-
"""下载 TCGA-LIHC 全基因组 RSEM 表达（供 TIDE 分析用）。
TIDE 需 ~22842 基因表达（Dysfunction/Exclusion/MSI 签名），远超 VI 基因集(2029)。
分页查询 cBioPortal molecular-data/fetch，组装基因 x 样本矩阵。
"""
import io, json, time, sys
import requests
import pandas as pd

PROFILE = "lihc_tcga_rna_seq_v2_mrna"
CBIO_URL = f"https://www.cbioportal.org/api/molecular-profiles/{PROFILE}/molecular-data/fetch"
BATCH = 400

print("[1] 从 tidepy model.pkl 提取 TIDE 所需基因 ...")
m = pd.read_pickle(io.BytesIO(open(
    'C:/Users/Restoreharmony/.workbuddy/binaries/python/envs/default/Lib/site-packages/tidepy/data/model.pkl', 'rb').read()))
genes = set()
for k in ['Dysfunction', 'Exclusion']:
    genes.update(int(x) for x in m['tide'][k].index)
genes.update(int(x) for x in m['msi'].index)
for k, v in m['biomarkers'].items():
    genes.update(int(x) for x in v.index)
genes = sorted(genes)
print(f"    基因总数: {len(genes)}")

def cbio_fetch(entrez_ids, retry=4):
    body = {"sampleListId": PROFILE, "entrezGeneIds": entrez_ids}
    for attempt in range(retry):
        try:
            r = requests.post(CBIO_URL, json=body, timeout=180)
            if r.status_code == 200:
                data = r.json()
                if isinstance(data, list):
                    return [(d.get("entrezGeneId"), d.get("patientId"), d.get("value")) for d in data]
                return []
            time.sleep(2 * (attempt + 1))
        except Exception as e:
            print(f"      批次失败 {type(e).__name__}: {e}, 重试 {attempt+1}")
            time.sleep(3)
    return []

print("[2] 分页查询 cBioPortal ...")
all_rows = []
n_batches = (len(genes) + BATCH - 1) // BATCH
for i in range(0, len(genes), BATCH):
    batch = genes[i:i + BATCH]
    rows = cbio_fetch(batch)
    all_rows.extend(rows)
    done = (i // BATCH) + 1
    if done % 10 == 0 or done == n_batches:
        print(f"    进度 {done}/{n_batches} 批, 累计 {len(all_rows)} 条记录")
    time.sleep(0.3)

print(f"[3] 组装表达矩阵 ... 共 {len(all_rows)} 条记录")
df = pd.DataFrame(all_rows, columns=["entrez", "patient", "value"])
df = df.dropna()
# 同一基因多患者 -> 透视
mat = df.pivot_table(index="entrez", columns="patient", values="value", aggfunc="mean")
print(f"    矩阵: {mat.shape[0]} 基因 x {mat.shape[1]} 样本")

out = "results/tcga_lihc_fullexpr.csv"
mat.to_csv(out)
print(f"[4] 已保存 {out}")
print("完成。")
