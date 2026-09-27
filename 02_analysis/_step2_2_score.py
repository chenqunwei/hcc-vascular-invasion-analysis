import json
import time
import urllib.request
import requests
import pandas as pd
import numpy as np
from scipy import stats

print("=" * 60)
print("Step 2.2  投影 TCGA-LIHC（cBioPortal RSEM）算血管侵犯倾向分数")
print("=" * 60)

PROFILE = "lihc_tcga_rna_seq_v2_mrna"
SAMPLE_LIST = "lihc_tcga_rna_seq_v2_mrna"
CBIO_URL = f"https://www.cbioportal.org/api/molecular-profiles/{PROFILE}/molecular-data/fetch"

# ---------- 1. 加载血管侵犯基因 ----------
up = pd.read_csv("results/vascular_invasion_up_genes.csv")["gene"].tolist()
dn = pd.read_csv("results/vascular_invasion_dn_genes.csv")["gene"].tolist()
print(f"\n[1] 血管侵犯基因: 上调 {len(up)} + 下调 {len(dn)} = {len(up)+len(dn)}")

# ---------- 2. Ensembl -> Entrez（mygene） ----------
def ensembl_to_entrez(ids, batch=1000):
    mapping = {}
    for i in range(0, len(ids), batch):
        b = ids[i:i+batch]
        body = json.dumps({"q": ",".join(b), "scopes": "ensembl.gene",
                           "fields": "entrezgene", "species": "human"}).encode()
        req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            res = json.load(r)
        for it in res:
            q = it.get("query")
            ent = it.get("entrezgene")
            if q and ent:
                mapping[q] = int(ent)
        time.sleep(0.3)
    return mapping

print("\n[2] Ensembl -> Entrez Gene ID（mygene）...")
all_genes = up + dn
ens2ent = ensembl_to_entrez(all_genes)
print(f"    成功映射 {len(ens2ent)}/{len(all_genes)} 个基因")

up_ent = [ens2ent[g] for g in up if g in ens2ent]
dn_ent = [ens2ent[g] for g in dn if g in ens2ent]

# ---------- 3. cBioPortal 批量 fetch ----------
def cbio_fetch(entrez_ids, retry=4):
    body = {"sampleListId": SAMPLE_LIST, "entrezGeneIds": entrez_ids}
    for attempt in range(retry):
        try:
            r = requests.post(CBIO_URL, json=body, timeout=120)
            if r.status_code == 200:
                data = r.json()
                if isinstance(data, list):
                    return [(d.get("entrezGeneId"), d.get("patientId"), d.get("value")) for d in data]
                print(f"    fetch 非列表: {json.dumps(data)[:200]}")
                return []
            else:
                print(f"    HTTP {r.status_code} (尝试{attempt+1}/{retry}): {r.text[:100]}")
                time.sleep(3)
        except Exception as e:
            print(f"    请求异常 (尝试{attempt+1}/{retry}): {type(e).__name__}: {e}")
            time.sleep(3)
    return []

print("\n[3] 从 cBioPortal 批量拉取表达（每批 150 基因，增量落盘）...")
all_ent = list(set(up_ent + dn_ent))
n_batch = (len(all_ent) + 149) // 150
raw_path = "results/tcga_lihc_vi_raw_records.csv"
# 清空旧文件
with open(raw_path, "w") as f:
    f.write("entrez,patient,value\n")
for i in range(0, len(all_ent), 150):
    batch = all_ent[i:i+150]
    rec = cbio_fetch(batch)
    if rec:
        df_batch = pd.DataFrame(rec, columns=["entrez", "patient", "value"])
        df_batch.to_csv(raw_path, mode="a", header=False, index=False)
    print(f"    批次 {i//150+1}/{n_batch}: 基因 {len(batch)}，记录 {len(rec)}")
    time.sleep(0.2)

print("\n[4] 读取记录并组装基因 x 患者矩阵 ...")
raw = pd.read_csv(raw_path)
raw["value"] = pd.to_numeric(raw["value"], errors="coerce")
print(f"    总记录数: {len(raw)}")
mat = raw.pivot_table(index="entrez", columns="patient", values="value", aggfunc="mean")
print(f"    矩阵: {mat.shape[0]} 基因 x {mat.shape[1]} 患者")

# ---------- 5. 算血管侵犯评分 ----------
print("\n[5] 计算血管侵犯倾向分数（VI score）...")
logmat = np.log2(mat + 1)
z = logmat.sub(logmat.mean(axis=1), axis=0).div(logmat.std(axis=1), axis=0)
z = z.replace([np.inf, -np.inf], np.nan).fillna(0)

up_in = [e for e in up_ent if e in z.index]
dn_in = [e for e in dn_ent if e in z.index]
print(f"    上调命中 {len(up_in)}/{len(up_ent)}，下调命中 {len(dn_in)}/{len(dn_ent)}")
up_score = z.loc[up_in].mean(axis=0)
dn_score = z.loc[dn_in].mean(axis=0)
vi_score = up_score - dn_score

# ---------- 6. 关联 MVI + 区分度 ----------
print("\n[6] 关联临床 MVI 状态，验证区分度 ...")
clin = pd.read_csv("results/tcga_lihc_clinical.csv")
clin = clin[clin["mvi"].isin([0, 1])]
merged = pd.DataFrame({"vi_score": vi_score}).join(
    clin.set_index("patient_id")[["mvi", "recurrence", "death", "dfs_months", "os_months"]],
    how="inner")
merged = merged.dropna(subset=["mvi"])
print(f"    有 VI 分数 + MVI 标注: {len(merged)} 例")

mvi_pos = merged[merged["mvi"] == 1]["vi_score"]
mvi_neg = merged[merged["mvi"] == 0]["vi_score"]
u, p = stats.mannwhitneyu(mvi_pos, mvi_neg, alternative="two-sided")
auc = u / (len(mvi_pos) * len(mvi_neg))
print(f"    MVI+ (Micro, n={len(mvi_pos)}): 均值 {mvi_pos.mean():.3f}")
print(f"    MVI- (None, n={len(mvi_neg)}): 均值 {mvi_neg.mean():.3f}")
print(f"    Mann-Whitney p = {p:.4e}")
print(f"    AUC = {auc:.3f}")

# ---------- 7. 保存 ----------
out = merged.reset_index().rename(columns={"index": "patient_id"})
out = out[["patient_id", "vi_score", "mvi", "recurrence", "death", "dfs_months", "os_months"]]
out.to_csv("results/tcga_lihc_vi_score.csv", index=False)
mat.to_csv("results/tcga_lihc_vi_genes_expr.csv")
print(f"\n[7] 已保存:")
print(f"    results/tcga_lihc_vi_score.csv（{len(out)} 例 VI 分数）")
print(f"    results/tcga_lihc_vi_genes_expr.csv（{mat.shape[0]} 基因 x {mat.shape[1]} 样本表达）")

print("\n评分分布摘要:")
print(out["vi_score"].describe().round(3).to_string())
print("\nStep 2.2 完成。")
