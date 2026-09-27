import json
import time
import urllib.request
import requests
import pandas as pd
import numpy as np
from scipy import stats

print("=" * 60)
print("Step 2.2 续跑：拉取剩余基因 + 组装矩阵 + 评分")
print("=" * 60)

PROFILE = "lihc_tcga_rna_seq_v2_mrna"
SAMPLE_LIST = "lihc_tcga_rna_seq_v2_mrna"
CBIO_URL = f"https://www.cbioportal.org/api/molecular-profiles/{PROFILE}/molecular-data/fetch"
RAW_PATH = "results/tcga_lihc_vi_raw_records.csv"

# ---------- 1. 加载基因 + Ensembl->Entrez ----------
up = pd.read_csv("results/vascular_invasion_up_genes.csv")["gene"].tolist()
dn = pd.read_csv("results/vascular_invasion_dn_genes.csv")["gene"].tolist()

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
            q = it.get("query"); ent = it.get("entrezgene")
            if q and ent:
                mapping[q] = int(ent)
        time.sleep(0.3)
    return mapping

print("[1] Ensembl -> Entrez ...")
ens2ent = ensembl_to_entrez(up + dn)
up_ent = [ens2ent[g] for g in up if g in ens2ent]
dn_ent = [ens2ent[g] for g in dn if g in ens2ent]
all_ent = list(set(up_ent + dn_ent))
print(f"    映射 {len(ens2ent)} 基因，up={len(up_ent)} dn={len(dn_ent)}")

# ---------- 2. 找剩余基因 ----------
print("[2] 检查已落盘数据，找剩余基因 ...")
done = set()
try:
    raw = pd.read_csv(RAW_PATH, usecols=["entrez"])
    done = set(raw["entrez"].unique())
except Exception:
    pass
remaining = [e for e in all_ent if e not in done]
print(f"    已覆盖 {len(done)}，剩余 {len(remaining)}")

# ---------- 3. fetch 剩余 ----------
def cbio_fetch(entrez_ids, retry=5):
    body = {"sampleListId": SAMPLE_LIST, "entrezGeneIds": entrez_ids}
    for attempt in range(retry):
        try:
            r = requests.post(CBIO_URL, json=body, timeout=120)
            if r.status_code == 200:
                data = r.json()
                if isinstance(data, list):
                    return [(d.get("entrezGeneId"), d.get("patientId"), d.get("value")) for d in data]
                print(f"    非列表: {json.dumps(data)[:150]}")
                return []
            print(f"    HTTP {r.status_code} (尝试{attempt+1}): {r.text[:80]}")
            time.sleep(4)
        except Exception as e:
            print(f"    异常 (尝试{attempt+1}): {type(e).__name__}: {e}")
            time.sleep(4)
    return []

if remaining:
    print(f"[3] 拉取剩余 {len(remaining)} 基因 ...")
    for i in range(0, len(remaining), 150):
        batch = remaining[i:i+150]
        rec = cbio_fetch(batch)
        if rec:
            pd.DataFrame(rec, columns=["entrez", "patient", "value"]).to_csv(
                RAW_PATH, mode="a", header=False, index=False)
        print(f"    剩余批次 {i//150+1}: 基因 {len(batch)}，记录 {len(rec)}")
        time.sleep(0.3)
else:
    print("[3] 无剩余基因，跳过 fetch")

# ---------- 4. 组装矩阵 ----------
print("[4] 组装矩阵 ...")
raw = pd.read_csv(RAW_PATH)
raw["value"] = pd.to_numeric(raw["value"], errors="coerce")
mat = raw.pivot_table(index="entrez", columns="patient", values="value", aggfunc="mean")
print(f"    矩阵: {mat.shape[0]} 基因 x {mat.shape[1]} 患者")

# ---------- 5. 评分 ----------
print("[5] 计算 VI score ...")
logmat = np.log2(mat + 1)
z = logmat.sub(logmat.mean(axis=1), axis=0).div(logmat.std(axis=1), axis=0)
z = z.replace([np.inf, -np.inf], np.nan).fillna(0)
up_in = [e for e in up_ent if e in z.index]
dn_in = [e for e in dn_ent if e in z.index]
print(f"    上调命中 {len(up_in)}/{len(up_ent)}，下调命中 {len(dn_in)}/{len(dn_ent)}")
vi_score = z.loc[up_in].mean(axis=0) - z.loc[dn_in].mean(axis=0)

# ---------- 6. 关联 MVI + AUC ----------
print("[6] 关联 MVI ...")
clin = pd.read_csv("results/tcga_lihc_clinical.csv")
clin = clin[clin["mvi"].isin([0, 1])]
merged = pd.DataFrame({"vi_score": vi_score}).join(
    clin.set_index("patient_id")[["mvi", "recurrence", "death", "dfs_months", "os_months"]], how="inner")
merged = merged.dropna(subset=["mvi"])
print(f"    有 VI 分数 + MVI 标注: {len(merged)} 例")
mvi_pos = merged[merged["mvi"] == 1]["vi_score"]
mvi_neg = merged[merged["mvi"] == 0]["vi_score"]
u, p = stats.mannwhitneyu(mvi_pos, mvi_neg, alternative="two-sided")
auc = u / (len(mvi_pos) * len(mvi_neg))
print(f"    MVI+ (n={len(mvi_pos)}): 均值 {mvi_pos.mean():.3f}")
print(f"    MVI- (n={len(mvi_neg)}): 均值 {mvi_neg.mean():.3f}")
print(f"    Mann-Whitney p = {p:.4e}")
print(f"    AUC = {auc:.3f}")

# ---------- 7. 保存 ----------
out = merged.reset_index().rename(columns={"index": "patient_id"})
out = out[["patient_id", "vi_score", "mvi", "recurrence", "death", "dfs_months", "os_months"]]
out.to_csv("results/tcga_lihc_vi_score.csv", index=False)
mat.to_csv("results/tcga_lihc_vi_genes_expr.csv")
print(f"\n[7] 已保存 results/tcga_lihc_vi_score.csv（{len(out)} 例）")
print("\n评分分布摘要:")
print(out["vi_score"].describe().round(3).to_string())
print("\nStep 2.2 完成。")
