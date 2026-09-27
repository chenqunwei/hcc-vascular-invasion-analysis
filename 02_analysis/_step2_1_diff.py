import tarfile, gzip, io, os, urllib.request, json, time
import pandas as pd
import numpy as np
from scipy import stats
from collections import Counter

BASE = "data/geo"
os.makedirs("results", exist_ok=True)

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

def load_gse69164(path):
    with gzip.open(path, "rb") as f:
        raw = pd.read_excel(io.BytesIO(f.read()), engine="openpyxl")
    raw = raw.set_index("Gene")
    sample_cols = [c for c in raw.columns if any(k in c for k in ["_ANT", "_HCC", "_PVTT"])]
    return raw[sample_cols].apply(pd.to_numeric, errors="coerce")

def symbol_to_ensembl(symbols):
    """批量 Symbol -> Ensembl ID（mygene API，分批次）"""
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
            q = item.get("query")
            ens = item.get("ensembl", {})
            if isinstance(ens, dict):
                g = ens.get("gene")
            elif isinstance(ens, list) and ens:
                g = ens[0].get("gene")
            else:
                g = None
            if q and g:
                mapping[q] = g
        time.sleep(0.3)
    return mapping

def gradient_analysis(mat, group_key):
    """三组梯度分析，返回 DataFrame（gene, pvalue, up_trend, dn_trend）"""
    logmat = np.log2(mat + 1)
    res = []
    for g in logmat.index:
        vals = logmat.loc[g]
        n_grp = [vals[c] for c in mat.columns if group_key(c) == "N"]
        t_grp = [vals[c] for c in mat.columns if group_key(c) == "T"]
        p_grp = [vals[c] for c in mat.columns if group_key(c) == "P"]
        if len(n_grp) < 2 or len(t_grp) < 2 or len(p_grp) < 2:
            continue
        try:
            h, p = stats.kruskal(n_grp, t_grp, p_grp)
        except Exception:
            continue
        mn, mt, mp = np.median(n_grp), np.median(t_grp), np.median(p_grp)
        up = (mn < mt < mp)
        dn = (mn > mt > mp)
        res.append((g, p, up, dn))
    df = pd.DataFrame(res, columns=["gene", "pvalue", "up_trend", "dn_trend"])
    return df

print("读取数据 ...")
gse77509 = load_gse77509(os.path.join(BASE, "GSE77509_GSE77509_RAW.tar"))
gse69164 = load_gse69164(os.path.join(BASE, "GSE69164_GSE69164_Whole_33_sample_FPKM.xlsx.gz"))
print(f"GSE77509 (Ensembl): {gse77509.shape}, GSE69164 (Symbol): {gse69164.shape}")

# ID 转换：GSE69164 Symbol -> Ensembl
print("\n转换 GSE69164 基因 Symbol -> Ensembl ID ...")
sym2ens = symbol_to_ensembl(gse69164.index)
gse69164_ens = gse69164.rename(index=lambda s: sym2ens.get(s))
# 去掉未映射的行，合并重复 Ensembl ID（取均值）
gse69164_ens = gse69164_ens[gse69164_ens.index.notna()]
gse69164_ens = gse69164_ens[~gse69164_ens.index.duplicated(keep="first")]
print(f"GSE69164 映射后: {gse69164_ens.shape}（未映射 {len(gse69164) - len(gse69164_ens)} 个 Symbol）")

# 梯度分析
g77509 = {c: c[0] for c in gse77509.columns}
g69164 = {c: "P" if "_PVTT" in c else ("N" if "_ANT" in c else "T") for c in gse69164_ens.columns}

r77509 = gradient_analysis(gse77509, lambda c: g77509[c])
r69164 = gradient_analysis(gse69164_ens, lambda c: g69164[c])

for name, r in [("GSE77509", r77509), ("GSE69164", r69164)]:
    sig = r[r["pvalue"] < 0.05]
    print(f"\n{name}: p<0.05 共 {len(sig)}（单调上调 {sig['up_trend'].sum()}，单调下调 {sig['dn_trend'].sum()}）")

# 交集（用 Ensembl ID 对齐）
up77509 = set(r77509[(r77509["pvalue"] < 0.05) & r77509["up_trend"]]["gene"])
dn77509 = set(r77509[(r77509["pvalue"] < 0.05) & r77509["dn_trend"]]["gene"])
up69164 = set(r69164[(r69164["pvalue"] < 0.05) & r69164["up_trend"]]["gene"])
dn69164 = set(r69164[(r69164["pvalue"] < 0.05) & r69164["dn_trend"]]["gene"])

common_up = up77509 & up69164
common_dn = dn77509 & dn69164
print(f"\n两数据集一致「随血管侵犯渐进上调」: {len(common_up)} 个")
print(f"两数据集一致「随血管侵犯渐进下调」: {len(common_dn)} 个")

r77509.to_csv("results/GSE77509_gradient.csv", index=False)
r69164.to_csv("results/GSE69164_gradient.csv", index=False)
pd.Series(sorted(common_up), name="gene").to_csv("results/vascular_invasion_up_genes.csv", index=False)
pd.Series(sorted(common_dn), name="gene").to_csv("results/vascular_invasion_dn_genes.csv", index=False)

print("\n共同上调基因（Ensembl ID，前 30 个）:")
print(sorted(common_up)[:30])
