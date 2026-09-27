import urllib.request, urllib.parse, json

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# 需要下载的 GEO 数据集
datasets = {
    "GSE77509": "PVTT 配对 RNA-seq（19 例正常/原发/PVTT）",
    "GSE69164": "PVTT 配对 RNA-seq（11 例）",
    "GSE14520": "MVI + 复发（225 例，芯片）",
    "GSE202069": "免疫治疗（17 例 anti-PD1 有响应）",
    "GSE149614": "单细胞（10 肿瘤 + 8 癌旁）",
}

def esummary_gds(accession):
    url = BASE + "esummary.fcgi?" + urllib.parse.urlencode({
        "db": "gds", "id": accession, "retmode": "json"})
    req = urllib.request.Request(url, headers={"User-Agent": "HCC-download/1.0"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)["result"]

for acc, desc in datasets.items():
    print("=" * 70)
    print(f"[{acc}] {desc}")
    try:
        res = esummary_gds(acc)
        uid = list(res.keys())[0]
        rec = res[uid]
        print("  标题:", rec.get("title", ""))
        print("  样本数:", rec.get("n_samples", "?"))
        print("  平台:", (rec.get("gpl", "") or "")[:60])
        ftp = rec.get("ftplink", "")
        if ftp:
            print("  FTP:", ftp)
        supp = rec.get("suppfile", "")
        if supp:
            print("  补充文件:", supp[:120])
        print("  状态:", rec.get("status", ""))
    except Exception as e:
        print("  获取失败:", e)
