import urllib.request, urllib.parse, json, time

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

queries = {
    "Q1 血管侵犯+HCC+分型/签名/预后": '(vascular invasion[tiab] OR "microvascular invasion"[tiab]) AND hepatocellular carcinoma[tiab] AND (subtype[tiab] OR classification[tiab] OR signature[tiab] OR prognostic[tiab])',
    "Q2 MVI+分子分型/预后签名": '"microvascular invasion"[tiab] AND hepatocellular carcinoma[tiab] AND (molecular subtype[tiab] OR "molecular classification"[tiab] OR "prognostic signature"[tiab])',
    "Q3 PVTT+分型/签名/预后模型": '("portal vein tumor thrombus"[tiab] OR "portal vein tumour thrombus"[tiab]) AND hepatocellular carcinoma[tiab] AND (signature[tiab] OR classification[tiab] OR subtype[tiab] OR "prognostic model"[tiab])',
    "Q4 血管侵犯+免疫/系统治疗响应": '(vascular invasion[tiab] OR "microvascular invasion"[tiab]) AND hepatocellular carcinoma[tiab] AND (immunotherapy[tiab] OR "immune checkpoint"[tiab] OR "systemic therapy"[tiab]) AND (response[tiab] OR predict[tiab])',
    "Q5 关键查重 MVI+PVTT整合(连续谱)": '("microvascular invasion"[tiab] OR MVI[tiab]) AND ("portal vein tumor thrombus"[tiab] OR PVTT[tiab]) AND hepatocellular carcinoma[tiab]',
    "Q6 血管侵犯+单细胞/空间组学": '(vascular invasion[tiab] OR "microvascular invasion"[tiab]) AND hepatocellular carcinoma[tiab] AND ("single-cell"[tiab] OR "spatial transcriptom"[tiab] OR scRNA[tiab])',
}

def esearch(term, retmax=12):
    url = BASE + "esearch.fcgi?" + urllib.parse.urlencode({
        "db": "pubmed", "term": term, "retmax": retmax,
        "retmode": "json", "sort": "relevance"})
    req = urllib.request.Request(url, headers={"User-Agent": "HCC-research-check/1.0"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)["esearchresult"]

def esummary(ids):
    if not ids:
        return {}
    url = BASE + "esummary.fcgi?" + urllib.parse.urlencode({
        "db": "pubmed", "id": ",".join(ids), "retmode": "json"})
    req = urllib.request.Request(url, headers={"User-Agent": "HCC-research-check/1.0"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)["result"]

for name, term in queries.items():
    try:
        res = esearch(term)
    except Exception as e:
        print("=" * 64)
        print(f"[{name}]  检索失败: {e}")
        continue
    count = res.get("count", "?")
    ids = res.get("idlist", [])
    print("=" * 64)
    print(f"[{name}]  命中 {count} 篇")
    print("检索式: " + term)
    summ = esummary(ids)
    for pid in ids:
        rec = summ.get(pid, {})
        title = rec.get("title", "")
        year = (rec.get("pubdate", "") or "")[:4]
        jr = rec.get("fulljournalname", "")
        print(f"  - PMID {pid} ({year}) {jr}: {title[:100]}")
    time.sleep(0.4)
