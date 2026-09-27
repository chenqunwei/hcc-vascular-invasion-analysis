import urllib.request, urllib.parse, json, time

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

queries = {
    "A PVTT+纯生信(生信/预后模型/签名/分型)": '("portal vein tumor thrombus"[tiab] OR "portal vein tumour thrombus"[tiab]) AND hepatocellular carcinoma[tiab] AND (bioinformatics[tiab] OR "prognostic model"[tiab] OR "gene signature"[tiab] OR "risk score"[tiab] OR "molecular subtype"[tiab] OR "consensus clustering"[tiab])',
    "B PVTT+单细胞/空间组学": '("portal vein tumor thrombus"[tiab] OR "portal vein tumour thrombus"[tiab]) AND hepatocellular carcinoma[tiab] AND ("single-cell"[tiab] OR scRNA[tiab] OR "spatial transcriptom"[tiab] OR "single cell"[tiab])',
    "C 血管侵犯+HCC+生信(泛)": '"vascular invasion"[tiab] AND hepatocellular carcinoma[tiab] AND (bioinformatics[tiab] OR "prognostic signature"[tiab] OR "gene signature"[tiab])',
    "D 血管侵犯连续谱 MVI+PVTT+生信(关键)": '("microvascular invasion"[tiab] OR MVI[tiab]) AND ("portal vein tumor thrombus"[tiab] OR PVTT[tiab]) AND hepatocellular carcinoma[tiab] AND (bioinformatics[tiab] OR signature[tiab] OR "gene set"[tiab] OR transcriptomic[tiab])',
}

def esearch(term, retmax=8):
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
    summ = esummary(ids)
    for pid in ids:
        rec = summ.get(pid, {})
        title = rec.get("title", "")
        year = (rec.get("pubdate", "") or "")[:4]
        jr = rec.get("fulljournalname", "")
        print(f"  - PMID {pid} ({year}) {jr}: {title[:95]}")
    time.sleep(0.4)
