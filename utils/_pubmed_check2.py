import urllib.request, urllib.parse, json, time

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

queries = {
    "Q1 血管侵犯/PVTT+基因签名+免疫治疗响应预测": '("vascular invasion"[tiab] OR "portal vein tumor thrombus"[tiab] OR "portal vein tumour thrombus"[tiab] OR PVTT[tiab]) AND hepatocellular carcinoma[tiab] AND (signature[tiab] OR "gene set"[tiab] OR "risk score"[tiab] OR "prognostic model"[tiab] OR biomarker[tiab]) AND (immunotherapy[tiab] OR "immune checkpoint"[tiab] OR "anti-PD"[tiab] OR "PD-1"[tiab] OR atezolizumab[tiab] OR nivolumab[tiab]) AND (response[tiab] OR predict[tiab] OR sensitivity[tiab] OR efficacy[tiab])',
    "Q2 血管侵犯+基因签名+系统治疗/药物敏感性": '("vascular invasion"[tiab] OR "portal vein tumor thrombus"[tiab] OR PVTT[tiab]) AND hepatocellular carcinoma[tiab] AND ("gene signature"[tiab] OR "gene set"[tiab] OR "prognostic signature"[tiab]) AND (treatment[tiab] OR therapy[tiab] OR "drug sensitivity"[tiab] OR sorafenib[tiab] OR lenvatinib[tiab])',
    "Q3 晚期HCC+血管侵犯+免疫治疗预测(宽)": '"advanced hepatocellular carcinoma"[tiab] AND ("vascular invasion"[tiab] OR "portal vein tumor thrombus"[tiab] OR PVTT[tiab]) AND (immunotherapy[tiab] OR "immune checkpoint inhibitor"[tiab] OR "systemic therapy"[tiab]) AND (biomarker[tiab] OR signature[tiab] OR predict[tiab])',
    "Q4 关键 PVTT/血管侵犯+分子签名+免疫/系统治疗": '("portal vein tumor thrombus"[tiab] OR "vascular invasion"[tiab]) AND hepatocellular carcinoma[tiab] AND ("gene set"[tiab] OR "gene signature"[tiab] OR "molecular signature"[tiab] OR "transcriptomic signature"[tiab]) AND ("immune checkpoint"[tiab] OR immunotherapy[tiab] OR "systemic treatment"[tiab])',
    "Q5 PVTT+免疫微环境/免疫治疗(竞争基线)": '("portal vein tumor thrombus"[tiab] OR PVTT[tiab]) AND hepatocellular carcinoma[tiab] AND (immunotherapy[tiab] OR "immune microenvironment"[tiab] OR "immune checkpoint"[tiab])',
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
