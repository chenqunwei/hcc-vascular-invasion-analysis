import urllib.request, json, sys
from collections import Counter

def get_clinical(attr_id):
    url = ("https://www.cbioportal.org/api/studies/lihc_tcga/clinical-data"
           f"?clinicalAttributeId={attr_id}&projection=DETAILED")
    req = urllib.request.Request(url, headers={"User-Agent": "HCC-tcga/1.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

for attr in ["VASCULAR_INVASION", "DFS_STATUS", "OS_STATUS"]:
    try:
        data = get_clinical(attr)
        vals = [d.get("value") for d in data]
        cnt = Counter(vals)
        n_total = len(vals)
        n_na = cnt.get("NA", 0) + cnt.get("N/A", 0) + cnt.get("", 0) + cnt.get(None, 0)
        print(f"\n=== {attr}（共 {n_total} 例，非 NA {n_total - n_na} 例）===")
        for k, v in cnt.most_common(20):
            print(f"  {k}: {v}")
    except Exception as e:
        print(f"\n=== {attr} 查询失败：{e} ===", file=sys.stderr)
