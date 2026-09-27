import urllib.request, urllib.parse, json
from collections import Counter

# GDC API 查询 TCGA-LIHC 的 diagnoses 字段（含血管侵犯/MVI）
BASE = "https://api.gdc.cancer.gov/cases"

fields = [
    "diagnoses.vascular_invasion_present",
    "diagnoses.vascular_tumor_cell_type",
    "diagnoses.ajcc_pathologic_stage",
    "diagnoses.tumor_stage",
    "diagnoses.days_to_last_follow_up",
]

def query(field_list, size=500):
    filters = {"op": "in", "content": {"field": "project.project_id", "value": ["TCGA-LIHC"]}}
    params = {
        "filters": json.dumps(filters),
        "fields": ",".join(field_list),
        "format": "json",
        "size": str(size),
    }
    url = BASE + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "HCC-tcga/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

res = query(fields)
data = res["data"]["hits"]
print(f"返回病例数：{len(data)}")

# 统计各字段
for f in fields:
    short = f.split(".")[-1]
    cnt = Counter()
    for case in data:
        diag = case.get("diagnoses")
        if diag:
            for d in diag:
                v = d.get(short)
                cnt[str(v)] += 1
        else:
            cnt["<无 diagnoses>"] += 1
    print(f"\n=== {short} ===")
    for k, v in cnt.most_common():
        print(f"  {k}: {v}")
