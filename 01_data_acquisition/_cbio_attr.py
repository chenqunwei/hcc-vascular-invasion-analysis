import urllib.request, json

# 查 cBioPortal 的 TCGA-LIHC (lihc_tcga) 临床属性
url = "https://www.cbioportal.org/api/studies/lihc_tcga/clinical-attributes"
req = urllib.request.Request(url, headers={"User-Agent": "HCC-tcga/1.0", "Accept": "application/json"})
with urllib.request.urlopen(req, timeout=60) as r:
    attrs = json.load(r)

print(f"临床属性总数：{len(attrs)}")
print("\n包含 vascular/invasion/micro 的字段：")
for a in attrs:
    name = a.get("displayName", "")
    desc = a.get("description", "")
    if any(k in name.lower() for k in ["vascular", "invasion", "micro", "mvi"]):
        print(f"  - {a.get('clinicalAttributeId')} | {name} | {desc[:80]}")

print("\n全部字段名（前 60 个）：")
for a in attrs[:60]:
    print(f"  - {a.get('clinicalAttributeId')} | {a.get('displayName')}")
