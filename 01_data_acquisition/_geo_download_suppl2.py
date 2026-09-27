import urllib.request, os

OUT = "data/geo"
os.makedirs(OUT, exist_ok=True)

files = [
    ("GSE202069", "GSE202nnn/GSE202069/suppl/GSE202069_gene_tpm_expression.txt.gz"),
    ("GSE149614", "GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz"),
    ("GSE149614", "GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"),
    ("GSE149614", "GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.normalized.txt.gz"),
]

for acc, path in files:
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{path}"
    fname = os.path.basename(path)
    outfile = os.path.join(OUT, f"{acc}_{fname}")
    if os.path.exists(outfile) and os.path.getsize(outfile) > 1000:
        print(f"[{acc}] 已存在，跳过 {fname}", flush=True)
        continue
    print(f"[{acc}] 下载 {fname} ...", flush=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "HCC-download/1.0"})
        with urllib.request.urlopen(req, timeout=1200) as r:
            data = r.read()
        with open(outfile, "wb") as f:
            f.write(data)
        print(f"   成功：{len(data)/1024/1024:.2f} MB -> {outfile}", flush=True)
    except Exception as e:
        print(f"   失败：{e}", flush=True)
print("全部完成", flush=True)
