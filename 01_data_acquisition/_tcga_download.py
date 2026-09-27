import urllib.request, os

OUT = "data/tcga"
os.makedirs(OUT, exist_ok=True)

files = {
    "TCGA-LIHC.star_tpm.tsv.gz": "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.star_tpm.tsv.gz",
    "TCGA-LIHC.GDC_phenotype.tsv.gz": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.GDC_phenotype.tsv.gz",
    "TCGA-LIHC.survival.tsv": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.survival.tsv",
}

for fname, url in files.items():
    outfile = os.path.join(OUT, fname)
    if os.path.exists(outfile) and os.path.getsize(outfile) > 1000:
        print(f"[跳过] 已存在 {fname}", flush=True)
        continue
    print(f"[下载] {fname} ...", flush=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "HCC-download/1.0"})
        with urllib.request.urlopen(req, timeout=1200) as r:
            data = r.read()
        with open(outfile, "wb") as f:
            f.write(data)
        print(f"   成功：{len(data)/1024/1024:.2f} MB -> {outfile}", flush=True)
    except Exception as e:
        print(f"   失败：{e}", flush=True)
print("完成", flush=True)
