import urllib.request, os

OUT = "data/geo"
os.makedirs(OUT, exist_ok=True)

# 补充文件（表达矩阵/元数据）下载清单
files = [
    ("GSE77509", "GSE77nnn/GSE77509/suppl/GSE77509_RAW.tar"),
    ("GSE69164", "GSE69nnn/GSE69164/suppl/GSE69164_Whole_33_sample_FPKM.xlsx.gz"),
    ("GSE202069", "GSE202nnn/GSE202069/suppl/GSE202069_gene_tpm_expression.txt.gz"),
    ("GSE149614", "GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz"),
    ("GSE149614", "GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"),
    ("GSE149614", "GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.normalized.txt.gz"),
]

for acc, path in files:
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{path}"
    fname = os.path.basename(path)
    outfile = os.path.join(OUT, f"{acc}_{fname}")
    print(f"[{acc}] 下载 {fname} ...", flush=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "HCC-download/1.0"})
        with urllib.request.urlopen(req, timeout=600) as r:
            data = r.read()
        with open(outfile, "wb") as f:
            f.write(data)
        print(f"   成功：{len(data)/1024/1024:.2f} MB -> {outfile}", flush=True)
    except Exception as e:
        print(f"   失败：{e}", flush=True)
