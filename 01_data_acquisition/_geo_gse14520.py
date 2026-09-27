import urllib.request, os

OUT = "data/geo"
os.makedirs(OUT, exist_ok=True)

# GSE14520 supplementary 文件（用 GEO download 接口）
files = [
    "GSE14520_Extra_Supplement.txt.gz",
    "GSE14520_README.txt",
]

for fname in files:
    url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE14520&format=file&file={fname}"
    outfile = os.path.join(OUT, fname)
    print(f"下载 {fname} ...", flush=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "HCC-download/1.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        with open(outfile, "wb") as f:
            f.write(data)
        print(f"   成功：{len(data)/1024:.1f} KB", flush=True)
    except Exception as e:
        print(f"   失败：{e}", flush=True)

# GSE14520 Series Matrix（表达数据）
url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520_series_matrix.txt.gz"
outfile = os.path.join(OUT, "GSE14520_series_matrix.txt.gz")
print("下载 GSE14520 Series Matrix ...", flush=True)
try:
    req = urllib.request.Request(url, headers={"User-Agent": "HCC-download/1.0"})
    with urllib.request.urlopen(req, timeout=300) as r:
        data = r.read()
    with open(outfile, "wb") as f:
        f.write(data)
    print(f"   成功：{len(data)/1024/1024:.2f} MB", flush=True)
except Exception as e:
    print(f"   失败：{e}", flush=True)
