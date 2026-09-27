import urllib.request, os, sys

OUT = "data/geo"
os.makedirs(OUT, exist_ok=True)

# GEO 数据集 -> Series Matrix FTP 路径（子目录规则：去掉最后3位）
def gse_subdir(acc):
    num = acc[3:]  # 去掉 "GSE"
    prefix = num[:-3]  # 去掉最后3位
    return f"GSE{prefix}nnn"

datasets = [
    "GSE77509",  # PVTT 配对 RNA-seq
    "GSE69164",  # PVTT 配对 RNA-seq
    "GSE14520",  # MVI + 复发
    "GSE202069", # 免疫治疗
    "GSE149614", # 单细胞
]

for acc in datasets:
    sub = gse_subdir(acc)
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{sub}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    outfile = os.path.join(OUT, f"{acc}_series_matrix.txt.gz")
    print(f"[{acc}] 下载 {url}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "HCC-download/1.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        with open(outfile, "wb") as f:
            f.write(data)
        size_mb = len(data) / 1024 / 1024
        print(f"   成功：{size_mb:.2f} MB -> {outfile}")
    except Exception as e:
        print(f"   失败：{e}")
