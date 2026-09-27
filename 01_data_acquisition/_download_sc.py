"""下载 GSE149614 单细胞矩阵（Phase 3 优先）+ 补完其他。"""
import os
import time
import requests

BASE = "data/geo"

# GSE149614 count 优先（Phase 3），其余补完（跳过 normalized，太大且不需要）
DOWNLOADS = [
    ("GSE149614", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz", "GSE149614_count.txt.gz"),
    ("GSE36376", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36376/matrix/GSE36376_series_matrix.txt.gz", "GSE36376_series_matrix.txt.gz"),
    ("GSE14520", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_RAW.tar", "GSE14520_GSE14520_RAW.tar"),
    ("GSE242889", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242889/suppl/GSE242889_RAW.tar", "GSE242889_RAW.tar"),
]


def download(url, out_path):
    print(f"\n>>> {out_path}", flush=True)
    t0 = time.time()
    resume = os.path.getsize(out_path) if os.path.exists(out_path) else 0
    headers = {"User-Agent": "Mozilla/5.0"}
    if resume > 0:
        headers["Range"] = f"bytes={resume}-"
        print(f"  断点续传 {resume} bytes")
    try:
        with requests.get(url, timeout=600, stream=True, headers=headers) as r:
            r.raise_for_status()
            total = r.headers.get("Content-Length")
            total = (int(total) + resume) if total else None
            mode = "ab" if resume > 0 else "wb"
            downloaded = resume
            last_pct = -1
            with open(out_path, mode) as f:
                for chunk in r.iter_content(4 * 1024 * 1024):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total:
                            pct = int(downloaded / total * 100)
                            if pct >= last_pct + 10:
                                print(f"  {pct}% ({downloaded/1048576:.1f}/{total/1048576:.1f} MB)", flush=True)
                                last_pct = pct
        sz = os.path.getsize(out_path)
        print(f"  完成 {sz/1048576:.1f} MB ({time.time()-t0:.0f}s)")
        return sz > 1000
    except Exception as e:
        print(f"  失败: {type(e).__name__}: {e}")
        return False


for acc, url, fname in DOWNLOADS:
    out_path = os.path.join(BASE, fname)
    # 判断是否已完整（GSE36376 完整约 102MB）
    if fname == "GSE36376_series_matrix.txt.gz" and os.path.exists(out_path) and os.path.getsize(out_path) > 100 * 1048576:
        print(f"\n[跳过] {fname} 已完整")
        continue
    if os.path.exists(out_path) and os.path.getsize(out_path) > 100000 and fname != "GSE36376_series_matrix.txt.gz":
        print(f"\n[跳过] {fname} 已存在 {os.path.getsize(out_path)/1048576:.1f} MB")
        continue
    download(url, out_path)
    time.sleep(3)

print("\n=== 下载结束，data/geo/ 现状 ===")
for f in sorted(os.listdir(BASE)):
    print(f"  {f}: {os.path.getsize(os.path.join(BASE,f))/1048576:.1f} MB")
