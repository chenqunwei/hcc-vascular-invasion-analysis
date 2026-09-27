"""下载剩余数据（用 requests，避开 Git Bash curl SSL 问题）。"""
import os
import time
import requests

BASE = "data/geo"

# 下载清单（基于 NCBI FTP 目录查询结果）
DOWNLOADS = [
    ("GSE36376", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36376/matrix/GSE36376_series_matrix.txt.gz", "GSE36376_series_matrix.txt.gz"),
    ("GSE14520", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_RAW.tar", "GSE14520_GSE14520_RAW.tar"),
    ("GSE149614", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz", "GSE149614_count.txt.gz"),
    ("GSE149614", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.normalized.txt.gz", "GSE149614_normalized.txt.gz"),
    ("GSE242889", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242889/suppl/GSE242889_RAW.tar", "GSE242889_RAW.tar"),
]


def download(url, out_path):
    """requests 流式下载，断点续传支持。"""
    print(f"\n>>> 下载: {out_path}", flush=True)
    t0 = time.time()
    try:
        # 检查已下载字节
        resume = os.path.getsize(out_path) if os.path.exists(out_path) else 0
        headers = {"User-Agent": "Mozilla/5.0"}
        if resume > 0:
            headers["Range"] = f"bytes={resume}-"
            print(f"  断点续传，已下载 {resume} bytes")

        with requests.get(url, timeout=600, stream=True, headers=headers) as r:
            r.raise_for_status()
            total = r.headers.get("Content-Length")
            if total:
                total = int(total) + resume
                print(f"  总大小: {total/1048576:.1f} MB")
            else:
                total = None

            mode = "ab" if resume > 0 else "wb"
            with open(out_path, mode) as f:
                chunk_size = 4 * 1024 * 1024  # 4MB
                downloaded = resume
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total:
                            pct = downloaded / total * 100
                            if int(downloaded / chunk_size) % 10 == 0:
                                speed = downloaded / (time.time() - t0) / 1048576
                                print(f"  {pct:.1f}% ({downloaded/1048576:.1f}/{total/1048576:.1f} MB, {speed:.2f} MB/s)", flush=True)

        sz = os.path.getsize(out_path)
        print(f"  完成: {sz/1048576:.1f} MB, 耗时 {time.time()-t0:.1f} s")
        return sz > 1000
    except Exception as e:
        print(f"  失败: {type(e).__name__}: {e}")
        return False


for acc, url, fname in DOWNLOADS:
    out_path = os.path.join(BASE, fname)
    if os.path.exists(out_path) and os.path.getsize(out_path) > 10000:
        print(f"\n[跳过已存在] {fname}: {os.path.getsize(out_path)/1048576:.1f} MB")
        continue
    ok = download(url, out_path)
    time.sleep(3)  # NCBI 限流保护
    if not ok:
        print(f"  !!! {acc} 下载失败")

print("\n=== 全部下载完成 ===")
print("当前 data/geo/ 文件:")
for f in sorted(os.listdir(BASE)):
    sz = os.path.getsize(os.path.join(BASE, f))
    print(f"  {f}: {sz/1048576:.1f} MB")