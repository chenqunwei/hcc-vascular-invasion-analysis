"""仅续传 GSE149614 单细胞 count 矩阵（断点续传 + 自动重试），不触碰其他文件。"""
import os
import time
import requests

BASE = "data/geo"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
OUT = os.path.join(BASE, "GSE149614_count.txt.gz")
TARGET = 157.7 * 1048576  # 目标完整大小约 157.7 MB
MAX_RETRY = 50  # 每轮失败最多重试次数


def do_resume():
    resume = os.path.getsize(OUT) if os.path.exists(OUT) else 0
    if resume >= TARGET - 1024 * 1024:
        print(f"[完成] 文件已达目标大小 {resume/1048576:.1f} MB", flush=True)
        return True
    headers = {"User-Agent": "Mozilla/5.0"}
    if resume > 0:
        headers["Range"] = f"bytes={resume}-"
        print(f"[续传] 从 {resume/1048576:.1f} MB 继续", flush=True)
    t0 = time.time()
    last_pct = -1
    with requests.get(URL, timeout=600, stream=True, headers=headers) as r:
        r.raise_for_status()
        total = r.headers.get("Content-Length")
        total = (int(total) + resume) if total else None
        downloaded = resume
        with open(OUT, "ab" if resume > 0 else "wb") as f:
            for chunk in r.iter_content(4 * 1024 * 1024):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = int(downloaded / total * 100)
                        if pct >= last_pct + 10:
                            print(f"  {pct}% ({downloaded/1048576:.1f}/{total/1048576:.1f} MB)", flush=True)
                            last_pct = pct
    sz = os.path.getsize(OUT)
    print(f"[完成] {sz/1048576:.1f} MB ({time.time()-t0:.0f}s)", flush=True)
    return sz >= TARGET - 1024 * 1024


attempt = 0
while attempt < MAX_RETRY:
    attempt += 1
    try:
        if do_resume():
            raise SystemExit(0)
    except SystemExit:
        raise
    except Exception as e:
        print(f"[重试 {attempt}/{MAX_RETRY}] 失败: {type(e).__name__}: {e}", flush=True)
        time.sleep(10)
print("[失败] 重试次数耗尽，仍未能完成下载", flush=True)
raise SystemExit(1)
