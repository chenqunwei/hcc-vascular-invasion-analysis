# -*- coding: utf-8 -*-
"""Step 4.2 前置：下载 GSE14520 series matrix + GPL 注释（断点续传版）"""
import requests, os, time

BASE = "https://ftp.ncbi.nlm.nih.gov/geo"
FILES = [
    (f"{BASE}/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL3921_series_matrix.txt.gz", "data/geo/GSE14520-GPL3921_series_matrix.txt.gz"),
    (f"{BASE}/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL571_series_matrix.txt.gz", "data/geo/GSE14520-GPL571_series_matrix.txt.gz"),
    (f"{BASE}/platforms/GPL3nnn/GPL3921/annot/GPL3921.annot.gz", "data/geo/GPL3921.annot.gz"),
    (f"{BASE}/platforms/GPLnnn/GPL571/annot/GPL571.annot.gz", "data/geo/GPL571.annot.gz"),
]
# 预期大小（MB，用于判断完整性；GPL 注释用 0 表示不校验大小）
EXPECT_MB = {
    "GSE14520-GPL3921_series_matrix.txt.gz": 20.3,
    "GSE14520-GPL571_series_matrix.txt.gz": 2.1,
}

os.makedirs("data/geo", exist_ok=True)

def is_done(out):
    if not os.path.exists(out):
        return False
    sz = os.path.getsize(out)
    if sz < 1000:
        return False
    name = os.path.basename(out)
    if name in EXPECT_MB and sz < EXPECT_MB[name] * 0.9 * 1048576:
        return False
    return True

def download(url, out):
    resume = os.path.getsize(out) if os.path.exists(out) else 0
    headers = {"User-Agent": "Mozilla/5.0"}
    if resume > 0:
        headers["Range"] = f"bytes={resume}-"
    r = requests.get(url, headers=headers, stream=True, timeout=60)
    if r.status_code == 200 or r.status_code == 206:
        mode = "ab" if resume > 0 else "wb"
        with open(out, mode) as f:
            for chunk in r.iter_content(256 * 1024):
                f.write(chunk)
        return True
    return False

for url, out in FILES:
    if is_done(out):
        print(f"已完整(跳过): {os.path.basename(out)} ({os.path.getsize(out)/1048576:.1f} MB)")
        continue
    print(f"下载: {os.path.basename(out)} ...")
    for attempt in range(3):
        try:
            if download(url, out):
                if is_done(out):
                    print(f"  完成 ({os.path.getsize(out)/1048576:.1f} MB)")
                    break
                else:
                    print(f"  未完整，重试 ...")
            else:
                print(f"  HTTP 异常，重试 ...")
        except Exception as e:
            print(f"  失败 {type(e).__name__}: {e}, 重试 ...")
        time.sleep(3)

print("\n=== 下载结果 ===")
for _, out in FILES:
    if os.path.exists(out):
        sz = os.path.getsize(out) / 1048576
        done = "✅" if is_done(out) else "⚠️未完整"
        print(f"  {os.path.basename(out)}: {sz:.1f} MB {done}")
    else:
        print(f"  {os.path.basename(out)}: 缺失")
