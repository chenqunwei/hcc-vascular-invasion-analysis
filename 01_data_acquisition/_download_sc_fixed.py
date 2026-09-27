"""重新下载 GSE149614 单细胞 count 矩阵（修正版）。

修正点：
1. 精确目标大小 165,349,783 字节（HEAD 实测），而非模糊的 157.7MB 估算。
2. 严格校验服务器是否接受 Range（返回 206）：若服务器忽略 Range 返回 200，
   则必须从头重写（wb），绝不能 append，否则会重复追加导致文件膨胀损坏。
3. 自动重试 + 断点续传。
"""
import os
import time
import requests

BASE = "data/geo"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
OUT = os.path.join(BASE, "GSE149614_count.txt.gz")
TARGET = 165349783  # 远程文件精确字节数（HEAD Content-Length）
MAX_RETRY = 200
SLEEP = 15


def do_resume():
    resume = os.path.getsize(OUT) if os.path.exists(OUT) else 0
    if resume >= TARGET:
        print(f"[完成] {resume} bytes 已达目标", flush=True)
        return True
    headers = {"User-Agent": "Mozilla/5.0"}
    mode = "wb"
    if resume > 0:
        headers["Range"] = f"bytes={resume}-"
        mode = "ab"
    t0 = time.time()
    with requests.get(URL, timeout=600, stream=True, headers=headers) as r:
        sc = r.status_code
        cl = r.headers.get("Content-Length")
        cr = r.headers.get("Content-Range")
        if resume > 0 and sc != 206:
            # 服务器忽略 Range，返回 200 完整内容 -> 从头重写，避免重复追加
            print(f"[警告] resume={resume} 但服务器返回 {sc}（忽略Range），从头重写", flush=True)
            mode = "wb"
            resume = 0
        else:
            print(f"[续传] status={sc} resume={resume} Content-Range={cr} Content-Length={cl}", flush=True)
        with open(OUT, mode) as f:
            for chunk in r.iter_content(4 * 1024 * 1024):
                if chunk:
                    f.write(chunk)
    sz = os.path.getsize(OUT)
    print(f"[本轮结束] {sz} bytes = {sz / 1048576:.1f} MB (耗时 {time.time() - t0:.0f}s)", flush=True)
    return sz >= TARGET


attempt = 0
while attempt < MAX_RETRY:
    attempt += 1
    try:
        if do_resume():
            print("[成功] count 矩阵下载完成", flush=True)
            raise SystemExit(0)
    except SystemExit:
        raise
    except Exception as e:
        print(f"[重试 {attempt}/{MAX_RETRY}] {type(e).__name__}: {e}", flush=True)
        time.sleep(SLEEP)

print("[失败] 重试次数耗尽，仍未能完成下载", flush=True)
raise SystemExit(1)
