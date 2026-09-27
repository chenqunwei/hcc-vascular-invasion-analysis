# -*- coding: utf-8 -*-
"""Inspect images embedded in v5_fixed.docx: order of blips, rId->target, and EMU size."""
import re, zipfile

DOCX = "HCC血管侵犯连续谱_手稿_v5_fixed.docx"
z = zipfile.ZipFile(DOCX)
names = z.namelist()
# rels
rels = z.read("word/_rels/document.xml.rels").decode("utf-8")
rid_target = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', rels))
# document.xml
xml = z.read("word/document.xml").decode("utf-8")

# blips in order
blips = re.findall(r'<a:blip[^>]*r:embed="(rId\d+)"', xml)
print("Number of blips:", len(blips))
# For each blip, extract the enclosing <a:ext cx=".." cy="..">
# find ext near each blip by scanning after the blip index
positions = [m.start() for m in re.finditer(r'<a:blip[^>]*r:embed="(rId\d+)"', xml)]
for i, (pos, rid) in enumerate(zip(positions, blips)):
    seg = xml[pos:pos+1200]
    ext = re.search(r'<a:ext cx="(\d+)" cy="(\d+)"', seg)
    tgt = rid_target.get(rid, "?")
    cx = int(ext.group(1)) if ext else 0
    cy = int(ext.group(2)) if ext else 0
    ar = round(cx/cy, 3) if cy else 0
    print(f"  blip#{i+1}: {rid} -> {tgt}  {cx}x{cy} EMU  AR={ar}")
