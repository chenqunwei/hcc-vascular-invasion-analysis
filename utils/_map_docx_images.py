# -*- coding: utf-8 -*-
"""Compare new Fig PNG pixel dims with embedded media dims to build the mapping."""
import struct, zipfile

def png_dims(path):
    with open(path, "rb") as f:
        f.read(8)  # signature
        chunk = f.read(8)
        w, h = struct.unpack(">II", chunk[4:12])
    return w, h

NEW = {
    "Fig2": "results/Fig2_molecular_subtyping.png",
    "Fig3": "results/Fig3_MVI_clinical.png",
    "Fig4": "results/Fig4_survival.png",
    "Fig5": "results/Fig5_immune.png",
}
print("NEW PNGs (w x h, AR=w/h):")
new_dims = {}
for k, p in NEW.items():
    w, h = png_dims(p)
    new_dims[k] = (w, h)
    print(f"  {k}: {w} x {h}  AR={w/h:.3f}")

# Embedded
z = zipfile.ZipFile("HCC血管侵犯连续谱_手稿_v5_fixed.docx")
emb = {}
for n in z.namelist():
    if n.startswith("word/media/"):
        data = z.read(n)
        w, h = struct.unpack(">II", data[16:24])
        emb[n] = (w, h)
print("\nEMBEDDED media (w x h, AR=w/h):")
for n in sorted(emb):
    w, h = emb[n]
    print(f"  {n}: {w} x {h}  AR={w/h:.3f}")
