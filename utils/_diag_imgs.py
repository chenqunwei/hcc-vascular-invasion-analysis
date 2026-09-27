# -*- coding: utf-8 -*-
"""Diagnostic: read true PNG pixel dims (new PNGs + embedded media), map each
embedded blip to its figure number via caption scan, and report aspect-ratio
compatibility so we can safely re-embed new Fig1-5."""
import re, zipfile, struct, os

DOCX = "HCC血管侵犯连续谱_手稿_v5_fixed.docx"

def png_dims_from_bytes(data):
    # signature 8B + IHDR: len(4)+type(4)+width(4)+height(4)
    w, h = struct.unpack(">II", data[16:24])
    return w, h

# ---- NEW PNGs on disk ----
NEW = {
    1: "results/Fig1_vascular_program.png",
    2: "results/Fig2_molecular_subtyping.png",
    3: "results/Fig3_MVI_clinical.png",
    4: "results/Fig4_survival.png",
    5: "results/Fig5_immune.png",
}
print("=== NEW PNGs (pixel w x h, AR=w/h) ===")
new_dims = {}
for fig, p in NEW.items():
    with open(p, "rb") as f:
        data = f.read(33)
    w, h = png_dims_from_bytes(data)
    new_dims[fig] = (w, h)
    print(f"  Fig{fig}: {w} x {h}  AR={w/h:.4f}  ({os.path.basename(p)})")

# ---- EMBEDDED media ----
z = zipfile.ZipFile(DOCX)
xml = z.read("word/document.xml").decode("utf-8")
rels = z.read("word/_rels/document.xml.rels").decode("utf-8")
rid_target = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', rels))

# blips in document order
blip_positions = [m.start() for m in re.finditer(r'<a:blip[^>]*r:embed="(rId\d+)"', xml)]
blips = re.findall(r'<a:blip[^>]*r:embed="(rId\d+)"', xml)

# Pre-compute all "Figure N" caption positions (search whole doc)
cap_re = re.compile(r'Figure\s+(\d+)', re.IGNORECASE)

print("\n=== EMBEDDED images in document order ===")
emb_fig = {}  # rId -> figure number
for i, (pos, rid) in enumerate(zip(blip_positions, blips)):
    seg = xml[pos:pos+1400]
    ext = re.search(r'<a:ext cx="(\d+)" cy="(\d+)"', seg)
    cx = int(ext.group(1)) if ext else 0
    cy = int(ext.group(2)) if ext else 0
    tgt = rid_target.get(rid, "?")
    # find closest preceding "Figure N" caption before this blip
    preceding = [(m.start(), int(m.group(1))) for m in cap_re.finditer(xml[:pos])]
    fig_no = preceding[-1][1] if preceding else None
    emb_fig[rid] = fig_no
    # pixel dims of embedded media; tgt is like "media/image1.png"
    full = "word/" + tgt
    try:
        data = z.read(full)
        ew, eh = png_dims_from_bytes(data)
    except Exception as e:
        ew, eh = (-1, -1)
    ar_px = (ew/eh) if eh > 0 else 0.0
    print(f"  blip#{i+1}: {rid} -> {tgt}  display={cx}x{cy} EMU  AR_display={cx/cy:.4f}  px={ew}x{eh} AR_px={ar_px:.4f}  captionFig={fig_no}")

# ---- Compatibility check: for each embedded blip whose caption is Fig 1..5,
# compare its pixel AR to the new PNG AR ----
print("\n=== Re-embed plan & AR compatibility ===")
for i, (pos, rid) in enumerate(zip(blip_positions, blips)):
    fig = emb_fig[rid]
    if fig in new_dims:
        tgt = rid_target.get(rid)
        full = "word/" + tgt
        data = z.read(full)
        ew, eh = png_dims_from_bytes(data)
        nw, nh = new_dims[fig]
        ar_emb = ew/eh
        ar_new = nw/nh
        ratio = ar_new/ar_emb
        flag = "OK" if abs(ratio-1) < 0.02 else ("WARN distort %.3f" % ratio)
        print(f"  {rid} (Fig{fig}): embedded px AR={ar_emb:.4f}  new px AR={ar_new:.4f}  -> {flag}")
    else:
        print(f"  {rid} (Fig{fig}): no new PNG (keep old)")

z.close()
