import re
from pathlib import Path
from docx import Document

SRC = Path(r"E:/workbuddy/HCC pure bioinformatics analysis/HCC血管侵犯连续谱_手稿_v5.docx")
doc = Document(str(SRC))
paras = doc.paragraphs

ref_start = None
for i,p in enumerate(paras):
    if re.match(r'^(references|参考文献)\s*$', p.text.strip(), re.I):
        ref_start = i; break

cit_re = re.compile(r'\[\d+\]')
print("=== BODY PARAGRAPHS CONTAINING CITATIONS ===")
for i,p in enumerate(paras):
    if ref_start is not None and i >= ref_start:
        continue
    t = p.text.strip()
    if cit_re.search(t):
        # mark the heading-ish style
        tag = ""
        if p.style and p.style.name:
            tag = f"[{p.style.name}]"
        print(f"\n--- P{i} {tag} ---")
        print(t)
