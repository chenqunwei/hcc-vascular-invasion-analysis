import re, sys
from pathlib import Path
from docx import Document

SRC = Path(r"E:/workbuddy/HCC pure bioinformatics analysis/HCC血管侵犯连续谱_手稿_v5.docx")
doc = Document(str(SRC))

paras = doc.paragraphs
print("TOTAL PARAGRAPHS:", len(paras))

# Heuristic: find the References section start
# Usually a heading "References" / "参考文献" near the end.
ref_start_idx = None
for i, p in enumerate(paras):
    t = p.text.strip()
    if re.match(r'^(references|参考文献|reference list)\s*$', t, re.I):
        ref_start_idx = i
        break
print("REFERENCES HEADING IDX:", ref_start_idx)

# In-text citation markers across whole doc
cit_re = re.compile(r'\[(\d+(?:\s*[-,]\s*\d+)?)\]')
all_citations = []
for i, p in enumerate(paras):
    for m in cit_re.finditer(p.text):
        all_citations.append((i, m.group(0), m.group(1)))

# Count of plain [n] single-numbered
single = [g for (_,_,g) in all_citations if re.fullmatch(r'\d+', g.replace(' ','')) and ',' not in g and '-' not in g]
ranges = [(lbl,grp) for (_,lbl,grp) in all_citations if (',' in grp or '-' in grp)]
print("TOTAL IN-TEXT CITATION MARKERS (incl ranges):", len(all_citations))
print("SINGLE-NUM MARKERS:", len(single))
print("RANGE/COMBO MARKERS:", len(ranges), ranges[:20])

# Body citations (before references heading)
body_cits = [(i,lbl,grp) for (i,lbl,grp) in all_citations if ref_start_idx is None or i < ref_start_idx]
print("BODY CITATION MARKERS:", len(body_cits))

# Reference entries: paragraphs after ref_start_idx
if ref_start_idx is not None:
    ref_paras = paras[ref_start_idx+1:]
    # A reference entry often starts with a number "1." or "[1]" or just "1 "
    entries = []
    buf = []
    for p in ref_paras:
        t = p.text.rstrip()
        if not t.strip():
            if buf:
                entries.append("\n".join(buf)); buf=[]
            continue
        # new entry if starts with a leading number pattern OR clearly a new ref
        if re.match(r'^\s*\[?\d+\]?[\.\)]?\s+', t) and buf:
            entries.append("\n".join(buf)); buf=[t]
        else:
            buf.append(t)
    if buf:
        entries.append("\n".join(buf))
    print("REFERENCE ENTRIES (heuristic):", len(entries))
    for j,e in enumerate(entries[:40],1):
        print(f"--- ENTRY {j} ---")
        print(e[:600])
