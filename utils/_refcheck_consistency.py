import re
from pathlib import Path
from docx import Document

SRC = Path(r"E:/workbuddy/HCC pure bioinformatics analysis/HCC血管侵犯连续谱_手稿_v5.docx")
doc = Document(str(SRC))
paras = doc.paragraphs

# References heading
ref_start = None
for i,p in enumerate(paras):
    if re.match(r'^(references|参考文献)\s*$', p.text.strip(), re.I):
        ref_start = i; break

cit_re = re.compile(r'\[(\d+)\]')
body_cits = []   # (para_idx, num)
for i,p in enumerate(paras):
    if ref_start is not None and i >= ref_start: 
        continue
    for m in cit_re.finditer(p.text):
        body_cits.append((i, int(m.group(1))))

nums = [n for _,n in body_cits]
print("Body citation count:", len(nums))
print("Body citation number set:", sorted(set(nums)))
print("Min/max cited:", min(nums), max(nums))
print("Any num >33 or <1:", [n for n in nums if n>33 or n<1])

# First-appearance order
first_seen = {}
for i,n in body_cits:
    first_seen.setdefault(n, i)
order = sorted(first_seen.items(), key=lambda kv: kv[1])
fa_order = [n for n,_ in order]
print("\nFirst-appearance order (as text progresses):")
print(fa_order)
# check ascending
asc_ok = all(fa_order[i] < fa_order[i+1] for i in range(len(fa_order)-1))
print("Strictly ascending first-appearance order:", asc_ok)
if not asc_ok:
    for i in range(len(fa_order)-1):
        if fa_order[i] > fa_order[i+1]:
            print(f"  OUT-OF-ORDER at pos {i}: {fa_order[i]} appears before {fa_order[i+1]}")

# Orphan check: refs 1..33 all cited in body?
missing = [n for n in range(1,34) if n not in set(nums)]
print("\nReference numbers NOT cited in body (orphans):", missing)

# Count per ref
from collections import Counter
c = Counter(nums)
print("\nCitation frequency per ref (ref:count):")
for n in range(1,34):
    print(f"  [{n}]: {c.get(n,0)}")
