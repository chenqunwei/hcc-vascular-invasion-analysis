# -*- coding: utf-8 -*-
"""Print key sentences from the FIXED v5 docx to visually confirm correctness."""
import docx
from docx.text.paragraph import Paragraph
from docx.table import Table

doc = docx.Document("HCC血管侵犯连续谱_手稿_v5_fixed.docx")

def iter_block_items(parent):
    from docx.oxml.ns import qn
    for child in parent.element.body.iterchildren():
        if child.tag == qn('w:p'):
            yield Paragraph(child, parent)
        elif child.tag == qn('w:tbl'):
            yield Table(child, parent)

KW = ["VI (n", "median VI score", "MVI-positive tumors increased", "median 1,",
      "Low (12", "Stage III/IV tumors tended", "Sex was modestly",
      "1.82-fold", "HR = 1.82", "HR = 1.38", "rose monotonically",
      "median 1,456 vs", "χ² = 8.44", "Table 1", "Table 2"]

seen = set()
for blk in iter_block_items(doc):
    if isinstance(blk, Paragraph):
        t = blk.text.strip()
        if t and any(k in t for k in KW):
            if t not in seen:
                seen.add(t)
                print("P:", t)
    else:
        for ri, row in enumerate(blk.rows):
            cells = [c.text.strip() for c in row.cells]
            line = " | ".join(cells)
            if any(k in line for k in ["Low VI", "VI score, per", "VI subtype", "Age,", "Male,", "Stage", "MVI positive", "AFP,", "Child-Pugh"]):
                print(f"T-R{ri}:", line)
