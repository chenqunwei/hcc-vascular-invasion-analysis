# -*- coding: utf-8 -*-
"""Dump v5.docx paragraph + table cell texts (ASCII-safe) for inspection."""
import docx
from docx.document import Document as _Doc
from docx.table import Table
from docx.text.paragraph import Paragraph

doc = docx.Document("HCC血管侵犯连续谱_手稿_v5.docx")

def iter_block_items(parent):
    from docx.oxml.ns import qn
    for child in parent.element.body.iterchildren():
        if child.tag == qn('w:p'):
            yield Paragraph(child, parent)
        elif child.tag == qn('w:tbl'):
            yield Table(child, parent)

import sys
out = []
tbl_idx = 0
for blk in iter_block_items(doc):
    if isinstance(blk, Paragraph):
        t = blk.text.strip()
        if t:
            out.append("P: " + t)
    else:
        tbl_idx += 1
        out.append(f"=== TABLE {tbl_idx} ===")
        for r, row in enumerate(blk.rows):
            cells = [c.text.strip().replace("\n"," / ") for c in row.cells]
            out.append(f"  R{r}: " + " | ".join(cells))

with open("_docx_dump.txt","w",encoding="ascii",errors="replace") as f:
    f.write("\n".join(out))
print("dumped", len(out), "blocks")
