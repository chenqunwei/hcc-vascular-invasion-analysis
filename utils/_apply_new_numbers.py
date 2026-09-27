# -*- coding: utf-8 -*-
"""Apply the NEW deterministic numbers to v3.md (plain text) and v5.docx (in place).
Changes ONLY numeric/label values; preserves all formatting and reference annotations in v5.docx.
"""
import re, io
import docx
from docx.document import Document as _Doc
from docx.table import Table
from docx.text.paragraph import Paragraph
from docx.oxml.ns import qn

EN = "\u2013"   # en dash
CHI2 = "\u03c7\u00b2"  # chi-squared symbol

DOCX = "HCC血管侵犯连续谱_手稿_v5.docx"
MD = "HCC血管侵犯连续谱_手稿_v3.md"

# ----------------------------------------------------------------------
# 1) v3.md plain-text replacements (full strings, safe)
# ----------------------------------------------------------------------
md_repl = [
    ("Low VI (n = 100), Intermediate VI (n = 128), and High VI (n = 71)",
     "Low VI (n = 99), Intermediate VI (n = 128), and High VI (n = 72)"),
    ("(499, 705, and 914, respectively)", "(497, 705, and 914, respectively)"),
    ("21.0%, 32.8%, and 42.3%", "21.2%, 32.8%, and 41.7%"),
    ("21.0% (Low), 32.8% (Intermediate), and 42.3% (High VI) (" + CHI2 + " = 9.06, P = 0.011)",
     "21.2% (Low), 32.8% (Intermediate), and 41.7% (High VI) (" + CHI2 + " = 8.44, P = 0.0147)"),
    ("median 1,412 ng/mL", "median 1,456 ng/mL"),
    ("Low (12.5 ng/mL)", "Low (12.0 ng/mL)"),
    ("28.4% vs 14.1% and 17.2%, P = 0.064", "27.9% vs 14.3% and 17.2%, P = 0.077"),
    ("50.7% vs 68.0% and 74.2%, P = 0.003", "51.4% vs 67.7% and 74.2%, P = 0.004"),
    ("1.97-fold higher recurrence risk", "1.82-fold higher recurrence risk"),
    ("HR = 1.97, 95% CI 1.16" + EN + "3.33, P = 0.012",
     "HR = 1.82, 95% CI 1.08" + EN + "3.05, P = 0.024"),
    ("HR = 1.43 per subtype level, 95% CI 1.12" + EN + "1.82, P = 0.0043",
     "HR = 1.38 per subtype level, 95% CI 1.08" + EN + "1.76, P = 0.0091"),
    ("rose monotonically from 21% in the Low VI subtype to 42% in the High VI subtype, and recurrence risk increased by 43% per subtype level (HR = 1.43, P = 0.0043)",
     "rose monotonically from 21.2% in the Low VI subtype to 41.7% in the High VI subtype, and recurrence risk increased by 38% per subtype level (HR = 1.38, P = 0.0091)"),
    ("median 1,412 vs 12.5 ng/mL", "median 1,456 vs 12.0 ng/mL"),
]

# docx run-level replacements (avoid χ² multi-run issue by matching the numeric tail)
docx_repl = [
    ("Low VI (n = 100), Intermediate VI (n = 128), and High VI (n = 71)",
     "Low VI (n = 99), Intermediate VI (n = 128), and High VI (n = 72)"),
    ("(499, 705, and 914, respectively)", "(497, 705, and 914, respectively)"),
    ("21.0%, 32.8%, and 42.3%", "21.2%, 32.8%, and 41.7%"),
    ("21.0% (Low), 32.8% (Intermediate), and 42.3% (High VI)",
     "21.2% (Low), 32.8% (Intermediate), and 41.7% (High VI)"),
    ("9.06, P = 0.011", "8.44, P = 0.0147"),   # covers both inline + Fig3 legend (χ² prefix untouched)
    ("median 1,412 ng/mL", "median 1,456 ng/mL"),
    ("Low (12.5 ng/mL)", "Low (12.0 ng/mL)"),
    ("28.4% vs 14.1% and 17.2%, P = 0.064", "27.9% vs 14.3% and 17.2%, P = 0.077"),
    ("50.7% vs 68.0% and 74.2%, P = 0.003", "51.4% vs 67.7% and 74.2%, P = 0.004"),
    ("1.97-fold higher recurrence risk", "1.82-fold higher recurrence risk"),
    ("HR = 1.97, 95% CI 1.16" + EN + "3.33, P = 0.012",
     "HR = 1.82, 95% CI 1.08" + EN + "3.05, P = 0.024"),
    ("HR = 1.43 per subtype level, 95% CI 1.12" + EN + "1.82, P = 0.0043",
     "HR = 1.38 per subtype level, 95% CI 1.08" + EN + "1.76, P = 0.0091"),
    ("rose monotonically from 21% in the Low VI subtype to 42% in the High VI subtype, and recurrence risk increased by 43% per subtype level (HR = 1.43, P = 0.0043)",
     "rose monotonically from 21.2% in the Low VI subtype to 41.7% in the High VI subtype, and recurrence risk increased by 38% per subtype level (HR = 1.38, P = 0.0091)"),
    ("median 1,412 vs 12.5 ng/mL", "median 1,456 vs 12.0 ng/mL"),
]

def apply_run_repl(p, repl):
    """Run-level replacement preserving formatting. Returns count of applied edits."""
    cnt = 0
    for old, new in repl:
        for run in p.runs:
            if old in run.text:
                run.text = run.text.replace(old, new)
                cnt += 1
    return cnt

# ----------------------------------------------------------------------
# 2) DOCX table rewrites (authoritative NEW values)
# ----------------------------------------------------------------------
def set_cell(cell, text):
    cell.text = ""
    p = cell.paragraphs[0]
    run = p.add_run(text)

# Table 1 new values
T1_HEADER = ["Variable", "Low VI (n=99)", "Intermediate VI (n=128)", "High VI (n=72)", "P"]
T1_ROWS = [
    ["Age, median (years)", "60.0", "64.0", "59.0", "0.065"],
    ["Male, n (%)", "67 (67.7%)", "95 (74.2%)", "37 (51.4%)", "0.004"],
    ["Stage III/IV, n (%)", "13 (14.3%)", "21 (17.2%)", "19 (27.9%)", "0.077"],
    ["MVI positive, n (%)", "21 (21.2%)", "42 (32.8%)", "30 (41.7%)", "0.015"],
    ["AFP, median (ng/mL)", "12.0", "7.0", "1456.0", "<0.0001"],
    ["Child-Pugh A, n (%)", "64 (94.1%)", "94 (88.7%)", "41 (91.1%)", "0.477"],
]
# Table 2 new values (rows 1-6)
T2_HEADER = ["Variable", "DFS HR (95% CI)", "DFS P", "OS HR (95% CI)", "OS P"]
T2_ROWS = [
    ["VI score, per 1 SD", "1.35 (1.12" + EN + "1.64)", "0.0020", "1.34 (1.07" + EN + "1.67)", "0.0096"],
    ["VI subtype, High vs Low", "1.82 (1.08" + EN + "3.05)", "0.024", "1.34 (0.76" + EN + "2.38)", "0.310"],
    ["VI subtype, per level", "1.38 (1.08" + EN + "1.76)", "0.0091", "1.14 (0.85" + EN + "1.52)", "0.383"],
    ["Advanced stage (III/IV vs I/II)", "2.16 (1.41" + EN + "3.30)", "0.0004", "2.11 (1.33" + EN + "3.34)", "0.0015"],
    ["Age", "1.00 (0.99" + EN + "1.01)", "0.979", "1.02 (1.00" + EN + "1.04)", "0.041"],
    ["Sex, Male", "0.91 (0.62" + EN + "1.33)", "0.618", "0.81 (0.52" + EN + "1.25)", "0.331"],
]

def rewrite_table(table, header, rows):
    for j, val in enumerate(header):
        set_cell(table.rows[0].cells[j], val)
    for i, row in enumerate(rows, start=1):
        for j, val in enumerate(row):
            set_cell(table.rows[i].cells[j], val)

def find_tables(doc):
    t1 = t2 = None
    for t in doc.tables:
        r0 = [c.text.strip() for c in t.rows[0].cells]
        if r0 and r0[0] == "Variable":
            t1 = t
        if any("DFS HR" in c for c in r0) or any("VI score, per 1 SD" in c for c in r0):
            t2 = t
    return t1, t2

# ----------------------------------------------------------------------
# EXECUTE
# ----------------------------------------------------------------------
# v3.md
txt = open(MD, encoding="utf-8").read()
md_before = txt
for old, new in md_repl:
    if old not in txt:
        print("[MD-WARN] not found:", repr(old[:40]))
    txt = txt.replace(old, new)
open(MD, "w", encoding="utf-8").write(txt)
print("[MD] applied", sum(1 for _ in md_repl), "replacements")

# v5.docx
doc = docx.Document(DOCX)
total_run_edits = 0
for p in doc.paragraphs:
    total_run_edits += apply_run_repl(p, docx_repl)
t1, t2 = find_tables(doc)
print("[DOCX] tables found: T1=", t1 is not None, "T2=", t2 is not None)
if t1:
    rewrite_table(t1, T1_HEADER, T1_ROWS)
    print("[DOCX] Table 1 rewritten")
if t2:
    rewrite_table(t2, T2_HEADER, T2_ROWS)
    print("[DOCX] Table 2 rewritten")
doc.save(DOCX)
print("[DOCX] run-level edits applied:", total_run_edits)
print("DONE")
