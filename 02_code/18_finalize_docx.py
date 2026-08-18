#!/usr/bin/env python
"""
18_finalize_docx.py
-------------------
Finalize the generated Manuscript_final.docx by ensuring every <w:tbl> carries a
<w:tblGrid> element. The flextable::body_add_flextable path used by
16_write_paper_final.R emits tables with autofit/pct layout but omits w:tblGrid,
which strict OOXML validators (and some Word versions) flag as corrupt. This
post-processor injects a w:tblGrid with one w:gridCol per first-row cell so the
document opens reliably everywhere.

Usage:
    python 18_finalize_docx.py [input.docx] [output.docx]

Defaults: input = 06_paper/Manuscript_final.docx, output = same file (in place).
A backup with suffix .bak_tblgrid is written next to the target.
"""
import sys, os, re, zipfile, shutil

def fix_table(t: str) -> str:
    if "<w:tblGrid>" in t:
        return t
    rows = re.findall(r"<w:tr\b.*?</w:tr>", t, re.S)
    if not rows:
        return t
    cells = re.findall(r"<w:tc\b.*?</w:tc>", rows[0], re.S)
    n = len(cells)
    widths = []
    for c in cells:
        m = re.search(r'<w:tcW\s+[^>]*?w:w="(\d+)"[^>]*?w:type="(dxa|pct|auto)"', c)
        widths.append(m.group(1) if (m and m.group(2) == "dxa") else "0")
    if len(widths) != n:
        widths = ["0"] * n
    grid = "<w:tblGrid>" + "".join(f'<w:gridCol w:w="{w}"/>' for w in widths) + "</w:tblGrid>"
    idx = t.find("</w:tblPr>")
    if idx == -1:
        return t
    at = idx + len("</w:tblPr>")
    return t[:at] + grid + t[at:]

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    inp = sys.argv[1] if len(sys.argv) > 1 else os.path.join(root, "06_paper", "Manuscript_final.docx")
    out = sys.argv[2] if len(sys.argv) > 2 else inp
    if not os.path.exists(inp):
        sys.exit(f"Input not found: {inp}")
    bak = inp + ".bak_tblgrid"
    if not os.path.exists(bak):
        shutil.copy(inp, bak)
    z = zipfile.ZipFile(inp)
    doc_xml = z.read("word/document.xml").decode("utf-8", "ignore")
    tbls = re.findall(r"<w:tbl\b.*?</w:tbl>", doc_xml, re.S)
    new_doc = doc_xml
    for t in reversed(tbls):
        new_doc = new_doc.replace(t, fix_table(t), 1)
    parts = {n: z.read(n) for n in z.namelist()}
    parts["word/document.xml"] = new_doc.encode("utf-8")
    tmp = out + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z2:
        for n, data in parts.items():
            z2.writestr(n, data)
    z.close()
    os.replace(tmp, out)
    # verify
    zv = zipfile.ZipFile(out)
    assert zv.testzip() is None, "ZIP integrity failed after rewrite"
    from docx import Document
    d = Document(out)
    print(f"OK: {out} | tables={len(d.tables)} | tblGrid injected into {len(tbls)} table(s)")

if __name__ == "__main__":
    main()
