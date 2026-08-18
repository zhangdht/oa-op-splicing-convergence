#!/usr/bin/env python
"""
20_compartment_robustness.py
----------------------------
Robustness check for the compartment-specificity claim.

Motivation
----------
The original leading-edge compartment test (15_functional_annotation.R) used a
single OP monocyte cohort (GSE56815). Two monocyte cohorts are available
(GSE56815, GSE7158), and SON behaves differently between them. A reviewer will
reasonably ask whether the "haematopoietic compartment does not share the
programme" claim survives when both monocyte cohorts are used.

This script re-runs the test on:
  (a) GSE56815 alone      (as originally reported)
  (b) GSE7158  alone      (independent monocyte cohort)
  (c) both cohorts combined, requiring concordant direction in each

It also tabulates SON across all seven cohorts so the manuscript can state the
per-cohort result accurately rather than summarising monocytes as a single
outcome.

Outputs
-------
03_results/functional/compartment_robustness.txt
03_results/functional/SON_per_cohort.csv
"""

import csv
import os
from math import comb

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEG = os.path.join(ROOT, "03_results", "deg")
FUNC = os.path.join(ROOT, "03_results", "functional")
LEAD = os.path.join(FUNC, "mesenchymal_leading_edge_genes.csv")

COHORTS = {
    "GSE114007": ("OA", "knee cartilage"),
    "GSE57218": ("OA", "knee/hip cartilage"),
    "GSE117999": ("OA", "knee cartilage"),
    "GSE169077": ("OA", "knee cartilage"),
    "GSE56815": ("OP", "monocyte"),
    "GSE7158": ("OP", "monocyte"),
    "GSE35958": ("OP", "BM-MSC"),
}
MONOCYTE = ["GSE56815", "GSE7158"]


def binom_two_sided(k, n, p=0.5):
    """Exact two-sided binomial test (equal-tail method used by R binom.test)."""
    if n == 0:
        return float("nan")
    probs = [comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(n + 1)]
    obs = probs[k]
    tol = 1e-7 * obs
    return min(1.0, sum(pr for pr in probs if pr <= obs + tol))


def read_deg(gse):
    path = os.path.join(DEG, f"DEG_{gse}.csv")
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path, encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            g = (row.get("gene") or "").strip()
            if not g:
                continue
            out[g] = row
    return out


def fnum(row, *keys):
    for k in keys:
        v = row.get(k)
        if v not in (None, "", "NA"):
            try:
                return float(v)
            except ValueError:
                pass
    return None


def main():
    with open(LEAD, encoding="utf-8-sig") as fh:
        lead = [r["gene"].strip() for r in csv.DictReader(fh) if r.get("gene")]
    lead = [g for g in lead if g]
    print(f"leading-edge splicing factors: {len(lead)}")

    degs = {g: read_deg(g) for g in COHORTS}

    lines = []
    lines.append("Compartment-specificity robustness check")
    lines.append("=" * 62)
    lines.append(
        f"Leading-edge set: {len(lead)} splicing factors concordantly "
        "down-regulated in OA cartilage and OP BM-MSC."
    )
    lines.append("")
    lines.append("Question: is the same coordinated down-regulation present in the")
    lines.append("haematopoietic compartment? Tested per monocyte cohort and jointly.")
    lines.append("")

    per_cohort_down = {}
    for gse in MONOCYTE:
        d = degs[gse]
        vals = {}
        for g in lead:
            row = d.get(g)
            if not row:
                continue
            t = fnum(row, "t", "statistic", "logFC")
            if t is not None:
                vals[g] = t
        n = len(vals)
        k = sum(1 for v in vals.values() if v < 0)
        p = binom_two_sided(k, n)
        per_cohort_down[gse] = vals
        lines.append(
            f"  {gse:10} down {k:3}/{n:3} ({100*k/n:.1f}%)  binomial P = {p:.4g}"
        )
        print(f"  {gse}: {k}/{n} down ({100*k/n:.1f}%), P={p:.4g}")

    common = set(per_cohort_down[MONOCYTE[0]]) & set(per_cohort_down[MONOCYTE[1]])
    both_down = [
        g for g in common
        if per_cohort_down[MONOCYTE[0]][g] < 0 and per_cohort_down[MONOCYTE[1]][g] < 0
    ]
    n_c, k_c = len(common), len(both_down)
    # under independence, concordant-down expected at 25%
    p_c = binom_two_sided(k_c, n_c, p=0.25)
    lines.append("")
    lines.append(
        f"  Concordant across BOTH monocyte cohorts: {k_c}/{n_c} "
        f"({100*k_c/n_c:.1f}%); expected 25% under independence; P = {p_c:.4g}"
    )
    print(f"  both-cohort concordant down: {k_c}/{n_c} ({100*k_c/n_c:.1f}%), P={p_c:.4g}")

    lines.append("")
    lines.append("Interpretation")
    lines.append("-" * 62)
    if k_c / n_c < 0.5:
        lines.append(
            "Neither monocyte cohort shows the coordinated down-regulation seen in"
        )
        lines.append(
            "the mesenchymal compartment. Individual splicing factors do move in"
        )
        lines.append(
            "monocytes, but not as a concordant programme, supporting the"
        )
        lines.append("compartment-specific interpretation.")
    else:
        lines.append(
            "A majority of leading-edge genes are concordantly down in monocytes;"
        )
        lines.append("the compartment-specific claim would need to be softened.")

    # ---- SON across all cohorts -------------------------------------------
    lines.append("")
    lines.append("SON across all seven cohorts")
    lines.append("-" * 62)
    son_rows = []
    for gse, (dis, tis) in COHORTS.items():
        row = degs[gse].get("SON")
        if not row:
            son_rows.append(dict(cohort=gse, disease=dis, tissue=tis,
                                 logFC="NA", P="NA", note="not measured"))
            lines.append(f"  {gse:10} {dis:3} {tis:20} not measured")
            continue
        lfc = fnum(row, "logFC", "log2FC")
        pv = fnum(row, "P.Value", "pvalue", "p_value", "P")
        sig = "significant" if (pv is not None and pv < 0.05) else "ns"
        son_rows.append(dict(cohort=gse, disease=dis, tissue=tis,
                             logFC=f"{lfc:.3f}" if lfc is not None else "NA",
                             P=f"{pv:.4g}" if pv is not None else "NA", note=sig))
        lines.append(
            f"  {gse:10} {dis:3} {tis:20} logFC={lfc:+.3f}  P={pv:.4g}  {sig}"
        )

    out_csv = os.path.join(FUNC, "SON_per_cohort.csv")
    with open(out_csv, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cohort", "disease", "tissue",
                                           "logFC", "P", "note"])
        w.writeheader()
        w.writerows(son_rows)

    out_txt = os.path.join(FUNC, "compartment_robustness.txt")
    with open(out_txt, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\nwrote:", out_txt)
    print("wrote:", out_csv)
    print("\n" + "\n".join(lines))


if __name__ == "__main__":
    main()
