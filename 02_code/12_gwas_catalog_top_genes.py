# -*- coding: utf-8 -*-
"""
12_gwas_catalog_top_genes.py

不下载 GB 级汇总统计，改用 Europe PMC 检索剪接轴/WGCNA 交集基因在骨病 GWAS 中
是否有已发表报道。仅用于论文的"遗传支持"描述，不做正式 MR。
"""
import os, csv, json, time, urllib.parse, urllib.request

ROOT = r"D:/bioinfo05"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

os.makedirs(os.path.join(ROOT, "03_results/genetics"), exist_ok=True)

GENES = ["SON", "MMP9", "PIK3CG", "HLA-DMB", "FCGR2A", "ADM", "BIRC5",
         "ALOX5AP", "ARHGAP15", "VNN2", "SNX10", "ATP8B4", "FSTL1", "KLF5",
         "EFTUD2", "SF3B1", "SF3B2", "SF3A1", "PRPF8", "SNRNP200", "U2AF1",
         "SRSF1", "HNRNPA1", "HNRNPA1L2"]

TRAITS = ["bone mineral density", "osteoporosis", "fracture", "osteoarthritis",
          "hip fracture", "BMD"]


def search(query, page_size=10):
    url = EPMC + "?" + urllib.parse.urlencode({
        "query": query, "format": "json", "resultType": "core",
        "pageSize": str(page_size), "sort": "CITED desc"})
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=45) as f:
                return json.load(f)
        except Exception:
            time.sleep(2)
    return None


rows = []
for gene in GENES:
    q = f'({gene} AND ("genome-wide association" OR "GWAS") AND ({" OR ".join(f"\"{t}\"" for t in TRAITS)})) AND SRC:MED'
    res = search(q, page_size=8)
    hits = res.get("resultList", {}).get("result", []) if res else []
    for it in hits:
        rows.append({
            "gene": gene,
            "pmid": it.get("pmid"),
            "title": it.get("title", "").rstrip("."),
            "journal": (it.get("journalInfo", {}).get("journal", {}) or {}).get("title", ""),
            "year": it.get("pubYear"),
            "query": q
        })
    print(f"{gene}: {len(hits)} hits")
    time.sleep(0.35)

out = os.path.join(ROOT, "03_results/genetics", "gene_trait_literature_hits.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["gene", "pmid", "title", "journal", "year", "query"])
    w.writeheader(); w.writerows(rows)

print(f"\nWrote {len(rows)} hits to {out}")
