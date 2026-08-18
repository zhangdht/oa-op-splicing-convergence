# -*- coding: utf-8 -*-
"""
08b_verify_references_round2.py
第二轮：补齐数据集来源文献、方法学引用，并剔除第一轮中主题不符的宽松命中。
"""
import json, time, urllib.parse, urllib.request, csv, os

ROOT = r"D:/bioinfo05"
OUT_DIR = os.path.join(ROOT, "06_paper")
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# 第一轮中主题不符、必须剔除的条目
DROP = {"tissue_heterogeneity_transcriptomics", "reproducibility_biomarker_signature",
        "benjamini_hochberg"}

QUERIES = [
    ("benjamini2001_fdr",
     'AUTH:"Benjamini Y" AND TITLE:"false discovery rate" AND SRC:MED', ["false discovery rate"]),
    ("ramos2014_gse57218",
     '(TITLE:"Genes involved in the osteoarthritis process identified through genome wide expression analysis in articular cartilage") AND SRC:MED',
     ["osteoarthritis", "cartilage"]),
    ("liu2005_gse56815",
     '((TITLE:"circulating monocytes" OR TITLE:"peripheral blood monocytes") AND TITLE:"bone mineral density") AND SRC:MED',
     ["monocyte"]),
    ("xiao2008_gse7158",
     '(TITLE:"peak bone mass" AND (TITLE:"gene expression" OR TITLE:"monocyte")) AND SRC:MED',
     ["bone"]),
    ("subramanian2005_gsea",
     'TITLE:"Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles" AND SRC:MED',
     ["Gene set enrichment analysis"]),
    ("ashburner2000_go",
     'TITLE:"Gene ontology: tool for the unification of biology" AND SRC:MED', ["Gene ontology"]),
    ("jassal2020_reactome",
     'TITLE:"The reactome pathway knowledgebase" AND SRC:MED', ["reactome"]),
    ("robin2011_proc",
     'TITLE:"pROC: an open-source package for R and S+ to analyze and compare ROC curves" AND SRC:MED',
     ["pROC"]),
    ("delong1988_auc",
     'TITLE:"Comparing the areas under two or more correlated receiver operating characteristic curves" AND SRC:MED',
     ["receiver operating characteristic"]),
    ("op_oa_inverse_meta",
     '((TITLE:"osteoarthritis" AND TITLE:"osteoporosis") AND (TITLE:"meta-analysis" OR TITLE:"systematic review")) AND SRC:MED',
     ["osteoarthritis", "osteoporosis"]),
    ("hip_fracture_oa_risk",
     '(TITLE:"osteoarthritis" AND TITLE:"fracture" AND TITLE:"risk") AND SRC:MED',
     ["osteoarthritis", "fracture"]),
    ("msc_senescence_splicing",
     '((TITLE:"mesenchymal" OR TITLE:"stem cell") AND TITLE:"senescence" AND (TITLE:"splicing" OR TITLE:"spliceosome")) AND SRC:MED',
     ["senescence"]),
    ("son_nuclear_speckle_stem",
     '(TITLE:"SON" AND (TITLE:"nuclear speckle" OR TITLE:"pluripotency" OR TITLE:"stem cell")) AND SRC:MED',
     ["SON"]),
    ("oa_cartilage_rnaseq_integrative",
     '((TITLE:"osteoarthritis") AND (TITLE:"integrative" OR TITLE:"integrated") AND (TITLE:"analysis")) AND SRC:MED',
     ["osteoarthritis"]),
    ("bone_marrow_msc_osteoporosis_transcriptome",
     '((TITLE:"bone marrow" AND TITLE:"osteoporosis") AND (TITLE:"transcriptom" OR TITLE:"gene expression")) AND SRC:MED',
     ["osteoporosis"]),
    ("cross_tissue_gene_expression_variation",
     'TITLE:"The Genotype-Tissue Expression (GTEx) project" AND SRC:MED', ["GTEx"]),
    ("wilkinson2016_fair",
     'TITLE:"The FAIR Guiding Principles for scientific data management and stewardship" AND SRC:MED',
     ["FAIR"]),
    ("hedges_effect_size",
     '(TITLE:"effect size" AND (TITLE:"meta-analysis" OR TITLE:"estimation")) AND SRC:MED',
     ["effect size"]),
    ("splicing_dysregulation_disease",
     '(TITLE:"alternative splicing" AND TITLE:"disease") AND SRC:MED', ["splicing"]),
    ("osteoblast_adipocyte_switch",
     '((TITLE:"osteoblast" AND TITLE:"adipocyte") AND (TITLE:"switch" OR TITLE:"lineage" OR TITLE:"balance")) AND SRC:MED',
     ["adipocyte"]),
]


def epmc_search(query, page_size=6):
    params = {"query": query, "format": "json", "resultType": "core",
              "pageSize": str(page_size), "sort": "CITED desc"}
    url = EPMC + "?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as f:
                return json.load(f)
        except Exception:
            if attempt == 3:
                return None
            time.sleep(3)
    return None


def pick(res, must):
    if not res:
        return None
    for it in res.get("resultList", {}).get("result", []):
        if it.get("pmid") and all(k.lower() in (it.get("title") or "").lower() for k in must):
            return it
    for it in res.get("resultList", {}).get("result", []):
        if it.get("pmid"):
            it["_loose"] = True
            return it
    return None


new_rows = []
for key, q, must in QUERIES:
    it = pick(epmc_search(q), must)
    if it is None:
        print("[MISS]", key)
        continue
    ji = it.get("journalInfo", {}) or {}
    row = {"key": key, "pmid": it.get("pmid", ""), "doi": it.get("doi", ""),
           "title": (it.get("title") or "").rstrip("."),
           "authors": it.get("authorString", ""),
           "journal": (ji.get("journal", {}) or {}).get("title", ""),
           "year": it.get("pubYear", ""), "volume": ji.get("volume", ""),
           "issue": ji.get("issue", ""), "pages": it.get("pageInfo", ""),
           "cited_by": it.get("citedByCount", ""),
           "loose_match": "YES" if it.get("_loose") else "", "query": q}
    new_rows.append(row)
    print("[OK]  %-42s PMID:%-9s %s (%s) %s" % (key, row["pmid"], row["journal"][:32],
                                                row["year"], row["title"][:58]))
    time.sleep(0.4)

# 合并第一轮
p1 = os.path.join(OUT_DIR, "references_verified.csv")
old = [r for r in csv.DictReader(open(p1, encoding="utf-8-sig")) if r["key"] not in DROP]
allr = old + new_rows
# 去重（同 PMID 只留一条）
seen, ded = set(), []
for r in allr:
    if r["pmid"] in seen:
        continue
    seen.add(r["pmid"]); ded.append(r)

out = os.path.join(OUT_DIR, "references_verified.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(ded[0].keys()))
    w.writeheader(); w.writerows(ded)
print("\n最终真实文献 %d 条（已剔除 %d 条主题不符）-> %s" % (len(ded), len(DROP), out))
