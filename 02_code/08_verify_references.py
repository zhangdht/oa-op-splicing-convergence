# -*- coding: utf-8 -*-
"""
08_verify_references.py

原则：论文中的每一条参考文献都必须是从 Europe PMC / PubMed 实时检索回来的真实记录，
      绝不允许凭记忆书写 PMID 或标题。

做法：给出「检索式」而非「PMID」。脚本执行检索，取回真实记录（PMID/DOI/期刊/年份/作者），
      人工（在报告里）核对返回条目是否确为所需文献；不匹配的条目会被标记，不会静默使用。

输出：06_paper/references_verified.csv
      06_paper/references_verification_log.txt
"""
import json
import time
import urllib.parse
import urllib.request
import csv
import os

ROOT = r"D:/bioinfo05"
OUT_DIR = os.path.join(ROOT, "06_paper")
os.makedirs(OUT_DIR, exist_ok=True)

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# (内部标识, 检索式, 期望关键词——用于自动校验返回条目是否对得上)
QUERIES = [
    ("morris2019_ebmd",
     'TITLE:"An atlas of genetic influences on osteoporosis in humans and mice" AND SRC:MED',
     ["atlas", "osteoporosis"]),
    ("tachmazidou2019_oa",
     'TITLE:"Identification of new therapeutic targets for osteoarthritis through genome-wide analyses of UK Biobank data" AND SRC:MED',
     ["osteoarthritis", "UK Biobank"]),
    ("hartley_bmd_oa_mr",
     '(TITLE:"osteoarthritis" AND TITLE:"bone mineral density" AND TITLE:"Mendelian randomization") AND SRC:MED',
     ["osteoarthritis", "bone mineral density"]),
    ("langfelder2008_wgcna",
     'TITLE:"WGCNA: an R package for weighted correlation network analysis" AND SRC:MED',
     ["WGCNA"]),
    ("hanzelmann2013_gsva",
     'TITLE:"GSVA: gene set variation analysis for microarray and RNA-seq data" AND SRC:MED',
     ["GSVA"]),
    ("ritchie2015_limma",
     'TITLE:"limma powers differential expression analyses for RNA-sequencing and microarray studies" AND SRC:MED',
     ["limma"]),
    ("liberzon2015_hallmark",
     'TITLE:"The Molecular Signatures Database Hallmark Gene Set Collection" AND SRC:MED',
     ["Hallmark", "Molecular Signatures"]),
    ("hao2021_seurat",
     'TITLE:"Integrated analysis of multimodal single-cell data" AND SRC:MED',
     ["multimodal single-cell"]),
    ("friedman2010_glmnet",
     'TITLE:"Regularization Paths for Generalized Linear Models via Coordinate Descent" AND SRC:MED',
     ["Regularization Paths"]),
    ("fisch2018_gse114007",
     '(TITLE:"Identification of transcription factors responsible for dysregulated networks in human osteoarthritis cartilage") AND SRC:MED',
     ["osteoarthritis cartilage", "transcription factor"]),
    ("benisch2012_gse35958",
     '(TITLE:"The transcriptional profile of mesenchymal stem cell populations in primary osteoporosis") AND SRC:MED',
     ["mesenchymal stem cell", "osteoporosis"]),
    ("son_ztkk_syndrome",
     '(TITLE:"SON" AND TITLE:"ZTTK") AND SRC:MED',
     ["SON"]),
    ("son_splicing_function",
     '(TITLE:"SON" AND (TITLE:"splicing" OR TITLE:"spliceosome")) AND SRC:MED',
     ["SON"]),
    ("splicing_ageing_harries",
     '(TITLE:"alternative splicing" AND (TITLE:"ageing" OR TITLE:"aging" OR TITLE:"senescence")) AND SRC:MED',
     ["splicing"]),
    ("splicing_factor_senescence",
     '(TITLE:"splicing factor" AND TITLE:"senescence") AND SRC:MED',
     ["splicing factor", "senescence"]),
    ("op_oa_comorbidity",
     '(TITLE:"osteoporosis" AND TITLE:"osteoarthritis" AND (TITLE:"inverse" OR TITLE:"association" OR TITLE:"relationship")) AND SRC:MED',
     ["osteoporosis", "osteoarthritis"]),
    ("subchondral_bone_oa",
     '(TITLE:"subchondral bone" AND TITLE:"osteoarthritis") AND SRC:MED',
     ["subchondral bone", "osteoarthritis"]),
    ("msc_osteoporosis_differentiation",
     '(TITLE:"mesenchymal stem cells" AND TITLE:"osteoporosis" AND (TITLE:"differentiation" OR TITLE:"adipogenic")) AND SRC:MED',
     ["mesenchymal", "osteoporosis"]),
    ("oa_chondrocyte_scrnaseq",
     '(TITLE:"single-cell RNA" AND TITLE:"osteoarthritis" AND TITLE:"cartilage") AND SRC:MED',
     ["single-cell", "osteoarthritis"]),
    ("chondrocyte_hypertrophy_oa",
     '(TITLE:"chondrocyte hypertrophy" AND TITLE:"osteoarthritis") AND SRC:MED',
     ["chondrocyte", "osteoarthritis"]),
    ("stouffer_meta_zscore",
     '(TITLE:"combining" AND TITLE:"p-values") AND SRC:MED',
     ["p-values"]),
    ("gwas_catalog",
     'TITLE:"The NHGRI-EBI GWAS Catalog" AND SRC:MED',
     ["GWAS Catalog"]),
    ("batch_effect_combat",
     'TITLE:"Adjusting batch effects in microarray expression data using empirical Bayes methods" AND SRC:MED',
     ["batch effects"]),
    ("benjamini_hochberg",
     'TITLE:"Controlling the false discovery rate" AND SRC:MED',
     ["false discovery rate"]),
    ("bone_cartilage_crosstalk",
     '(TITLE:"crosstalk" AND TITLE:"bone" AND TITLE:"cartilage") AND SRC:MED',
     ["bone", "cartilage"]),
    ("spliceosome_bone_osteoblast",
     '((TITLE:"splicing" OR TITLE:"spliceosome") AND (TITLE:"osteoblast" OR TITLE:"bone")) AND SRC:MED',
     ["splic"]),
    ("srsf_chondrocyte",
     '((TITLE:"SRSF" OR TITLE:"serine/arginine") AND (TITLE:"chondrocyte" OR TITLE:"cartilage")) AND SRC:MED',
     ["chondro"]),
    ("nuclear_speckle",
     'TITLE:"nuclear speckles" AND SRC:MED',
     ["nuclear speckle"]),
    ("geo_database",
     'TITLE:"NCBI GEO: archive for functional genomics data sets" AND SRC:MED',
     ["GEO"]),
    ("oa_transcriptome_meta",
     '(TITLE:"osteoarthritis" AND TITLE:"transcriptome" AND (TITLE:"meta-analysis" OR TITLE:"integrative")) AND SRC:MED',
     ["osteoarthritis"]),
    ("op_monocyte_transcriptome",
     '(TITLE:"monocyte" AND TITLE:"osteoporosis" AND (TITLE:"gene expression" OR TITLE:"transcriptom")) AND SRC:MED',
     ["monocyte", "osteoporosis"]),
    ("osteoclast_monocyte_precursor",
     '(TITLE:"osteoclast" AND TITLE:"monocyte" AND TITLE:"precursor") AND SRC:MED',
     ["osteoclast", "monocyte"]),
    ("tissue_heterogeneity_transcriptomics",
     '(TITLE:"tissue" AND TITLE:"heterogeneity" AND TITLE:"gene expression") AND SRC:MED',
     ["heterogeneity"]),
    ("reproducibility_biomarker_signature",
     '(TITLE:"gene signature" AND (TITLE:"reproducib" OR TITLE:"validation")) AND SRC:MED',
     ["signature"]),
]


def epmc_search(query, page_size=5):
    params = {
        "query": query,
        "format": "json",
        "resultType": "core",
        "pageSize": str(page_size),
        "sort": "CITED desc",
    }
    url = EPMC + "?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as f:
                return json.load(f)
        except Exception as e:
            if attempt == 3:
                print("  FAILED:", e)
                return None
            time.sleep(3)
    return None


def pick(res, must_have):
    if not res:
        return None
    items = res.get("resultList", {}).get("result", [])
    for it in items:
        if not it.get("pmid"):
            continue
        title = (it.get("title") or "").lower()
        if all(k.lower() in title for k in must_have):
            return it
    # 放宽：只要有 PMID 就返回首条，但标记为需人工确认
    for it in items:
        if it.get("pmid"):
            it["_loose"] = True
            return it
    return None


rows = []
log = []
for key, q, must in QUERIES:
    print("[QUERY]", key)
    res = epmc_search(q)
    it = pick(res, must)
    if it is None:
        log.append("MISS  %-32s  %s" % (key, q))
        print("   -> NOT FOUND")
        continue
    ji = it.get("journalInfo", {}) or {}
    j = (ji.get("journal", {}) or {}).get("title", "")
    row = {
        "key": key,
        "pmid": it.get("pmid", ""),
        "doi": it.get("doi", ""),
        "title": it.get("title", "").rstrip("."),
        "authors": it.get("authorString", ""),
        "journal": j,
        "year": it.get("pubYear", ""),
        "volume": ji.get("volume", ""),
        "issue": ji.get("issue", ""),
        "pages": it.get("pageInfo", ""),
        "cited_by": it.get("citedByCount", ""),
        "loose_match": "YES" if it.get("_loose") else "",
        "query": q,
    }
    rows.append(row)
    flag = " [LOOSE-需人工确认]" if it.get("_loose") else ""
    line = "OK    %-32s  PMID:%-9s %s (%s) %s%s" % (
        key, row["pmid"], row["journal"], row["year"], row["title"][:70], flag)
    log.append(line)
    print("   ->", row["pmid"], "|", row["journal"], row["year"], "|", row["title"][:70], flag)
    time.sleep(0.4)

csv_path = os.path.join(OUT_DIR, "references_verified.csv")
with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

with open(os.path.join(OUT_DIR, "references_verification_log.txt"), "w", encoding="utf-8") as f:
    f.write("Europe PMC 实时检索核查日志\n")
    f.write("生成时间: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
    f.write("共检索 %d 条，成功取回 %d 条真实记录\n\n" % (len(QUERIES), len(rows)))
    f.write("\n".join(log))

print("\n共取回真实文献 %d / %d 条 -> %s" % (len(rows), len(QUERIES), csv_path))
