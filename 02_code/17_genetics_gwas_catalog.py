#!/usr/bin/env python
"""
17_genetics_gwas_catalog.py
---------------------------
Literature-level genetic support from the GWAS Catalog (alt-full associations TSV).

Replaces the earlier false-positive literature scrape. We query the OFFICIAL
GWAS Catalog for our candidate genes -- 49 mesenchymal-axis splicing factors
(03_results/functional/mesenchymal_leading_edge_genes.csv), SON, and the 21
nominal WGCNA-shared genes -- and report only records whose trait is bone/OA
related (osteoarthritis, bone mineral density / eBMD, osteoporosis, fracture).

Outputs (03_results/genetics/):
  gene_trait_hits_bone_OA.csv  -- real GWAS Catalog hits for our genes on bone/OA traits
  gene_all_gwas_hits.csv       -- every GWAS Catalog trait our genes were mapped to (context)
  genetics_summary.txt         -- counts and interpretation notes
"""
import csv, os, re, sys, zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TSV  = os.path.join(ROOT, "01_data", "gwas", "gwas-catalog-download-associations-alt-full.tsv")
ZIP  = os.path.join(ROOT, "01_data", "gwas", "gwas_assoc_full.zip")
LEDGE= os.path.join(ROOT, "03_results", "functional", "mesenchymal_leading_edge_genes.csv")
OUTD = os.path.join(ROOT, "03_results", "genetics")
os.makedirs(OUTD, exist_ok=True)

WGCNA21 = ["HLA-DMB","MMP9","PIK3CG","BIRC5","ADM","VNN2","SEMA6B","SNX10","ARHGAP15",
           "JUN","SCNN1A","ST8SIA4","ALOX5AP","ADAMTS7","TACSTD2","CACNA1I","ATP8B4",
           "HLX","FSTL1","FCGR2A","KLF5"]

# bone/OA trait regex (case-insensitive, on DISEASE/TRAIT + MAPPED_TRAIT)
TRAIT_RE = re.compile(r"osteoarthritis|bone mineral density|osteoporos|fracture", re.I)

def load_genes():
    genes = set()
    with open(LEDGE, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            g = (row.get("gene") or "").strip()
            if g: genes.add(g.upper())
    genes.add("SON")
    genes.update(g.upper() for g in WGCNA21)
    return genes

def split_genes(s):
    if not s: return []
    return [t.strip().upper() for t in re.split(r"[;,/\s]+", s) if t.strip()]

def ensure_tsv():
    """The 701 MB associations TSV is fully derivable from the archived zip, so
    only the zip is kept in 01_data/gwas. Extract on demand."""
    if os.path.exists(TSV):
        return
    if not os.path.exists(ZIP):
        sys.exit(f"Neither TSV nor source zip found.\n  TSV: {TSV}\n  ZIP: {ZIP}\n"
                 "Download it from "
                 "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
                 "gwas-catalog-associations_ontology-annotated-full.zip")
    print(f"TSV absent; extracting from {os.path.basename(ZIP)} ...")
    with zipfile.ZipFile(ZIP) as z:
        member = z.namelist()[0]
        z.extract(member, os.path.dirname(TSV))
    print(f"extracted -> {TSV}")


def main():
    ensure_tsv()
    genes = load_genes()
    print(f"candidate genes tracked: {len(genes)}")

    bone_oa = []          # (gene, trait, mapped_trait, pval, pmid, snp, risk, orb)
    all_hits = {}         # gene -> set(trait)
    n_rows = 0
    with open(TSV, encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        gi, ri = idx["MAPPED_GENE"], idx["REPORTED GENE(S)"]
        di, mi = idx["DISEASE/TRAIT"], idx["MAPPED_TRAIT"]
        pi, pmi = idx["P-VALUE"], idx["PUBMEDID"]
        si, rsi, oi = idx["STRONGEST SNP-RISK ALLELE"], idx["RISK ALLELE FREQUENCY"], idx["OR or BETA"]
        for line in f:
            n_rows += 1
            if n_rows % 1000000 == 0:
                print(f"  scanned {n_rows:,} rows ...", flush=True)
            row = line.rstrip("\n").split("\t")
            if len(row) <= max(gi, ri, di, mi):
                continue
            mapped = row[gi]; reported = row[ri]
            trait = (row[di] + " | " + row[mi]).lower()
            if not TRAIT_RE.search(trait):
                continue
            gset = set(split_genes(mapped)) | set(split_genes(reported))
            hit = gset & genes
            if not hit:
                continue
            pval = row[pi]
            pmid = row[pmi]
            for g in hit:
                bone_oa.append((g, row[di], row[mi], pval, pmid,
                               row[si], row[rsi], row[oi]))
                all_hits.setdefault(g, set()).add(row[di])

    # write bone/OA hits
    bo_path = os.path.join(OUTD, "gene_trait_hits_bone_OA.csv")
    with open(bo_path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene","disease_trait","mapped_trait","p_value","pubmed_id",
                    "strongest_snp_risk_allele","risk_allele_freq","or_or_beta"])
        w.writerows(bone_oa)
    print(f"bone/OA hits written: {len(bone_oa)} rows -> {bo_path}")

    # write all-trait context (count per gene)
    all_path = os.path.join(OUTD, "gene_all_gwas_hits.csv")
    with open(all_path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene","n_distinct_traits","example_traits"])
        for g in sorted(all_hits, key=lambda x: -len(all_hits[x])):
            w.writerow([g, len(all_hits[g]), " ; ".join(sorted(all_hits[g])[:5])])
    print(f"all-trait context written: {len(all_hits)} genes with >=1 GWAS hit -> {all_path}")

    # summary
    genes_with_bone_oa = sorted({r[0] for r in bone_oa})
    summary = []
    summary.append(f"GWAS Catalog genetics support (alt-full, {n_rows:,} association rows scanned)")
    summary.append(f"candidate genes tracked: {len(genes)}")
    summary.append(f"genes with >=1 GWAS hit on bone/OA traits: {len(genes_with_bone_oa)}")
    summary.append("  " + ", ".join(genes_with_bone_oa) if genes_with_bone_oa else "  (none)")
    summary.append(f"total bone/OA hit rows: {len(bone_oa)}")
    summary.append(f"genes with any GWAS hit (any trait): {len(all_hits)}")
    sum_path = os.path.join(OUTD, "genetics_summary.txt")
    with open(sum_path, "w", encoding="utf-8") as f:
        f.write("\n".join(summary) + "\n")
    print("\n".join(summary))
    print(f"\nsummary -> {sum_path}")

if __name__ == "__main__":
    main()
