###############################################################################
# 16_write_paper_final.R
# 取代 13_write_paper.R + 14_combine_paper.R
# 修复：
#   1. 正文引用从裸 PMID 改为按首次出现顺序编号的 [n]
#   2. 表格直接写入主文档（原先用 body_add_docx 生成 altChunk，
#      [Content_Types].xml 缺 docx 声明，Word 可能提示文件损坏）
#   3. 参考文献重复句点、作者串尾部句点
#   4. 修正 WGCNA 21 基因“富集于髓系/炎症”的说法——正式 GO 检验无显著通路
#   5. 补入 49 个间充质轴一致下调剪接因子的新结果
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(officer); library(flextable); library(magrittr)
})
ROOT <- "D:/bioinfo05"; setwd(ROOT)

refs_all <- fread("06_paper/references_verified.csv")

# ---------------------------------------------------------------------------
# 引用系统：两遍渲染
#   第一遍 CITE_MAP = NULL，只记录 key 首次出现顺序
#   第二遍按该顺序编号，正文写 [n]，文末参考文献同序排列
# ---------------------------------------------------------------------------
CITE_ENV <- new.env()
CITE_ENV$order <- character(0)
CITE_MAP <- NULL

cite <- function(...) {
  keys <- unlist(list(...))
  keys <- keys[keys %in% refs_all$key]
  if (!length(keys)) return("")
  for (k in keys) if (!(k %in% CITE_ENV$order)) CITE_ENV$order <- c(CITE_ENV$order, k)
  if (is.null(CITE_MAP)) return("[#]")
  n <- sort(unique(unlist(CITE_MAP[keys])))
  paste0("[", paste(n, collapse = ","), "]")
}

# 拼句子并清理标点前的空格（引用标记后接句点时会产生 " ."）
para <- function(...) {
  s <- paste(Filter(nzchar, c(...)), collapse = " ")
  s <- gsub("\\s+([.,;:)])", "\\1", s)
  s <- gsub("\\(\\s+", "(", s)
  s <- gsub("\\s{2,}", " ", s)
  trimws(s)
}

# ---------------------------------------------------------------------------
# 排版辅助
# ---------------------------------------------------------------------------
FIG_W <- 6.3
fp_cap <- fp_text(font.size = 9, font.family = "Calibri")
fp_cap_b <- fp_text(font.size = 9, bold = TRUE, font.family = "Calibri")
fp_sup   <- fp_text(font.size = 9, font.family = "Calibri", vertical.align = "superscript")
fp_sup_b <- fp_text(font.size = 9, bold = TRUE, font.family = "Calibri", vertical.align = "superscript")

add_fig <- function(doc, path, title, legend) {
  h <- FIG_W * 0.85
  doc %>%
    body_add_img(src = path, width = FIG_W, height = h) %>%
    body_add_fpar(fpar(ftext(title, fp_cap_b), ftext(paste0(" ", legend), fp_cap)))
}

mk_table <- function(doc, dt, title, note = NULL, widths = NULL) {
  ft <- flextable(as.data.frame(dt)) %>%
    theme_booktabs() %>%
    fontsize(size = 8.5, part = "all") %>%
    bold(part = "header") %>%
    align(align = "left", part = "all") %>%
    padding(padding.top = 1.5, padding.bottom = 1.5, part = "all")
  ft <- if (is.null(widths)) set_table_properties(ft, layout = "autofit", width = 1) else width(ft, width = widths)
  doc <- doc %>% body_add_fpar(fpar(ftext(title, fp_cap_b))) %>% flextable::body_add_flextable(ft)
  if (!is.null(note)) doc <- doc %>% body_add_fpar(fpar(ftext(note, fp_cap)))
  doc %>% body_add_par("", style = "Normal")
}

# 数值列统一有效位数，避免出现 0.14200000000000002 这类机器味
fmt_dt <- function(dt, digits = 3) {
  dt <- copy(as.data.table(dt))
  for (j in names(dt)) if (is.numeric(dt[[j]])) {
    v <- dt[[j]]
    dt[[j]] <- ifelse(is.na(v), "",
                      ifelse(abs(v) < 1e-3 & v != 0, formatC(v, format = "e", digits = 2),
                             formatC(round(v, digits), format = "fg", digits = digits + 1)))
  }
  dt
}

t1  <- fread("05_tables/Table1_cohorts_formatted.csv")
t2a <- fread("05_tables/Table2_convergence_controls.csv")
t2b <- fread("05_tables/Table2_family_level_convergence.csv")
t3  <- fread("05_tables/Table3_top_splicing_pathways.csv")
# Table 3: force 2-decimal display for Z columns so -3.4 renders as -3.40
# (consistency with -4.86 / -4.54; avoids a lone single-decimal outlier)
for (zcol in c("Z OA", "Z BM-MSC")) if (zcol %in% names(t3)) set(t3, j = zcol, value = sprintf("%.2f", t3[[zcol]]))
lead <- if (file.exists("03_results/functional/mesenchymal_leading_edge_genes.csv"))
  fread("03_results/functional/mesenchymal_leading_edge_genes.csv") else NULL

###############################################################################
# 文档构建
###############################################################################
build_paper <- function() {

doc <- read_docx() %>%
  body_add_par("Compartment-specific convergence of the RNA splicing programme in osteoarthritis and osteoporosis", style = "heading 1") %>%
  body_add_par("Transcriptomic sharing follows tissue lineage, not disease label", style = "heading 2") %>%
  body_add_fpar(fpar(
    ftext("Wei Yuwei", fp_cap_b), ftext("1", fp_sup_b),
    ftext(", Li Peng", fp_cap_b), ftext("1", fp_sup_b),
    ftext(", Chen Kai", fp_cap_b), ftext("1", fp_sup_b),
    ftext(", Zhang Dong", fp_cap_b), ftext("2", fp_sup_b), ftext("*", fp_sup_b)
  )) %>%
  body_add_fpar(fpar(ftext("1", fp_sup), ftext(" Department of Orthopedics and Traumatology, Sui Xi Hospital of Traditional Chinese Medicine, No. 1 Baiyang Road, Economic Development Zone, Suixi County, Huaibei City, Anhui Province, 235000, P. R. China", fp_cap))) %>%
  body_add_fpar(fpar(ftext("2", fp_sup), ftext(" Orthopedic Center, Chongqing Hospital of Jiangsu Province Hospital, No. 54 Tuowan Branch Road, Gunan Street, Qijiang District, Chongqing 401420, China", fp_cap))) %>%
  body_add_fpar(fpar(ftext("* Corresponding author: Zhang Dong, MD, Chief Physician. Email: zhangdht@126.com. ORCID: 0009-0006-9000-273X", fp_cap)))

# --------------------------------------------------------------------- Abstract
doc <- doc %>% body_add_par("Abstract", style = "heading 1")
doc <- doc %>% body_add_par(para(
  "Background. Osteoarthritis (OA) and osteoporosis (OP) are usually described as opposing diseases of cartilage and bone, yet their molecular relationship is still disputed.",
  "Most transcriptomic comparisons pool heterogeneous tissues under a single disease label. We asked whether that convention hides a real signal."
), style = "Normal")
doc <- doc %>% body_add_par(para(
  "Methods. Seven publicly available microarray and RNA-sequencing cohorts (303 samples) covering OA cartilage, OP circulating monocytes and OP bone-marrow mesenchymal stem cells (BM-MSCs) were harmonised to a common gene universe.",
  "We performed differential expression, weighted gene co-expression network analysis, gene-set variation analysis (GSVA) and Stouffer meta-analysis at pathway level.",
  "Convergence was assessed against three complementary null models: pathway-label permutation conditional on marginal significance, within-cohort case/control label permutation, and gene-set family clustering to remove pathway redundancy.",
  "Findings were examined in an independent single-cell atlas of human knee cartilage."
), style = "Normal")
doc <- doc %>% body_add_par(para(
  "Results. At gene level OA and OP shared essentially nothing: 10 overlapping genes out of 10,373 tested, fewer than expected by chance, and a genome-wide meta z-score correlation of -0.020.",
  "No cross-disease co-expression module pair survived correction (maximum Jaccard 0.032; all FDR > 0.99).",
  "Pooling all OP tissues also showed no pathway convergence (enrichment 0.62, permutation P = 0.93).",
  "Splitting OP by tissue compartment changed the picture. OA cartilage and OP BM-MSCs, both mesenchymal, converged significantly (family-level enrichment 1.20, permutation P = 0.004), whereas OA cartilage and OP monocytes fell below the null expectation (enrichment 0.55, P = 0.998).",
  "RNA splicing and mRNA metabolism drove the mesenchymal convergence: 5 of 66 convergent gene-set families were splicing-related against a background of 12 of 2,168 (13.7-fold, P = 1.5e-05), and all were down-regulated.",
  "Forty-nine core splicing factors, including SF3B1, HNRNPU, SRRM1, PRPF8 and CDC5L, were concordantly down-regulated in both mesenchymal tissues.",
  "SON, a spliceosome scaffold recently nominated as the causal hub of this comorbidity and screened as a drug target, tracked the programme but did not explain it: removing SON left the correlations unchanged, identifying it as a reporter of the module rather than its controller."
), style = "Normal")
doc <- doc %>% body_add_par(para(
  "Conclusions. The transcriptomic link between OA and OP is compartment-specific rather than disease-specific.",
  "A shared down-regulation of RNA splicing machinery unites articular cartilage and bone-marrow stroma, while the haematopoietic compartment follows a separate and partly opposite trajectory.",
  "Pooling OP tissues cancels this signal, which may explain why previous searches for a shared OA-OP hub have been inconsistent."
), style = "Normal")
doc <- doc %>% body_add_par("Keywords: osteoarthritis; osteoporosis; RNA splicing; tissue compartment; transcriptomics; mesenchymal stem cells", style = "Normal")

# ----------------------------------------------------------------- Introduction
doc <- doc %>% body_add_par("Introduction", style = "heading 1")
intro <- list(
  para("Osteoarthritis (OA) and osteoporosis (OP) are the two most common age-related musculoskeletal disorders and together account for a large share of chronic disability and health-care expenditure", cite("op_oa_comorbidity"), "."),
  para("Clinically the two have long been considered inversely related: high bone mineral density (BMD) is associated with an increased risk of knee and hip OA, whereas low BMD predicts fracture but appears protective against OA", cite("hartley_bmd_oa_mr"), ". This bone-cartilage paradox has recently been traced to subchondral bone compliance and Wnt-axis signalling", cite("peng2026_bone_cartilage_paradox"), "."),
  para("Mendelian randomisation studies support a causal component to this inverse relationship, indicating that genetic liability to higher BMD increases OA risk", cite("morris2019_ebmd", "tachmazidou2019_oa"), "."),
  para("What connects cartilage degradation to bone-marrow failure at the transcriptomic level remains unclear. Several integrative studies have searched for a shared molecular hub, and each has returned a different answer. An analysis of mesenchymal stem cells from both diseases proposed ten hub genes centred on matrix and differentiation biology, among them COL9A3, MMP3 and PTH1R", cite("chen2024_bmsc_hub"), ". A multimodal study intersecting OA and OP differential expression implicated ferroptosis regulators such as TXNIP and SLC2A3", cite("wang2025_ferroptosis"), ". A recent multi-omics effort combining co-expression networks, Mendelian randomisation and single-cell data nominated the nuclear speckle scaffold SON as the key protein linking the two diseases, and carried it forward into drug screening", cite("jiang2026_son_target"), "."),
  para("These candidate hubs do not overlap with one another. Such divergence is usually attributed to differences in cohort or platform, but there is a simpler possibility that has not been tested directly: the studies sample different tissues."),
  para("A recurring obstacle is tissue heterogeneity. OA studies almost always profile articular chondrocytes, whereas OP studies have used circulating monocytes, bone-marrow mesenchymal stem cells (BM-MSCs), osteoblasts or whole bone biopsies", cite("benisch2012_gse35958"), "."),
  para("These cell types belong to different developmental compartments. Chondrocytes and BM-MSCs derive from mesenchyme; monocytes and osteoclast precursors belong to the haematopoietic lineage", cite("osteoclast_monocyte_precursor"), "."),
  para("Pooling them under one disease label assumes that the disease, not the tissue, determines the transcriptomic programme. If the opposite is true, averaging opposing tissue effects will mask genuine biology rather than reveal it."),
  para("RNA splicing is a plausible candidate for a lineage-level shared programme. Splicing-factor expression declines with replicative and stress-induced senescence, and restoring it can partially reverse senescent phenotypes", cite("splicing_ageing_harries", "splicing_factor_senescence"), "."),
  para("In bone, commitment of BM-MSCs to the osteogenic rather than adipogenic fate depends on precise alternative splicing of lineage transcripts", cite("spliceosome_bone_osteoblast"), "."),
  para("We therefore hypothesised that OA cartilage and OP BM-MSCs, both mesenchymal derivatives, share a splicing programme that becomes invisible once OP monocytes are folded into the same comparison."),
  para("Here we test this compartment-specific convergence hypothesis across seven cohorts. We then return to SON as a worked example. Because it has already been nominated as the causal hub of this comorbidity and advanced as a druggable target", cite("jiang2026_son_target"), ", it offers an unusually concrete test of whether hub-centric reasoning survives a pathway-level analysis of the same biology.")
)
for (p in intro) doc <- doc %>% body_add_par(p, style = "Normal")

# ---------------------------------------------------------------------- Results
doc <- doc %>% body_add_par("Results", style = "heading 1")

doc <- doc %>% body_add_par("Seven cohorts spanning three tissue compartments", style = "heading 2")
doc <- doc %>% body_add_par(para("We assembled seven transcriptomic cohorts comprising 303 samples: four OA cartilage cohorts (GSE114007, GSE57218, GSE117999, GSE169077)", cite("fisch2018_gse114007"), ", two OP circulating-monocyte cohorts (GSE56815, GSE7158) and one OP BM-MSC cohort (GSE35958)", cite("benisch2012_gse35958"), "."), style = "Normal")
doc <- doc %>% body_add_par(para("Cohort size ranged from 9 to 80 samples across RNA sequencing and Affymetrix microarray platforms (Table 1). After batch correction and within-cohort standardisation, 10,373 genes were quantified in every cohort and used as the common universe for all downstream analyses. Figure 1 summarises the study design and cohort landscape."), style = "Normal")
doc <- add_fig(doc, "04_figures/Fig1_design_and_landscape.png",
  "Figure 1. Study design and data landscape.",
  para("(a) Compartment model: transcriptomic sharing is hypothesised to follow tissue lineage (mesenchymal versus haematopoietic) rather than disease label.",
       "(b) The seven cohorts and their tissue assignment.",
       "(c) Pairwise Spearman correlation of per-cohort pathway z-scores; individual cohorts correlate only weakly, motivating meta-analysis rather than direct pooling."))

doc <- doc %>% body_add_par("OA and OP share no gene-level signature", style = "heading 2")
doc <- doc %>% body_add_par(para("Meta-analysis across the four OA cartilage cohorts identified 1,124 differentially expressed genes at FDR < 0.05; the three OP cohorts yielded 231 (Figure 2)."), style = "Normal")
doc <- doc %>% body_add_par(para("Only ten genes appeared in both lists, fewer than the number expected under independence (hypergeometric P = 0.67), and the genome-wide correlation of meta z-scores was -0.020 (Spearman P = 0.16)."), style = "Normal")
doc <- doc %>% body_add_par(para("Weighted gene co-expression network analysis", cite("langfelder2008_wgcna"), "identified disease-associated modules within each condition, but no cross-disease module pair overlapped more than expected by chance (maximum Jaccard 0.032; all FDR > 0.99)."), style = "Normal")
doc <- doc %>% body_add_par(para("Twenty-one genes were nominally shared between at least one OA module and one OP module. Several are myeloid or inflammatory (HLA-DMB, FCGR2A, ALOX5AP, MMP9, VNN2), but a formal Gene Ontology test against the 10,373-gene universe returned no enriched biological process at FDR < 0.05. We therefore treat this list as descriptive only and do not interpret it as a shared inflammatory programme."), style = "Normal")
doc <- add_fig(doc, "04_figures/Fig2_gene_level_no_overlap.png",
  "Figure 2. Gene-level analyses reveal no shared OA-OP signature.",
  para("(a,b) Volcano plots of meta-analytic differential expression for OA cartilage (a) and OP (b).",
       "(c) OA versus OP meta z-scores; only ten genes reached FDR < 0.05 in both diseases.",
       "(d) Overlap of differentially expressed genes.",
       "(e) WGCNA cross-disease module correspondence; no module pair survived correction."))

doc <- doc %>% body_add_par("Pathway convergence appears only after splitting OP by compartment", style = "heading 2")
doc <- doc %>% body_add_par(para("We next moved to pathway level using GSVA", cite("hanzelmann2013_gsva"), "and Stouffer-weighted meta-analysis across 5,059 filtered gene sets from Hallmark, KEGG, Reactome and Gene Ontology biological process."), style = "Normal")
doc <- doc %>% body_add_par(para("Treated as two pooled diseases, OA and OP did not converge: six pathways were concordantly significant against a permutation mean of 9.65 (enrichment 0.62, P = 0.93), with a pathway-level correlation of -0.020."), style = "Normal")
doc <- doc %>% body_add_par(para("Before concluding that no convergence exists, we verified that the pipeline can detect one. Two OA cohorts profiled on different platforms converged strongly (enrichment 3.06, P < 0.001, rho = 0.328), confirming adequate sensitivity (Table 2)."), style = "Normal")
doc <- doc %>% body_add_par(para("Splitting OP by tissue then revealed a signal that pooling had erased. OA cartilage and OP BM-MSCs converged with an enrichment of 3.33 (P = 0.030, rho = 0.142) at pathway level, comparable in magnitude to the cross-platform OA positive control."), style = "Normal")
doc <- doc %>% body_add_par(para("The opposite pairing behaved in the opposite way: OA cartilage versus OP monocytes fell below the null expectation (enrichment 0.44, P = 0.997, rho = -0.065). The two compartments are not merely unrelated; they trend against each other, which is precisely why averaging them produces a null result."), style = "Normal")
doc <- doc %>% body_add_par(para("Because pathway databases contain heavily overlapping gene sets, a hypergeometric test over individual pathways overstates the evidence. We collapsed the 5,059 pathways into 2,168 gene-set families by Jaccard hierarchical clustering and repeated the analysis on family representatives. Mesenchymal convergence survived (66 observed versus 55.1 expected families, permutation P = 0.004), the pooled comparison did not (P = 0.93), and the OA-monocyte comparison remained below null (P = 0.998). Figure 3 presents these comparisons."), style = "Normal")
doc <- add_fig(doc, "04_figures/Fig3_compartment_convergence.png",
  "Figure 3. Pathway-level convergence is compartment-specific.",
  para("(a) Convergence enrichment (observed over permutation-expected) for seven comparisons, including positive and negative controls.",
       "(b,c) Pathway z-scores for OA cartilage versus OP BM-MSC (b) and versus OP monocytes (c).",
       "(d) Family-level convergence after collapsing redundant gene sets; only the mesenchymal pair remains significant."))

doc <- doc %>% body_add_par("RNA splicing is the shared mesenchymal axis", style = "heading 2")
doc <- doc %>% body_add_par(para("We then asked which biology drives mesenchymal convergence. Of the 87 pathways concordant between OA cartilage and OP BM-MSCs, 13 were RNA splicing or mRNA metabolism terms, a 17.1-fold enrichment over background (hypergeometric P = 1.3e-10), and every one was down-regulated in disease."), style = "Normal")
doc <- doc %>% body_add_par(para("At family level, where redundancy is controlled, 5 of the 66 convergent families were splicing-related against a background of 12 among 2,168 families (13.7-fold, P = 1.5e-05). Representative families were REGULATION OF RNA SPLICING (z = -4.86 in OA, -2.94 in BM-MSC), MRNA SPLICE SITE SELECTION (-4.68, -3.05) and MRNA TRANSPORT (-3.37, -2.20) (Table 3; Figure 4)."), style = "Normal")
doc <- doc %>% body_add_par(para("The enrichment was specific to the mesenchymal axis. No splicing pathway appeared among the pathways concordant between OA cartilage and OP monocytes, nor among the 128 pathways significant in monocytes alone."), style = "Normal")
if (!is.null(lead) && nrow(lead) > 0) {
doc <- doc %>% body_add_par(para(sprintf("Resolving the axis to individual transcripts, %d core splicing factors were down-regulated in both OA cartilage and OP BM-MSCs, including the branch-point recognition factor SF3B1, the heterogeneous nuclear ribonucleoprotein HNRNPU, the serine/arginine repetitive matrix protein SRRM1, the tri-snRNP component PRPF8 and the spliceosome-activating factor CDC5L (Table 4).", nrow(lead)),
  "Gene Ontology testing of these 49 genes returned 68 enriched biological processes, all centred on mRNA splicing and processing, confirming that the signal is not an artefact of gene-set annotation breadth."), style = "Normal")
doc <- doc %>% body_add_par(para("In OP monocytes the same 49 genes showed no coordinated behaviour. Across the two monocyte cohorts only 11 of 49 (22%) were concordantly down-regulated, indistinguishable from chance (binomial P = 0.74), and neither cohort alone reached significance (GSE56815: 22/49 down, P = 0.57; GSE7158: 27/49 down, P = 0.57). The mesenchymal down-regulation is therefore not a general disease or ageing effect visible in any tissue."), style = "Normal")
}
doc <- doc %>% body_add_par(para("Sample-level scores told the same story. A splicing programme score was reduced in OA cartilage (meta z = -2.53, P = 0.011) and strongly reduced in OP BM-MSCs (Hedges' g = -1.01), though the latter rests on nine samples and must be read as an effect-size estimate rather than a confirmatory test."), style = "Normal")
doc <- doc %>% body_add_par(para("OP monocytes moved in the opposite direction (GSE56815 g = +0.53, P = 0.019). Control programmes for osteoblast, chondrocyte and ribosomal biology showed no consistent cross-compartment convergence, indicating that the splicing result is not a generic transcriptional-activity artefact."), style = "Normal")
doc <- add_fig(doc, "04_figures/Fig4_splicing_axis.png",
  "Figure 4. RNA splicing is the mesenchymal shared axis.",
  para("(a) Splicing-pathway enrichment across compartment-defined gene sets; only the mesenchymal axis is enriched.",
       "(b) Top splicing pathways in the mesenchymal axis, all down-regulated in both tissues.",
       "(c) Splicing programme score by cohort (Hedges' g, case versus control).",
       "(d) Meta-analytic scores for four programmes; only splicing is consistently down-regulated in the mesenchymal compartment."))

doc <- doc %>% body_add_par("Single-cell data do not localise the splicing shift to a chondrocyte subtype", style = "heading 2")
doc <- doc %>% body_add_par(para("We examined an independent single-cell atlas of human knee cartilage (GSE324993; 37,453 cells from three healthy and three OA donors)", cite("oa_chondrocyte_scrnaseq"), "to ask whether the splicing shift is confined to a particular chondrocyte state."), style = "Normal")
doc <- doc %>% body_add_par(para("After pseudobulk aggregation by donor and subtype, which avoids the pseudoreplication that inflates single-cell P-values, the splicing programme was negative in six of eight subtypes but reached significance in none (all P > 0.45)."), style = "Normal")
doc <- doc %>% body_add_par(para("Cell-type composition did not differ between healthy and OA cartilage, and composition-adjusted bulk scores remained non-significant (P = 0.70). Effect sizes for the chondrocyte, osteoblast and hypertrophy control programmes were numerically larger than for splicing."), style = "Normal")
doc <- doc %>% body_add_par(para("We read this as a negative result at the resolution available. With three donors per group the analysis has little power, and a diffuse shift spread across all subtypes is exactly the pattern least likely to survive subtype-wise testing. SON was absent from the filtered single-cell matrix, so cell-type-resolved SON analysis was not possible. Figure 5 shows the single-cell landscape."), style = "Normal")
doc <- add_fig(doc, "04_figures/Fig5_single_cell_landscape.png",
  "Figure 5. Single-cell cartilage atlas shows no subtype-specific splicing shift.",
  para("(a,b) UMAP of 37,453 chondrocytes from three healthy and three OA donors, coloured by subtype (a) and condition (b).",
       "(c) Cell-type composition; no subtype differs significantly between groups.",
       "(d) Splicing programme score by subtype after pseudobulk aggregation; all FDR > 0.45."))

doc <- doc %>% body_add_par("SON reports the splicing programme but does not drive it", style = "heading 2")
doc <- doc %>% body_add_par(para("SON organises nuclear speckles and is required for efficient spliceosome assembly", cite("son_splicing_function"), ", which made it an attractive candidate hub before the analysis began. We tested that expectation explicitly."), style = "Normal")
doc <- doc %>% body_add_par(para("SON was down-regulated in OA cartilage (GSE114007 logFC = -0.40, P = 0.046) and in OP BM-MSCs (logFC = -0.38, P = 0.015). In the two OP monocyte cohorts it behaved inconsistently: significantly down in GSE7158 (logFC = -0.35, P = 0.012) but unchanged in GSE56815 (logFC = 0.02, P = 0.65), and the coordinated splicing programme was recovered in neither monocyte cohort (only 11 of 49 leading-edge factors concordantly down, P = 0.74). SON expression correlated with the splicing programme score in five of seven cohorts."), style = "Normal")
doc <- doc %>% body_add_par(para("That correlation is partly circular, because SON belongs to the splicing gene set. Removing SON from the set left the correlations essentially unchanged (for example 0.590 to 0.577), and benchmarking against size-matched random gene sets showed that SON exceeded the random background in only four of seven cohorts."), style = "Normal")
doc <- doc %>% body_add_par(para("SON therefore behaves as a faithful reporter of a coordinated splicing module rather than its upstream controller."), style = "Normal")
doc <- doc %>% body_add_par(para("This bears directly on the recent proposal that SON is the key protein linking OP and OA", cite("jiang2026_son_target"), ". Our data reproduce the observation on which that proposal rests: SON is genuinely down-regulated, and in the same direction. What the compartment analysis adds is context. SON is one of 49 splicing factors moving together in mesenchymal tissue, it carries no information beyond that module once the module is accounted for, and the co-ordination is restricted to the mesenchymal compartment. The earlier study evaluated the OP side almost entirely in peripheral blood and circulating monocytes (GSE56815, GSE48556, GSE7429, GSE7158), which is precisely the compartment in which we find the programme is not shared, and it took cartilage (GSE114007) as the OA side. A single-node interpretation is therefore difficult to reconcile with the multi-gene, compartment-restricted structure of the signal, even though the individual measurement is not in dispute."), style = "Normal")
doc <- doc %>% body_add_par(para("As an out-of-sample test we trained LASSO, random forest and linear-kernel SVM classifiers on OA cartilage using the splicing genes. Within OA the classifiers were moderately accurate (AUC 0.70-0.73), but transfer failed: AUC 0.30 to OP BM-MSCs (n = 9), near 0.50 to OP monocytes, and 0.32-0.60 in the reverse direction."), style = "Normal")
doc <- doc %>% body_add_par(para("Given nine samples in the BM-MSC cohort, the empirical null distribution of transfer AUC had a standard deviation of 0.245, so this test has almost no power and cannot distinguish a genuinely absent signature from an undetectable one. The honest conclusion is that population-level pathway convergence does not by itself deliver an individual-level biomarker. Figure 6 summarises the SON and classifier analyses."), style = "Normal")
doc <- add_fig(doc, "04_figures/Fig6_SON_ML_genetics.png",
  "Figure 6. SON tracks the splicing programme; classifier transfer fails.",
  para("(a) SON differential expression across cohorts.",
       "(b) Correlation between SON expression and splicing programme score, before and after removing SON from the gene set.",
       "(c) Cross-compartment classifier transfer AUCs with empirical null distributions.",
       "(d) Genetic convergence: GWAS Catalog loci for bone and joint traits overlap the candidate gene set (see Discussion)."))

# ------------------------------------------------------------------- Discussion
doc <- doc %>% body_add_par("Discussion", style = "heading 1")
disc <- list(
  para("Our analysis reframes a long-standing question. Rather than asking whether OA and OP share a gene signature, we asked whether they share a programme once their constituent tissues are assigned to the correct developmental compartment."),
  para("The answer has three layers. At gene level there is no overlap. At pathway level there is none either, as long as OP tissues are pooled. Within the mesenchymal compartment, however, OA cartilage and OP BM-MSCs share a reproducible down-regulation of RNA splicing and mRNA metabolism machinery."),
  para("This structure reconciles a literature that has often looked contradictory. Comparisons of OA cartilage with OP whole blood or monocytes see little sharing because the haematopoietic compartment carries an independent, partly opposite signal", cite("jiang2026_son_target"), ". Studies that happened to use bone-marrow stroma may have detected mesenchymal changes without recognising them as a lineage-level rather than disease-level phenomenon", cite("chen2024_bmsc_hub"), "."),
  para("It also explains why our gene-level result reads as more negative than published ones. Intersecting nominally significant differentially expressed genes from two cohorts routinely recovers tens of shared genes and a compact hub list", cite("chen2024_bmsc_hub", "wang2025_ferroptosis"), ". We recover such genes as well: 21 appear in co-expression modules of both diseases. The difference lies in what follows. Formal over-representation testing of those 21 genes against the 10,373-gene universe returned no enriched biological process at FDR < 0.05, and they did not reproduce in cross-cohort meta-analysis. Intersections of this size are expected under the null when cohorts are small and thresholds are nominal, which is consistent with the observation that successive studies keep reporting hub lists that do not overlap each other."),
  para("The biology is coherent with what is known about mesenchymal ageing. Splicing-factor expression falls with senescence", cite("splicing_ageing_harries", "splicing_factor_senescence"), "and MSC lineage choice depends on alternative splicing of osteogenic transcripts", cite("spliceosome_bone_osteoblast"), ". A shared decline in splicing capacity across cartilage and marrow stroma is therefore more plausibly a common mesenchymal ageing or stress response than a disease-transmission mechanism between joint and bone."),
  para("The SON analysis shows the practical cost of candidate-hub reasoning. SON is mechanistically attractive, it is down-regulated in the right tissues, and it correlates with the programme; on that basis it has already been proposed as the molecular link between OP and OA and taken into target screening", cite("jiang2026_son_target"), ". We do not dispute the underlying measurement, and our own data agree with it. The disagreement concerns what the measurement licenses. Removing SON from the gene set barely changed the correlations, which is the signature of a passenger rather than a driver, and the module it reports is coherent only within the mesenchymal compartment. Single-node intervention against a distributed module of this kind is unlikely to restore it, and a target nominated from haematopoietic or pooled tissue may not act in the compartment where the programme actually holds together. The wider point is that co-expression hubs are selected partly for connectivity, so a gene can rank first in a network and still be the least informative member of its own module."),
  para("Two results were negative and we report them as such. The single-cell atlas did not localise the splicing shift to any chondrocyte subtype, and cross-compartment classifiers did not transfer. Both analyses are underpowered, with three donors per group and nine BM-MSC samples respectively, so they constrain interpretation rather than refute the main finding."),
  para("Independent genetic evidence, although observational, converges with the transcriptomic result. Querying the GWAS Catalog (www.ebi.ac.uk/gwas) for the 70 candidate genes against bone and joint traits returned 114 genome-wide significant association records at 17 loci, including SON and ten additional core splicing factors (CDC5L, SF3A3, RBM39, SART1, QKI, PRPF38B, DHX15, CTNNBL1, KHDRBS1, NCBP1); the strongest signals were with estimated BMD at the heel, where the CDC5L locus reached P < 1e-100. These are tract-wide associations rather than proof of causal splicing dysregulation and do not establish direction, but they place the mesenchymal splicing programme within the inherited architecture of bone and joint liability."),
  para("Two recent genetic studies of splicing in joint tissue support both the plausibility of the mechanism and, independently, its tissue specificity. Splicing quantitative trait loci mapped across cartilage and synovium show a tissue-specific regulatory architecture for joint-related traits", cite("tian2026_cartilage_sqtl"), ", and response splicing quantitative trait loci in primary human chondrocytes identify OA risk genes that act through isoform usage rather than total transcript abundance", cite("byun2025_chondrocyte_sqtl"), ". Together these establish that splicing regulation in joint tissue is heritable, disease-relevant and not uniform across tissues, which is the pattern our transcriptomic result would predict. Neither study profiles bone-marrow stroma, so the mesenchymal-axis extension we propose here remains an open prediction rather than a confirmed one."),
  para("The remaining gap is causal direction. We could not complete a Mendelian randomisation analysis: full GWAS summary statistics for estimated BMD, fracture and OA were not retrievable at a workable download rate within this project, and partial files cannot support valid instrument selection. We therefore cannot say whether reduced splicing capacity precedes disease, follows it, or reflects shared upstream ageing. Extending the existing cartilage and synovium splicing-QTL maps to bone-marrow stroma, and combining them with the large BMD and OA GWAS", cite("morris2019_ebmd", "tachmazidou2019_oa"), ", would be the natural way to settle this."),
  para("Other limitations deserve emphasis. Only one BM-MSC cohort exists in the public domain, and it contains nine samples, so the mesenchymal axis rests asymmetrically on a small dataset even though the pathway-level statistics account for it. Cohorts differ in platform and in the anatomical site sampled. Sex and age distributions could not be harmonised across all datasets. Finally, transcript abundance of splicing factors is an indirect proxy: we did not measure splicing events themselves, and confirming altered isoform usage would require junction-level RNA-sequencing analysis in matched tissues."),
  para("In summary, the transcriptomic convergence between OA and OP is real but compartment-specific. Articular cartilage and bone-marrow stroma are linked by a down-regulated RNA splicing programme, while the haematopoietic compartment follows its own trajectory. Study designs that pool tissues by disease label will keep cancelling this signal, and biomarker or therapeutic work should be anchored to the tissue compartment rather than the diagnosis.")
)
for (p in disc) doc <- doc %>% body_add_par(p, style = "Normal")

# ---------------------------------------------------------------------- Methods
doc <- doc %>% body_add_par("Methods", style = "heading 1")
methods <- list(
  para("Data collection and preprocessing. Seven publicly available GEO datasets were retrieved and processed with a common pipeline. Expression matrices were log2-transformed, corrected for batch effects with ComBat", cite("batch_effect_combat"), ", and filtered to remove genes with zero variance or missing values. Case/control status, donor identifiers and available clinical metadata were harmonised across cohorts. Analyses were restricted to the 10,373 genes measured in all seven cohorts."),
  para("Differential expression. limma-voom was used for RNA-sequencing data and limma for microarrays", cite("ritchie2015_limma"), ". Paired designs were modelled with duplicateCorrelation using donor as the blocking variable. Cohort-level statistics were combined by Stouffer's weighted z-score method with weights proportional to the square root of sample size, and multiplicity was controlled by the Benjamini-Hochberg procedure", cite("benjamini1995_bh"), "."),
  para("Co-expression network analysis. WGCNA", cite("langfelder2008_wgcna"), "was run separately per disease with a soft-threshold power chosen by scale-free topology fit. Module-trait associations were computed against case/control status, and cross-disease module correspondence was tested by hypergeometric overlap with Benjamini-Hochberg correction."),
  para("Pathway analysis. GSVA scores were computed with Hallmark, KEGG, Reactome and Gene Ontology biological process gene sets from MSigDB, retaining sets with 15 to 300 genes present in the common universe (5,059 sets). Differential pathway activity was assessed with limma and combined across cohorts by the same Stouffer procedure. Convergence was defined as concordant direction with FDR < 0.05 in both arms of a comparison."),
  para("Null models. Three complementary nulls were used. First, pathway labels were permuted conditional on the marginal significant sets of each disease (2,000 permutations), which tests convergence while preserving the observed number of significant pathways. Second, case/control labels were permuted within cohorts (500 permutations) and the entire pipeline was rerun; this preserves the correlation structure among pathways but destroys all disease signal, so it establishes that signal exists rather than discriminating between competing convergence models. Third, gene sets were clustered into families by Jaccard distance (hierarchical clustering, cut height 0.7, yielding 2,168 families) and the analysis was repeated on family representatives to remove the pseudo-independence created by shared genes."),
  para("Functional annotation. Gene Ontology biological process over-representation was tested with clusterProfiler using the 10,373-gene common universe as background rather than the whole genome, and Benjamini-Hochberg correction at FDR < 0.05."),
  para("Single-cell analysis. The GSE324993 dataset was processed in Seurat", cite("hao2021_seurat"), "with log-normalisation, principal component analysis and UMAP. Programme scores were computed with AddModuleScore against control gene sets of matched expression. To avoid pseudoreplication, cell-level scores were aggregated to donor-by-subtype pseudobulk values before testing with two-sided t-tests and Hedges' g. Composition effects were assessed by comparing subtype proportions between groups and by recomputing bulk scores weighted to a common composition."),
  para("Machine learning. LASSO", cite("friedman2010_glmnet"), ", random forest and linear-kernel SVM classifiers were trained on within-cohort z-standardised expression of the splicing gene set in one compartment and evaluated in another. Discrimination was summarised by area under the ROC curve", cite("robin2011_proc"), ". Empirical null distributions were generated from 200 size-matched random gene sets."),
  para("Genetic evidence. The GWAS Catalog (www.ebi.ac.uk/gwas, alt-full associations release) was queried for the 70 candidate genes (49 mesenchymal-axis splicing factors, including SON, plus 21 nominal WGCNA-shared genes) against bone and joint traits (osteoarthritis, bone mineral density, osteoporosis, fracture). A record was retained when its mapped or reported gene matched a candidate and its trait matched a bone or joint term."),
  para("Software and reproducibility. All analyses were performed in R 4.6.1 with a fixed random seed. Analysis code, intermediate objects, figures and tables are deposited in the project repository, organised so that each numbered script can be rerun independently.")
)
for (p in methods) doc <- doc %>% body_add_par(p, style = "Normal")

# ----------------------------------------------------------------- Declarations
doc <- doc %>% body_add_par("Declarations", style = "heading 1")
doc <- doc %>% body_add_par(para("Ethics approval. Not applicable. This study analysed only publicly available, de-identified datasets obtained from the Gene Expression Omnibus. Each original study obtained its own ethical approval and informed consent."), style = "Normal")
doc <- doc %>% body_add_par(para("Data availability. All datasets are publicly available from the Gene Expression Omnibus under accessions GSE114007, GSE57218, GSE117999, GSE169077, GSE56815, GSE7158, GSE35958 and GSE324993."), style = "Normal")
doc <- doc %>% body_add_par(para("Code availability. Analysis code and intermediate results are deposited in a public GitHub repository (URL: https://github.com/[ACCOUNT]/[REPO]) and archived at Zenodo (DOI: 10.5281/zenodo.[ID]); the repository is set to public and released under an open-source licence."), style = "Normal")
doc <- doc %>% body_add_par(para("Funding. This work was supported by the Self-initiated Science and Technology Plan Project of Huaibei Municipal Bureau of Science and Technology (No. 2025HK044)."), style = "Normal")
doc <- doc %>% body_add_par(para("Author contributions. Wei Yuwei: data curation, formal analysis, investigation, visualisation, writing — review and editing. Li Peng: methodology, software, validation, writing — review and editing. Chen Kai: resources, data curation, validation, writing — review and editing. Zhang Dong: conceptualization, methodology, formal analysis, writing — original draft, supervision, funding acquisition. All authors read and approved the final manuscript."), style = "Normal")
doc <- doc %>% body_add_par(para("Corresponding author. Zhang Dong, MD, Chief Physician. Email: zhangdht@126.com. ORCID: 0009-0006-9000-273X."), style = "Normal")
doc <- doc %>% body_add_par(para("Competing interests. The authors declare no competing interests."), style = "Normal")

# --------------------------------------------------------------------- Tables
doc <- doc %>% body_add_break() %>% body_add_par("Tables", style = "heading 1")

doc <- mk_table(doc, t1, "Table 1. Transcriptomic cohorts included in the analysis.",
  "BM-MSC, bone-marrow mesenchymal stem cell. Gene counts refer to genes retained after per-cohort filtering; all cross-cohort analyses used the 10,373 genes common to every dataset.")

doc <- mk_table(doc, fmt_dt(t2a), "Table 2a. Pathway-level convergence and control comparisons.",
  "Enrichment is the observed number of concordant significant pathways divided by the permutation-derived expectation. P values come from 2,000 pathway-label permutations conditional on marginal significance.")

doc <- mk_table(doc, fmt_dt(t2b), "Table 2b. Family-level convergence after collapsing redundant gene sets.",
  "The 5,059 pathways were clustered into 2,168 families by Jaccard distance (cut height 0.7). Analysis was repeated on family representatives to remove pseudo-independence from shared genes.")

doc <- mk_table(doc, fmt_dt(t3), "Table 3. Splicing pathways of the mesenchymal axis.",
  "Z-scores are Stouffer-weighted meta-analytic values for OA cartilage and limma moderated statistics for the OP BM-MSC cohort. All listed pathways are down-regulated in both tissues.")

if (!is.null(lead) && nrow(lead) > 0) {
  lt <- copy(lead)
  setnames(lt, old = intersect(names(lt), c("gene","z_oa","t_msc","p_msc","t_mono")),
           new = c("Gene","Z (OA cartilage)","t (OP BM-MSC)","P (OP BM-MSC)","t (OP monocyte)")[
             seq_along(intersect(names(lt), c("gene","z_oa","t_msc","p_msc","t_mono")))])
  doc <- mk_table(doc, fmt_dt(lt), 
    sprintf("Table 4. Core splicing factors down-regulated in both OA cartilage and OP BM-MSCs (n = %d).", nrow(lt)),
    "Genes are ordered by combined effect. The final column shows the same genes in OP circulating monocytes, where 11 of 49 (22%) were concordantly down-regulated across both cohorts (GSE56815: 22/49 down; GSE7158: 27/49 down); none of these departures differed from chance (binomial P = 0.74 for concordance).")
}

# ------------------------------------------------------------------ References
doc <- doc %>% body_add_break() %>% body_add_par("References", style = "heading 1")
if (!is.null(CITE_MAP)) {
  ord <- CITE_ENV$order
  used <- refs_all[match(ord, key)]
  used <- used[!is.na(key)]
  for (i in seq_len(nrow(used))) {
    r <- used[i]
    au <- sub("\\.\\s*$", "", trimws(r$authors))
    ti <- sub("\\.\\s*$", "", trimws(r$title))
    jr <- sub("\\.\\s*$", "", trimws(r$journal))
    id <- if (!is.na(r$doi) && nzchar(r$doi)) paste0("doi:", r$doi) else paste0("PMID:", r$pmid)
    txt <- sprintf("%d. %s. %s. %s. %s. %s.", i, au, ti, jr, r$year, id)
    doc <- doc %>% body_add_par(txt, style = "Normal")
  }
  # 未被正文引用但已核实的文献单列，供投稿时补充
  unused <- refs_all[!key %in% ord]
  if (nrow(unused) > 0) {
    doc <- doc %>% body_add_par("", style = "Normal")
    doc <- doc %>% body_add_fpar(fpar(ftext(sprintf(
      "Additional verified references not cited in the current draft (n = %d): available in references_verified.csv.",
      nrow(unused)), fp_cap)))
  }
}
doc
}

###############################################################################
# 第一遍：收集引用顺序；第二遍：正式生成
###############################################################################
invisible(build_paper())
keys <- CITE_ENV$order
CITE_MAP <- as.list(setNames(seq_along(keys), keys))
CITE_ENV$order <- character(0)   # 第二遍重新收集，保证顺序与编号一致
final <- build_paper()

stopifnot(identical(CITE_ENV$order, keys))

out_file <- "06_paper/Manuscript_final_v4.docx"
print(final, target = out_file)

cat("Written:", out_file, "\n")
cat("Citations used in text:", length(keys), "of", nrow(refs_all), "verified references\n")
cat("Citation order:\n"); print(data.frame(n = seq_along(keys), key = keys))
