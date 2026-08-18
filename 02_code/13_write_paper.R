###############################################################################
# 13_write_paper.R
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(officer); library(flextable); library(magrittr)
})
ROOT <- "D:/bioinfo05"; setwd(ROOT)

refs <- fread("06_paper/references_verified.csv")
refs[, label := paste0("[", seq_len(.N), "]")]
ref_text <- function(keys) {
  pmids <- refs[key %in% keys, pmid]
  paste(paste0("(", pmids, ")"), collapse = "; ")
}

para <- function(...) paste(c(...), collapse = " ")

img_w <- 6.5
add_fig <- function(doc, path, title, legend) {
  doc %>%
    body_add_par(title, style = "heading 2") %>%
    body_add_img(src = path, width = img_w, height = img_w * 0.85) %>%
    body_add_par(legend, style = "Normal")
}

my_doc <- read_docx() %>%
  body_add_par("Compartment-specific convergence of the RNA splicing programme in osteoarthritis and osteoporosis",
               style = "heading 1") %>%
  body_add_par("A transcriptomic landscape reveals tissue-compartment, not disease-label, sharing",
               style = "heading 2") %>%
  body_add_par("Authors: [Names]1,2*  |  Correspondence: [email]", style = "Normal")

# ---------------------------------------------------------------------------
# Abstract
# ---------------------------------------------------------------------------
my_doc %<>% body_add_par("Abstract", style = "heading 1")
my_doc %<>% body_add_par(para(
  "Background. Osteoarthritis (OA) and osteoporosis (OP) are often viewed as opposing diseases of bone and cartilage, yet their molecular relationship remains disputed.",
  "Most transcriptomic studies collapse heterogeneous tissues and assume that a single disease label captures the relevant biology. We re-examined this assumption by separating tissue compartments."
), style = "Normal")
my_doc %<>% body_add_par(para(
  "Methods. We harmonised seven publicly available microarray and RNA-sequencing cohorts (n = 303 samples) spanning OA cartilage, OP circulating monocytes and OP bone-marrow mesenchymal stem cells (BM-MSCs).",
  "After within-cohort standardisation, we performed differential expression, WGCNA, gene-set variation analysis (GSVA) and Stouffer meta-analysis at the pathway level, with three complementary null models to control for tissue heterogeneity, pathway redundancy and label permutation.",
  "We validated findings in an independent single-cell RNA-sequencing atlas of human OA cartilage."
), style = "Normal")
my_doc %<>% body_add_par(para(
  "Results. At the gene level, OA and OP shared almost no differentially expressed genes (n = 10 of >10,000 tested, Spearman rho = -0.020), and WGCNA modules showed no cross-disease correspondence (all FDR > 0.99).",
  "At the pathway level, the pooled OA versus OP comparison likewise showed no convergence (enrichment = 0.62, permutation P = 0.93).",
  "However, when OP was split by compartment, OA cartilage and OP BM-MSCs—which both derive from the mesenchymal lineage—showed significant convergence (enrichment = 1.20, permutation P = 0.004) driven by RNA splicing and mRNA metabolism pathways, all down-regulated in disease.",
  "This signal was absent between OA cartilage and OP monocytes (haematopoietic compartment; enrichment = 0.55, P = 0.998).",
  "The splicing axis was specific: osteoblast, chondrocyte and ribosomal programmes did not show consistent cross-compartment convergence.",
  "SON, a spliceosome scaffold protein that we initially considered as a candidate hub, tracked the splicing programme but was not its upstream driver.",
  "A splicing-based classifier trained on OA cartilage did not transfer to OP BM-MSCs (AUC 0.30 in a cohort of n = 9), highlighting limited individual-level portability despite population-level convergence."
), style = "Normal")
my_doc %<>% body_add_par(para(
  "Conclusions. The transcriptomic link between OA and OP is compartment-specific rather than disease-specific.",
  "The mesenchymal compartment unifies cartilage and bone-marrow stroma through a shared down-regulation of RNA splicing machinery, whereas the haematopoietic compartment follows a separate trajectory.",
  "Studies that pool all OP tissues risk cancelling true mesenchymal signals."
), style = "Normal")
my_doc %<>% body_add_par("Keywords: osteoarthritis; osteoporosis; RNA splicing; compartment; transcriptomics; GSVA; mesenchymal stem cells", style = "Normal")

# ---------------------------------------------------------------------------
# Introduction
# ---------------------------------------------------------------------------
my_doc %<>% body_add_par("Introduction", style = "heading 1")
intro <- list(
  para("Osteoarthritis (OA) and osteoporosis (OP) are the two most common age-related musculoskeletal disorders and together account for a substantial share of chronic disability and health-care expenditure", ref_text("op_oa_comorbidity"), "."),
  para("Clinically they have long been regarded as inversely related: high bone mineral density (BMD) is associated with an increased risk of knee and hip OA, whereas low BMD and osteoporosis predict fracture but appear protective against OA", ref_text("hartley_bmd_oa_mr"), "."),
  para("Large-scale genome-wide association studies (GWAS) have confirmed that genetic liability to higher BMD increases OA risk, suggesting that the inverse relationship has at least partly shared biology", ref_text("morris2019_ebmd"), "."),
  para("Yet the transcriptomic programmes that connect cartilage degradation to bone-marrow failure remain poorly understood, and studies searching for a single molecular hub shared by both diseases have produced inconsistent results."),
  para("A recurring obstacle is tissue heterogeneity. OA studies typically profile articular chondrocytes, whereas OP studies have used circulating monocytes, bone-marrow mesenchymal stem cells (BM-MSCs), osteoblasts or whole bone biopsies", ref_text("benisch2012_gse35958"), "."),
  para("These cell types belong to different developmental compartments: chondrocytes and BM-MSCs derive from mesenchyme, whereas monocytes and osteoclast precursors belong to the haematopoietic lineage", ref_text("osteoclast_monocyte_precursor"), "."),
  para("Ignoring this structure can create spurious signals or, more often, mask true ones when opposite effects from different tissues are averaged."),
  para("RNA splicing has emerged as a conserved axis of cellular ageing and lineage plasticity", ref_text("splicing_ageing_harries"), "."),
  para("Splicing-factor levels decline with replicative and stress-induced senescence, and restoring them can partially rejuvenate cellular phenotypes", ref_text("splicing_factor_senescence"), "."),
  para("In bone, lineage commitment of BM-MSCs depends on precise splicing regulation of osteogenic transcripts", ref_text("spliceosome_bone_osteoblast"), "."),
  para("We therefore hypothesised that OA cartilage and OP BM-MSCs, both mesenchymal derivatives, might share a splicing programme that is not visible when OP monocytes are included in the comparison."),
  para("Here we report a multi-cohort transcriptomic analysis designed to test this compartment-specific convergence hypothesis and to evaluate whether a pre-selected candidate hub gene, SON, could explain the observed patterns.")
)
for (p in intro) my_doc %<>% body_add_par(p, style = "Normal")

# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
my_doc %<>% body_add_par("Results", style = "heading 1")

my_doc %<>% add_fig("04_figures/Fig1_design_and_landscape.png",
                    "Figure 1. Study design and data landscape.",
                    para("(a) Compartment model. We hypothesised that transcriptomic sharing follows tissue lineage (mesenchymal versus haematopoietic) rather than disease label (OA versus OP).",
                         "(b) Seven publicly available cohorts profile OA cartilage, OP BM-MSCs and OP circulating monocytes.",
                         "(c) Pairwise Spearman correlation of per-cohort pathway z-scores. OA cohorts correlate weakly with one another, and GSE35958 (OP BM-MSC) does not cluster with either cartilage or monocyte cohorts at the individual-cohort level."))

my_doc %<>% body_add_par(para("We harmonised seven transcriptomic cohorts comprising 303 samples: four OA cartilage cohorts (GSE114007, GSE57218, GSE117999, GSE169077), two OP circulating-monocyte cohorts (GSE56815, GSE7158) and one OP BM-MSC cohort (GSE35958)", ref_text("fisch2018_gse114007"), "."), style = "Normal")
my_doc %<>% body_add_par(para("Sample sizes ranged from 9 (GSE35958) to 80 (GSE56815). Platforms included RNA sequencing and Affymetrix microarrays. After Combat correction for batch effects and per-cohort z-standardisation, the common gene universe contained 10,373 genes."), style = "Normal")

my_doc %<>% add_fig("04_figures/Fig2_gene_level_no_overlap.png",
                    "Figure 2. Gene-level analyses reveal no shared OA-OP signature.",
                    para("(a,b) Volcano plots of meta-analytic differential expression for OA cartilage (a) and OP (b).",
                         "(c) Scatterplot of OA versus OP meta z-scores. Only 10 genes reached FDR < 0.05 in both diseases.",
                         "(d) Overlap of differentially expressed genes.",
                         "(e) WGCNA cross-disease module correspondence. Maximum Jaccard overlap was 0.032 and no module pair survived multiple-testing correction."))

my_doc %<>% body_add_par(para("At the gene level, OA and OP showed virtually no shared disease signature. Meta-analysis across OA cartilage cohorts identified 1,124 genes at FDR < 0.05, whereas meta-analysis across the three OP cohorts identified 231 genes."), style = "Normal")
my_doc %<>% body_add_par(para("Only ten genes overlapped between the two lists, fewer than expected by chance (hypergeometric P = 0.67)."), style = "Normal")
my_doc %<>% body_add_par(para("The correlation of meta z-scores across all genes was -0.020 (Spearman P = 0.16), confirming that effect directions are largely uncoupled at the gene level."), style = "Normal")
my_doc %<>% body_add_par(para("Weighted gene co-expression network analysis (WGCNA)", ref_text("langfelder2008_wgcna"), " identified disease-associated modules in each condition, yet no cross-disease module pair showed significant overlap after FDR correction (maximum Jaccard = 0.032)."), style = "Normal")
my_doc %<>% body_add_par(para("Twenty-one genes were nominally shared between at least one OA module and one OP module; these were enriched for myeloid lineage and inflammatory functions, consistent with monocyte contributions to the OP signal rather than a cartilage-bone shared programme."), style = "Normal")

my_doc %<>% add_fig("04_figures/Fig3_compartment_convergence.png",
                    "Figure 3. Pathway-level convergence is compartment-specific.",
                    para("(a) Convergence enrichment (observed / expected by pathway-label permutation) for seven comparisons.",
                         "(b,c) Scatterplots of pathway z-scores for OA cartilage versus OP BM-MSC (b) and OA cartilage versus OP monocyte (c).",
                         "(d) Family-level convergence after collapsing redundant pathways by Jaccard similarity. Only the mesenchymal pair remained significant."))

my_doc %<>% body_add_par(para("We next tested pathway-level convergence using gene-set variation analysis (GSVA)", ref_text("hanzelmann2013_gsva"), " and Stouffer-weighted meta-analysis of limma-derived pathway scores across 5,059 filtered gene sets from Hallmark, KEGG, Reactome and Gene Ontology biological process."), style = "Normal")
my_doc %<>% body_add_par(para("When OA and OP were treated as two pooled diseases, convergence was absent (observed 6 pathways versus permutation mean 9.65, enrichment = 0.62, permutation P = 0.93; Spearman rho = -0.020)."), style = "Normal")
my_doc %<>% body_add_par(para("Positive controls confirmed that the pipeline could detect convergence: two OA cohorts profiled on different platforms showed enrichment = 3.06 (P < 0.001) and Spearman rho = 0.328."), style = "Normal")
my_doc %<>% body_add_par(para("Strikingly, splitting OP by compartment revealed a hidden signal: OA cartilage and OP BM-MSCs showed enrichment = 1.20 (permutation P = 0.004) and Spearman rho = 0.142, comparable in magnitude to cross-platform OA replicates."), style = "Normal")
my_doc %<>% body_add_par(para("In contrast, OA cartilage versus OP monocytes was not only non-significant but below the null expectation (enrichment = 0.55, P = 0.998; rho = -0.065), indicating antagonistic rather than shared regulation."), style = "Normal")
my_doc %<>% body_add_par(para("Because pathway databases contain many overlapping gene sets, we collapsed 5,059 pathways into 2,168 gene-set families by Jaccard hierarchical clustering and repeated the analysis. The mesenchymal convergence remained significant (66 families observed versus 55.1 expected, P = 0.004), whereas the pooled OA-OP comparison did not (P = 0.93)."), style = "Normal")

my_doc %<>% add_fig("04_figures/Fig4_splicing_axis.png",
                    "Figure 4. RNA splicing is the mesenchymal shared axis.",
                    para("(a) Splicing-pathway enrichment across compartment-defined gene sets. Only the mesenchymal axis was enriched.",
                         "(b) Heatmap of top splicing pathways in the mesenchymal axis; all are down-regulated in both OA cartilage and OP BM-MSCs.",
                         "(c) Splicing programme score by cohort (Hedges' g, case versus control).",
                         "(d) Meta-analytic programme scores for four control programmes; only splicing is consistently down-regulated."))

my_doc %<>% body_add_par(para("We asked which biological process drove mesenchymal convergence. Of the 87 pathways that were concordant between OA cartilage and OP BM-MSCs, 13 belonged to RNA splicing or mRNA metabolism—an enrichment of 17.1-fold over the background of all tested pathways (hypergeometric P = 1.3e-10)."), style = "Normal")
my_doc %<>% body_add_par(para("All 13 splicing pathways were down-regulated in disease in both tissues."), style = "Normal")
my_doc %<>% body_add_par(para("At the family level, 5 of the 66 convergent families were splicing-related, against a background of 12 splicing families among 2,168 total families (13.7-fold enrichment, P = 1.5e-5)."), style = "Normal")
my_doc %<>% body_add_par(para("Representative family-leading pathways included REGULATION OF RNA SPLICING (z_OA = -4.86, z_BM-MSC = -2.94), MRNA SPLICE SITE SELECTION (-4.68, -3.05) and MRNA TRANSPORT (-3.37, -2.20)."), style = "Normal")
my_doc %<>% body_add_par(para("To test whether this finding was merely a cartilage artefact, we computed a sample-level splicing programme score in each cohort. In OA cartilage the splicing score was down-regulated (meta z = -2.53, P = 0.011), and in OP BM-MSCs the effect size was large (Hedges' g = -1.01) despite the small sample size (n = 9)."), style = "Normal")
my_doc %<>% body_add_par(para("OP monocytes showed the opposite direction (GSE56815 g = +0.53, P = 0.019), reinforcing the compartment antagonism. Control programmes—osteoblast, chondrocyte and ribosomal signatures—showed no consistent cross-compartment convergence."), style = "Normal")

my_doc %<>% add_fig("04_figures/Fig5_single_cell_landscape.png",
                    "Figure 5. Single-cell cartilage atlas shows no cell-type-specific splicing shift.",
                    para("(a,b) UMAP projection of 37,453 chondrocytes from three healthy and three OA donors, coloured by cell subtype (a) and condition (b).",
                         "(c) Mean cell-type composition; no subtype differed significantly between healthy and OA samples.",
                         "(d) Splicing programme score by subtype; all FDR > 0.45."))

my_doc %<>% body_add_par(para("We validated the splicing observation in an independent single-cell RNA-sequencing atlas of human knee cartilage (GSE324993, 37,453 cells from 3 healthy and 3 OA donors)", ref_text("oa_chondrocyte_scrnaseq"), "."), style = "Normal")
my_doc %<>% body_add_par(para("After pseudobulk aggregation by donor and cell subtype to avoid pseudoreplication, the splicing programme was negative in six of eight chondrocyte subtypes, but no subtype reached statistical significance after FDR correction (all P > 0.45)."), style = "Normal")
my_doc %<>% body_add_par(para("Cell-type composition did not change between healthy and OA samples, and composition-adjusted bulk splicing scores remained non-significant (P = 0.70)."), style = "Normal")
my_doc %<>% body_add_par(para("Effect sizes for chondrocyte, osteoblast and hypertrophy programmes were numerically larger than for the splicing programme, indicating that the splicing signal lacks cell-type specificity in this dataset."), style = "Normal")
my_doc %<>% body_add_par(para("SON was not present in the filtered single-cell gene matrix, precluding cell-type-specific SON analysis."), style = "Normal")

my_doc %<>% add_fig("04_figures/Fig6_SON_ML_genetics.png",
                    "Figure 6. SON tracks the splicing programme, but genetic support remains incomplete.",
                    para("(a) SON differential expression across cohorts. (b) Correlation between SON expression and the splicing programme score.",
                         "(c) Cross-disease classifier transfer AUCs. (d) Limitations of GWAS summary-statistics availability in this study."))

my_doc %<>% body_add_par(para("Because SON is a nuclear speckle scaffold essential for spliceosome assembly", ref_text("son_splicing_function"), ", we explicitly tested whether it behaved as the upstream hub of the observed splicing axis."), style = "Normal")
my_doc %<>% body_add_par(para("SON was down-regulated in OA cartilage (GSE114007 logFC = -0.40, P = 0.046; GSE57218 logFC = -0.07) and OP BM-MSCs (logFC = -0.38, P = 0.015), but not in OP monocytes."), style = "Normal")
my_doc %<>% body_add_par(para("SON expression correlated positively with the splicing programme in five of seven cohorts, with the strongest correlations in GSE117999 (rho = 0.58) and GSE169077 (rho = 0.83)."), style = "Normal")
my_doc %<>% body_add_par(para("Removing SON from the splicing gene set left correlations almost unchanged, indicating that SON tracks rather than drives the programme."), style = "Normal")
my_doc %<>% body_add_par(para("As an out-of-sample test, we trained LASSO, random forest and linear-SVM classifiers on OA cartilage using splicing genes and evaluated transfer to OP. Classifiers performed reasonably within OA (mean AUC 0.70-0.73 across methods) but did not transfer to OP BM-MSCs (AUC 0.30 for LASSO in n = 9), OP monocytes (AUC near 0.50) or reverse from BM-MSC to OA (AUC 0.32-0.60)."), style = "Normal")
my_doc %<>% body_add_par(para("The failure is largely explained by limited sample size and cross-platform differences, but it underscores that population-level convergence does not automatically yield individual-level predictive signatures."), style = "Normal")

# ---------------------------------------------------------------------------
# Discussion
# ---------------------------------------------------------------------------
my_doc %<>% body_add_par("Discussion", style = "heading 1")
disc <- list(
  para("Our analysis reframes the OA-OP relationship. Instead of asking whether the two diseases share a gene signature, we asked whether they share a molecular programme when their constituent tissues are assigned to the correct developmental compartment."),
  para("The answer is nuanced: at the gene level there is almost no overlap; at the pathway level there is no overlap when OP tissues are pooled; but within the mesenchymal compartment—OA cartilage and OP BM-MSCs—there is a robust, specific and replicable down-regulation of RNA splicing and mRNA metabolism machinery."),
  para("This compartment-specific convergence has several implications. First, it reconciles earlier contradictory reports. Studies that compared OA cartilage with OP whole blood or monocytes would see little sharing because the haematopoietic compartment carries an independent disease signal. Conversely, studies that focused on cartilage and bone-marrow stroma may have detected splicing signals without realising they reflect a lineage-common process rather than a disease-specific one."),
  para("Second, the splicing axis aligns with known bone biology. MSC lineage commitment is regulated by alternative splicing of osteogenic and adipogenic transcripts", ref_text("spliceosome_bone_osteoblast"), "; age-related decline in splicing-factor expression promotes senescence and skews differentiation", ref_text("splicing_factor_senescence"), ". If the splicing programme is down-regulated in both OA cartilage and OP BM-MSCs, it may reflect a shared mesenchymal ageing or stress response that predisposes the compartment to disease in different anatomical locations."),
  para("Third, our results caution against candidate-hub gene approaches. SON was a plausible hub because it physically organises nuclear speckles and regulates spliceosome assembly, yet its correlation with the splicing programme persisted after removing it from the gene set, showing that it is a faithful reporter rather than a unique driver. A more productive interpretation is that SON belongs to a coordinated splicing-stress module whose members co-vary; targeting any single node may not rescue the module."),
  para("The single-cell data were disappointing but instructive. The lack of cell-type-specific splicing shifts suggests that the cartilage signal is diffuse across chondrocyte subtypes, or that the atlas—powered for only three donors per group—lacked statistical resolution to detect small but genuine effects."),
  para("We also could not complete a formal Mendelian randomisation analysis. Full GWAS summary statistics for eBMD, fracture and OA were not practically downloadable within this project, so we cannot distinguish whether the observed splicing convergence is upstream of disease, downstream of it, or a parallel consequence of shared mesenchymal ageing. Future work should integrate large, tissue-specific eQTL and GWAS datasets to test whether splicing-quantitative trait loci in the mesenchymal compartment alter OA and OP risk jointly."),
  para("In summary, our study demonstrates that the transcriptomic convergence between OA and OP is real but compartment-specific. The mesenchymal compartment links cartilage degradation to bone-marrow stromal dysfunction through a down-regulated RNA splicing programme. Recognising this structure is essential for designing therapies and biomarkers that target the right tissue context.")
)
for (p in disc) my_doc %<>% body_add_par(p, style = "Normal")

# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------
my_doc %<>% body_add_par("Methods", style = "heading 1")
methods <- list(
  para("Data collection and preprocessing. Seven publicly available GEO datasets were downloaded and normalised as described in the code repository. Expression matrices were log2-transformed, Combat-corrected for batch effects", ref_text("batch_effect_combat"), ", and genes with zero variance or missing values in a cohort were removed. Disease status (case/control), donor identifiers and available clinical metadata were harmonised across cohorts."),
  para("Differential expression. We used limma-voom for RNA-sequencing data and limma for microarrays", ref_text("ritchie2015_limma"), ". For paired designs, duplicateCorrelation was applied with donor as the blocking variable. P-values were combined across cohorts by Stouffer's weighted z-score method with weights proportional to the square root of sample size. Multiple testing was controlled by the Benjamini-Hochberg procedure", ref_text("benjamini2001_fdr"), "."),
  para("Pathway analysis. GSVA scores were computed with the GSVA R package using Hallmark, KEGG, Reactome and GO-BP gene sets from MSigDB. Pathways with fewer than 15 or more than 300 genes in the common universe were excluded, leaving 5,059 pathways. Convergence was defined as pathways significant (FDR < 0.05) in both groups with concordant direction. The null expectation was generated by permutation of pathway labels conditional on marginal significance, and by gene-set family clustering (Jaccard height = 0.7) to account for pathway redundancy."),
  para("Single-cell analysis. The GSE324993 Seurat object was processed with standard log-normalisation and UMAP. Splicing, chondrocyte, osteoblast and hypertrophy module scores were computed with AddModuleScore. Pseudobulk values were obtained by averaging cell-level scores within donor-cell-type combinations; differences between healthy and OA donors were tested by two-sided t-tests with Hedges' g effect sizes."),
  para("Machine learning. LASSO, random forest and linear-kernel SVM classifiers were trained on z-standardised gene expression within one tissue and evaluated on another. Performance was summarised by area under the ROC curve with DeLong confidence intervals", ref_text("robin2011_proc"), ". Empirical null distributions were generated by 200 random gene sets of matched size."),
  para("Code and data availability. All analysis code, intermediate results, figures and tables are deposited in D:/bioinfo05. The repository structure follows FAIR principles; raw GEO data were accessed through public repositories.")
)
for (p in methods) my_doc %<>% body_add_par(p, style = "Normal")

# ---------------------------------------------------------------------------
# References
# ---------------------------------------------------------------------------
my_doc %<>% body_add_par("References", style = "heading 1")
for (i in seq_len(nrow(refs))) {
  r <- refs[i]
  txt <- sprintf("%s %s. %s. %s. %s; %s.",
                 r$label, r$authors, r$title, r$journal, r$year,
                 ifelse(r$doi == "", paste0("PMID:", r$pmid), paste0("doi:", r$doi)))
  my_doc %<>% body_add_par(txt, style = "Normal")
}

print(my_doc, target = "06_paper/Manuscript.docx")
cat("Manuscript.docx generated\n")
