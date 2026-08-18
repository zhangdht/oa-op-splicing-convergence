###############################################################################
# 11_tables_word.R
# 生成投稿级 Word 表格（officer/flextable），同时输出 CSV。
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(officer); library(flextable); library(magrittr)
})
ROOT <- "D:/bioinfo05"; setwd(ROOT)
dir.create("06_paper", showWarnings = FALSE)

ft_theme <- function(ft, fs = 7) {
  ft %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = fs, part = "all") %>%
    bold(part = "header") %>%
    border_outer(border = fp_border(color = "black", width = .5)) %>%
    border_inner_h(border = fp_border(color = "grey70", width = .25)) %>%
    border_inner_v(border = fp_border(color = "grey70", width = .25)) %>%
    padding(padding = 2, part = "all") %>%
    set_table_properties(layout = "autofit")
}

# -----------------------------------------------------------------------------
# Table 1: cohorts
# -----------------------------------------------------------------------------
t1 <- fread("05_tables/Table1_cohorts.csv")
setnames(t1, c("Cohort","Disease","Tissue","Platform","Control","Case","Genes"))
t1_out <- data.table(
  Cohort = t1$Cohort,
  Disease = t1$Disease,
  Tissue = t1$Tissue,
  Platform = t1$Platform,
  `n Control` = t1$Control,
  `n Case` = t1$Case,
  `n Total` = t1$Control + t1$Case,
  Genes = t1$Genes
)
fwrite(t1_out, "05_tables/Table1_cohorts_formatted.csv")
ft1 <- ft_theme(flextable(t1_out))

# -----------------------------------------------------------------------------
# Table 2: convergence controls + family-level
# -----------------------------------------------------------------------------
c2 <- fread("03_results/gsva/convergence_controls.csv")
comp_labs <- c("OA vs OP (pooled)",
               "OA RNA-seq vs OA array (positive ctrl)",
               "GSE114007 vs GSE57218",
               "GSE57218 vs GSE169077",
               "OP monocyte vs OP BM-MSC",
               "OA cartilage vs OP BM-MSC",
               "OA cartilage vs OP monocyte")
t2 <- data.table(
  Comparison = comp_labs,
  `n pathways` = c2$n_pathway,
  `n sig A` = c2$nsig_A,
  `n sig B` = c2$nsig_B,
  Observed = c2$observed,
  `Perm mean` = c2$perm_mean,
  Enrichment = round(c2$enrichment, 2),
  `Spearman rho` = round(c2$spearman_rho, 3),
  `Perm P` = c2$perm_P
)
fwrite(t2, "05_tables/Table2_convergence_controls.csv")
ft2 <- ft_theme(flextable(t2))

fam <- fread("03_results/gsva/family_level_convergence.csv")
comp_labs2 <- c("OA array vs OA RNA-seq (positive ctrl)",
                "OA cartilage vs OP BM-MSC (mesenchymal)",
                "OA cartilage vs OP monocyte (cross-compartment)",
                "OA vs OP (pooled)")
t2b <- data.table(
  Comparison = comp_labs2,
  `n families` = fam$n_family,
  `n sig A` = fam$sigA,
  `n sig B` = fam$sigB,
  Observed = fam$observed,
  Expected = round(fam$expected, 2),
  Enrichment = round(fam$fold, 2),
  `Perm P` = fam$P_perm,
  `Splice families` = fam$n_splice_conv
)
fwrite(t2b, "05_tables/Table2_family_level_convergence.csv")
ft2b <- ft_theme(flextable(t2b))

# -----------------------------------------------------------------------------
# Table 3: splicing axis
# -----------------------------------------------------------------------------
t3a <- fread("05_tables/Table3_splicing_axis.csv")
t3a_out <- data.table(
  Pathway = gsub("_", " ", t3a$Pathway),
  Source = t3a$Source,
  `Z OA` = round(t3a$Z_OA_cartilage, 2),
  `FDR OA` = formatC(t3a$FDR_OA, format = "e", digits = 1),
  `Z BM-MSC` = round(t3a$Z_OP_BMMSC, 2),
  `FDR BM-MSC` = formatC(t3a$FDR_OP, format = "e", digits = 1),
  Direction = t3a$Direction
)
fwrite(t3a_out, "05_tables/Table3_top_splicing_pathways.csv")
ft3a <- ft_theme(flextable(t3a_out))

sp_meta <- fread("03_results/splicing/programme_meta_by_disease.csv")
sp_meta[, programme := tools::toTitleCase(programme)]
t3b <- data.table(
  Programme = sp_meta$programme,
  Disease = sp_meta$disease,
  `n cohorts` = sp_meta$k,
  `Meta Z` = round(sp_meta$z_meta, 2),
  `Meta P` = formatC(sp_meta$P_meta, format = "e", digits = 1),
  `Mean Hedges g` = round(sp_meta$g_mean, 2)
)
fwrite(t3b, "05_tables/Table3_programme_meta.csv")
ft3b <- ft_theme(flextable(t3b))

# -----------------------------------------------------------------------------
# Word 文档
# -----------------------------------------------------------------------------
my_doc <- read_docx() %>%
  body_add_par("Table 1. Study cohorts", style = "heading 1") %>%
  body_add_flextable(ft1) %>%
  body_add_par("Table 2. Pathway-level convergence: gene-set enrichment controls (top) and family-level convergence (bottom)", style = "heading 1") %>%
  body_add_flextable(ft2) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft2b) %>%
  body_add_par("Table 3. RNA splicing programme: representative pathways (top) and meta-analytic programme scores (bottom)", style = "heading 1") %>%
  body_add_flextable(ft3a) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft3b)
print(my_doc, target = "06_paper/Tables.docx")

cat("Tables.docx generated\n")
