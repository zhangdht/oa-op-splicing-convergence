###############################################################################
# 10_figures_4_6.R  —— 发表级图 4-6
#   Fig4  RNA 剪接程序作为间充质共享轴（富集 / 热图 / 样本水平 / 特异性）
#   Fig5  单细胞软骨细胞定位（UMAP / 组成 / 程序分）
#   Fig6  SON 定位 + ML 迁移 + 遗传支持（或局限）
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(Seurat)
})
ROOT <- "D:/bioinfo05"; setwd(ROOT)
set.seed(1)
BASE <- 7
th <- theme_bw(base_size = BASE, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = .2, colour = "grey92"),
        panel.border = element_rect(linewidth = .45, colour = "grey25"),
        axis.text  = element_text(colour = "grey15", size = BASE - .5),
        axis.title = element_text(colour = "black", size = BASE + .5),
        plot.title = element_text(face = "bold", size = BASE + 1.5, hjust = 0),
        plot.subtitle = element_text(size = BASE - .5, colour = "grey35"),
        legend.key.size = unit(.32, "cm"),
        legend.text = element_text(size = BASE - 1),
        legend.title = element_text(size = BASE - .5, face = "bold"),
        strip.background = element_rect(fill = "grey94", colour = "grey25", linewidth = .4),
        strip.text = element_text(size = BASE - .5, face = "bold"))
theme_set(th)
COL_OA <- "#C0392B"; COL_OP <- "#2471A3"; COL_MES <- "#B7791F"; COL_HEM <- "#1F8A70"
tag_th <- theme(plot.tag = element_text(face = "bold", size = BASE + 4))

save_fig <- function(p, name, w, h) {
  ggsave(sprintf("04_figures/%s.png", name), p, width = w, height = h, dpi = 300, units = "in", bg = "white")
  ggsave(sprintf("04_figures/%s.tiff", name), p, width = w, height = h, dpi = 300, units = "in", bg = "white", compression = "lzw")
  cat("saved", name, "\n")
}

## ============================================================ 数据 =========
sp_enrich <- fread("03_results/splicing/splicing_enrichment_by_set.csv")
sp_enrich <- sp_enrich[c(1, 2, 4, 3, 5)]
sp_enrich[, set_short := c("A: mesenchymal\nshared", "B: antagonistic\n(OA vs mono)",
                             "OA only", "OP BM-MSC\nonly", "OP monocyte\nonly")]
sp_enrich[, set_short := factor(set_short, levels = set_short)]

pmOA <- fread("03_results/gsva/pathway_meta_OA.csv")
pmMSC <- fread("03_results/gsva/pathway_meta_OP.csv")
setnames(pmMSC, c("pathway","k","z_meta","delta_mean","concord","P_meta","FDR_meta"))
A <- fread("03_results/gsva/compartment_pathways_A_mesenchymal.csv")
setnames(A, "FDR_MSC", "f_MSC")
sp_pathways <- fread("03_results/splicing/A_axis_splicing_pathways.csv")$pathway
SPLICE_RX <- "SPLIC|SPLICEOSOME|MRNA_PROCESS|MRNA_METABOL|MRNA_CATABOL|MRNA_TRANSPORT|NUCLEAR_SPECK|RNA_SPLIC|SNRNP|EXON_JUNCTION"
prog <- fread("03_results/splicing/programme_scores_percohort.csv")
prog[, cohort := factor(cohort)]
t1 <- fread("05_tables/Table1_cohorts.csv")
setnames(t1, c("cohort","disease","tissue","platform","n_control","n_case","genes"))
prog <- merge(prog, t1[, .(cohort, n_control, n_case)], by = "cohort", all.x = TRUE)
prog[, n := n_control + n_case]

son_de <- fread("03_results/splice/SON_per_cohort_DE.csv")
son_corr <- fread("03_results/splicing/SON_correlation_deconfounded.csv")
son_de <- merge(son_de, unique(son_corr[, .(cohort, tissue)]), by = "cohort", all.x = TRUE)
son_de[, t := sign(logFC) * qt(P/2, df = pmax(1L, n - 2), lower.tail = FALSE)]
son_de[, se := abs(logFC) / pmax(1e-6, abs(t))]
son_de[, xlo := logFC - 1.96 * se]
son_de[, xhi := logFC + 1.96 * se]
ml <- fread("03_results/ml/transfer_auc_full.csv")
ml_null <- fread("03_results/ml/random_geneset_null.csv")

## ============================================================ FIG 4 ========
f4a <- ggplot(sp_enrich, aes(set_short, fold, fill = n_splice > 0)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = .35, colour = "grey40") +
  geom_col(width = .6, colour = "grey20", linewidth = .25) +
  geom_text(aes(label = sprintf("%d/%d\nP=%.2g", n_splice, n_set, P_hyper)),
            vjust = -.25, size = 1.9, lineheight = .9) +
  scale_fill_manual(values = c("TRUE" = COL_MES, "FALSE" = "grey75"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .32))) +
  labs(x = NULL, y = "Splicing-pathway enrichment\n(fold over background)",
       title = "Splicing programmes are over-represented in the mesenchymal axis") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = BASE - 2, lineheight = .85))

## 4b 热图：A 轴剪接通路 z-score（OA vs MSC）
A[, splicing := grepl(SPLICE_RX, pathway)]
A[, short_name := gsub("^GOBP_", "", pathway)]
A[, short_name := gsub("_", " ", short_name)]
A[, short_name := substring(short_name, 1, 40)]
A_top <- A[splicing == TRUE][order(score)][1:min(13, .N)]
hm <- melt(A_top[, .(short_name, z_OA, z_MSC)], id.vars = "short_name",
           variable.name = "tissue", value.name = "z")
hm[, tissue := factor(tissue, levels = c("z_OA", "z_MSC"),
                      labels = c("OA cartilage", "OP BM-MSC"))]
f4b <- ggplot(hm, aes(tissue, short_name, fill = z)) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = sprintf("%.1f", z)), size = 1.9,
            colour = ifelse(abs(hm$z) > 3.5, "white", "grey15")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = expression(italic(Z))) +
  labs(x = NULL, y = NULL, title = "Top splicing pathways in the mesenchymal axis",
       subtitle = "all concordantly down-regulated") +
  theme(axis.text.y = element_text(size = BASE - 1.5), panel.grid = element_blank())

## 4c 样本水平 programme 评分（每个队列的 Hedges' g）
prog_spl <- prog[programme == "splicing"]
prog_spl[, cohort := factor(cohort, levels = prog_spl[, .(g = -hedges_g), by = cohort][order(g), cohort])]
prog_spl[, cls := ifelse(disease == "OA", COL_OA, COL_OP)]
prog_spl[, compartment := ifelse(tissue == "monocyte", "Haematopoietic", "Mesenchymal")]

f4c <- ggplot(prog_spl, aes(x = hedges_g, y = cohort, colour = cls, shape = compartment)) +
  geom_vline(xintercept = 0, linewidth = .35, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = hedges_g - 1.96 * hedges_g / abs(t) * sqrt(2 / n),
                     xmax = hedges_g + 1.96 * hedges_g / abs(t) * sqrt(2 / n)),
                 height = 0, linewidth = .55) +
  geom_point(size = 3) +
  scale_colour_identity(guide = "none") +
  scale_shape_manual(values = c("Mesenchymal" = 16, "Haematopoietic" = 17), name = NULL) +
  labs(x = "Hedges' g (Case vs Control)", y = NULL,
       title = "Splicing programme score in OA and OP cohorts",
       subtitle = "negative = down in disease") +
  theme(legend.position = "top")

## 4d 特异性对照：4 类程序在 A 轴的富集
ctrl_enr <- fread("03_results/splicing/splicing_enrichment_by_set.csv")
prog_enr <- prog[,.(P_meta=min(P)), by=.(programme)]
prog_enr[, programme := tools::toTitleCase(programme)]
prog_enr[, sig := ifelse(P_meta < .05, "P<0.05", "ns")]
f4d <- ggplot(prog_enr, aes(programme, -log10(P_meta), fill = sig)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = .3, colour = "grey50") +
  geom_col(width = .5, colour = "grey20", linewidth = .25) +
  scale_fill_manual(values = c("P<0.05" = COL_MES, "ns" = "grey78"), guide = "none") +
  labs(x = NULL, y = expression(-log[10]~italic(P)[meta]),
       title = "Specificity: only the splicing programme is consistently down") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = BASE - 1.5))

FIG4 <- (f4a | f4b) / (f4c | f4d) +
  plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "a") & tag_th
save_fig(FIG4, "Fig4_splicing_axis", 7.2, 6.8)

## ============================================================ FIG 5 ========
seu <- readRDS("03_results/intermediate/seurat_GSE324993.rds")
meta <- as.data.table(seu@meta.data, keep.rownames = "cell")
meta[, UMAP_1 := seu@reductions$umap@cell.embeddings[, 1]]
meta[, UMAP_2 := seu@reductions$umap@cell.embeddings[, 2]]
meta[, sample := factor(gsub("Sample_", "", sample))]
meta[, group := factor(group, levels = c("Healthy", "OA"))]

ss <- min(20000, nrow(meta))
idx <- sample(nrow(meta), ss)
f5a <- ggplot(meta[idx], aes(UMAP_1, UMAP_2, colour = celltype_short)) +
  geom_point(size = .15, alpha = .35) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  labs(colour = "Cell type", title = "Chondrocyte subtypes") +
  theme(legend.position = "right")

f5b <- ggplot(meta[idx], aes(UMAP_1, UMAP_2, colour = group)) +
  geom_point(size = .15, alpha = .35) +
  scale_colour_manual(values = c("Healthy" = "#7FB3D5", "OA" = COL_OA), name = NULL) +
  labs(title = "Condition")

## 组成条形图
comp <- fread("03_results/singlecell/celltype_composition.csv")
comp_long <- melt(comp, id.vars = "celltype_short")
comp_long[, sample := gsub("HC", "HC", gsub("OA", "OA", variable))]
comp_long[, condition := gsub("[0-9]", "", variable)]
comp_long[, pct := value * 100]
comp_summary <- comp_long[, .(mean_pct = mean(pct), sd_pct = sd(pct)),
                          by = .(celltype_short, condition)]
f5c <- ggplot(comp_summary, aes(celltype_short, mean_pct, fill = condition)) +
  geom_col(position = position_dodge(width = .7), width = .6, colour = "grey20", linewidth = .25) +
  geom_errorbar(aes(ymin = pmax(0, mean_pct - sd_pct), ymax = mean_pct + sd_pct),
                position = position_dodge(width = .7), width = .2, linewidth = .3) +
  scale_fill_manual(values = c("HC" = "#7FB3D5", "OA" = COL_OA), name = NULL) +
  labs(x = NULL, y = "Mean composition (%)", title = "Cell-type composition (3 HC vs 3 OA)",
       subtitle = "no significant shifts by chi-square test") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE - 1.5))

## 程序分（按细胞类型）
prog_sc <- fread("03_results/singlecell/programme_by_celltype_test.csv")
prog_sc[, sig := ifelse(FDR < .05, "FDR<0.05", "ns")]
prog_sc[, celltype_short := factor(celltype_short,
                                    levels = prog_sc[programme == "Splicing"][order(hedges_g), celltype_short])]
f5d <- ggplot(prog_sc[programme == "Splicing"], aes(hedges_g, celltype_short, fill = sig)) +
  geom_vline(xintercept = 0, linewidth = .3, linetype = 2, colour = "grey50") +
  geom_col(width = .6, colour = "grey20", linewidth = .25) +
  scale_fill_manual(values = c("FDR<0.05" = COL_MES, "ns" = "grey78"), guide = "none") +
  labs(x = "Hedges' g (OA vs HC)", y = NULL,
       title = "Splicing score by chondrocyte subtype",
       subtitle = "none survive FDR correction") +
  theme(axis.text.y = element_text(size = BASE - 1.5))

FIG5 <- (f5a | f5b) / (f5c | f5d) +
  plot_layout(heights = c(1.2, 1)) + plot_annotation(tag_levels = "a") & tag_th
save_fig(FIG5, "Fig5_single_cell_landscape", 7.2, 6.6)

## ============================================================ FIG 6 ========
## 6a SON 差异表达森林图
son_de[, compartment := ifelse(tissue == "monocyte", "Haematopoietic", "Mesenchymal")]
son_de[, cohort := factor(cohort, levels = son_de[order(logFC), cohort])]
son_de[, sig := ifelse(FDR < .05, "FDR<0.05", ifelse(P < .05, "P<0.05", "ns"))]
f6a <- ggplot(son_de, aes(logFC, cohort, colour = sig)) +
  geom_vline(xintercept = 0, linewidth = .3, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = xlo, xmax = xhi),
                 height = 0, linewidth = .55) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = c("FDR<0.05" = COL_OA, "P<0.05" = "#D35400", "ns" = "grey70"),
                      name = NULL) +
  labs(x = expression(SON~log[2]~FC), y = NULL,
       title = "SON is down-regulated in OA cartilage and OP BM-MSC") +
  theme(legend.position = "top")

## 6b SON-剪接相关 vs 随机零分布
son_corr[, compartment := ifelse(tissue == "monocyte", "Haematopoietic", "Mesenchymal")]
son_corr[, cohort := factor(cohort, levels = son_corr[order(rho_with_SON_in_set), cohort])]
son_corr[, sig := ifelse(emp_P < .05, "P<0.05", "ns")]
f6b <- ggplot(son_corr, aes(rho_with_SON_in_set, cohort, fill = sig)) +
  geom_vline(xintercept = 0, linewidth = .3, linetype = 2, colour = "grey50") +
  geom_col(width = .55, colour = "grey20", linewidth = .25) +
  scale_fill_manual(values = c("P<0.05" = COL_MES, "ns" = "grey75"), guide = "none") +
  labs(x = expression(Spearman~rho~(SON~vs~splicing~programme)), y = NULL,
       title = "SON tracks the splicing programme") +
  theme(axis.text.y = element_text(size = BASE - 1.5))

## 6c ML 迁移 AUC 点图
ml_sum <- ml[, .(AUC = sum(AUC * n_test, na.rm = TRUE) / sum(n_test[!is.na(AUC)]),
                 n = sum(n_test[!is.na(AUC)])), by = .(feature_set, scenario, method)]
ml_sum[, scenario2 := factor(scenario,
  levels = c("OA internal (LOCO)",
             "OA cartilage -> OP BM-MSC (same compartment)",
             "OA cartilage -> OP monocyte (cross compartment)",
             "OP BM-MSC -> OA cartilage (reverse)",
             "OP monocyte internal (positive ctrl)"))]
ml_sum[, short := c("OA internal", "OA -> BM-MSC", "OA -> monocyte",
                    "BM-MSC -> OA", "OP mono internal")[as.integer(scenario2)]]
ml_sum <- ml_sum[!is.na(AUC)]
ml_sum[, scenario3 := factor(short, levels = c("OA internal", "OA -> BM-MSC", "BM-MSC -> OA",
                                               "OA -> monocyte", "OP mono internal"))]
f6c <- ggplot(ml_sum[feature_set == "Splicing"],
              aes(scenario3, AUC, colour = method, shape = method)) +
  geom_hline(yintercept = .5, linetype = 2, linewidth = .3, colour = "grey50") +
  geom_point(position = position_dodge(width = .5), size = 2.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  scale_colour_manual(values = c("LASSO" = COL_MES, "RF" = "#5D6D7E", "SVM" = "#1F8A70"),
                      name = NULL) +
  scale_shape_manual(values = c("LASSO" = 16, "RF" = 17, "SVM" = 15), name = NULL) +
  labs(x = NULL, y = "AUC", title = "Splicing classifier transfers weakly to OP") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = BASE - 2, lineheight = .85),
        legend.position = "top")

## 6d 遗传支持 / 研究局限（GWAS sumstats 未完整下载）
lim_txt <- paste0(
  "GWAS summary-statistics limitations\n",
  "• Full-scale MR was not completed: full GWAS summary statistics\n",
  "  for eBMD, fracture and OA exceeded practical download limits.\n",
  "• Genetic evidence is therefore limited to published catalogues\n",
  "  and previously reported OP-OA MR relationships (see text).\n",
  "• Findings are supported at the transcriptomic but not\n",
  "  at the germline-causal level in this study."
)
f6d <- ggplot() + annotate("text", x = 1, y = 1, size = 2.5, hjust = .5, vjust = .5,
    label = lim_txt, lineheight = 1.15, colour = "grey20") +
  theme_void(base_size = BASE) +
  xlim(0, 2) + ylim(0, 2) +
  labs(title = "Genetic evidence caveat")

FIG6 <- (f6a | f6b) / (f6c | f6d) +
  plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "a") & tag_th
save_fig(FIG6, "Fig6_SON_ML_genetics", 7.2, 6.4)

cat("FIG 4-6 done\n")
