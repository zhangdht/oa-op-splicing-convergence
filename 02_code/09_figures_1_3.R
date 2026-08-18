###############################################################################
# 09_figures_1_3.R  —— 发表级图 1-3
#   Fig1  研究设计与数据全景（分区模型 / 队列构成 / 队列间差异特征相关性）
#   Fig2  基因层面：两病无共享特征（火山 / 一致性散点 / 重叠 / WGCNA 模块对应）
#   Fig3  通路层面：分区特异收敛（对照条形 / A轴散点 / B轴散点 / 家族级 / 模型示意）
# 输出 04_figures/Fig1.tiff|png ... 300 dpi
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(ggrepel); library(scales); library(grid)
})
ROOT <- "D:/bioinfo05"; setwd(ROOT)
dir.create("04_figures", showWarnings = FALSE)
set.seed(1)

## --------------------------------------------------------- 全局主题 --------
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

COL_OA   <- "#C0392B"; COL_OP <- "#2471A3"
COL_MES  <- "#B7791F"; COL_HEM <- "#1F8A70"
COL_NEU  <- "grey55"
tag_th <- theme(plot.tag = element_text(face = "bold", size = BASE + 4))

save_fig <- function(p, name, w, h) {
  ggsave(sprintf("04_figures/%s.png", name), p, width = w, height = h,
         dpi = 300, units = "in", bg = "white")
  ggsave(sprintf("04_figures/%s.tiff", name), p, width = w, height = h,
         dpi = 300, units = "in", bg = "white", compression = "lzw")
  cat("saved", name, "\n")
}

## ============================================================ 数据 =========
tab1  <- fread("05_tables/Table1_cohorts.csv")
deOA  <- fread("03_results/deg/META_OA.csv")
deOP  <- fread("03_results/deg/META_OP.csv")
wov   <- fread("03_results/wgcna/cross_disease_module_overlap.csv")
pde   <- fread("03_results/gsva/pathway_DE_percohort.csv")
ctrl  <- fread("03_results/gsva/convergence_controls.csv")
fam   <- fread("03_results/gsva/family_level_convergence.csv")
pmOA  <- fread("03_results/gsva/pathway_meta_OA.csv")
pmOP  <- fread("03_results/gsva/pathway_meta_OP.csv")

OA_G   <- c("GSE114007","GSE57218","GSE117999","GSE169077")
MSC_G  <- "GSE35958"; MONO_G <- c("GSE56815","GSE7158")

meta_set <- function(dt) {
  dt <- dt[is.finite(P) & P > 0 & P < 1]
  dt[, z := qnorm(P/2, lower.tail = FALSE) * sign(delta)]
  dt[, w := sqrt(n)]
  out <- dt[, .(z = sum(z*w)/sqrt(sum(w^2))), by = pathway]
  out[, FDR := p.adjust(2*pnorm(-abs(z)), "BH")][]
}
m_oa <- meta_set(pde[cohort %in% OA_G]); m_msc <- meta_set(pde[cohort %in% MSC_G])
m_mono <- meta_set(pde[cohort %in% MONO_G])

###############################################################################
## FIGURE 1
###############################################################################
## --- 1a 分区模型示意 -------------------------------------------------------
seg <- data.table(
  x = c(1, 1, 1), xend = c(2.6, 2.6, 2.6),
  y = c(3, 2, 1), yend = c(3, 2, 1)
)
boxes <- data.table(
  x = c(0.5, 0.5, 0.5, 3.4, 3.4),
  y = c(3, 2, 1, 2.5, 1),
  lab = c("OA articular\ncartilage", "OP bone-marrow\nMSC", "OP circulating\nmonocyte",
          "MESENCHYMAL\ncompartment", "HAEMATOPOIETIC\ncompartment"),
  type = c("t","t","t","c","c"),
  col = c(COL_OA, COL_OP, COL_OP, COL_MES, COL_HEM)
)
lk <- data.table(x = c(1.15, 1.15, 1.15), xend = c(2.75, 2.75, 2.75),
                 y = c(3, 2, 1), yend = c(2.75, 2.35, 1))
f1a <- ggplot() +
  geom_curve(data = lk[1:2], aes(x = x, y = y, xend = xend, yend = yend),
             curvature = .12, linewidth = .5, colour = COL_MES,
             arrow = arrow(length = unit(.07, "in"), type = "closed")) +
  geom_curve(data = lk[3], aes(x = x, y = y, xend = xend, yend = yend),
             curvature = -.12, linewidth = .5, colour = COL_HEM,
             arrow = arrow(length = unit(.07, "in"), type = "closed")) +
  geom_label(data = boxes[type == "t"], aes(x, y, label = lab, colour = col),
             size = 2.05, label.size = .35, fill = "white", lineheight = .95,
             fontface = "bold", label.r = unit(.08, "lines")) +
  geom_label(data = boxes[type == "c"], aes(x, y, label = lab, fill = col),
             size = 2.05, colour = "white", label.size = 0, lineheight = .95,
             fontface = "bold", label.r = unit(.08, "lines")) +
  scale_colour_identity() + scale_fill_identity() +
  annotate("text", x = 2.1, y = 0.32,
           label = "Hypothesis: sharing follows tissue compartment, not disease label",
           size = 1.95, fontface = "italic", colour = "grey25") +
  coord_cartesian(xlim = c(0, 4.3), ylim = c(0.15, 3.6), clip = "off") +
  theme_void(base_size = BASE) +
  labs(title = "Compartment model")

## --- 1b 队列构成 ----------------------------------------------------------
t1 <- melt(tab1, id.vars = c("Cohort","Disease","Tissue","Platform"),
           measure.vars = c("Control","Case"), variable.name = "grp", value.name = "n")
t1[, comp := fifelse(Tissue == "monocyte", "Haematopoietic", "Mesenchymal")]
t1[, lab := sprintf("%s\n%s | %s", Cohort, Tissue, Platform)]
t1[, lab := factor(lab, levels = unique(lab[order(Disease, -n)]))]
f1b <- ggplot(t1, aes(x = lab, y = n, fill = interaction(Disease, grp))) +
  geom_col(width = .68, colour = "grey20", linewidth = .25) +
  scale_fill_manual(values = c("OA.Control" = "#F5B7B1", "OA.Case" = COL_OA,
                               "OP.Control" = "#AED6F1", "OP.Case" = COL_OP),
                    labels = c("OA control","OP control","OA case","OP case"),
                    name = NULL) +
  facet_grid(~ comp, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Samples (n)", title = "Study cohorts") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = BASE - 2, lineheight = .85),
        legend.position = "top")

## --- 1c 队列间差异特征相关性热图 -------------------------------------------
pde[, z := qnorm(P/2, lower.tail = FALSE) * sign(delta)]
W <- dcast(pde, pathway ~ cohort, value.var = "z")
M <- as.matrix(W[, -1]); rownames(M) <- W$pathway
CM <- cor(M, use = "pairwise.complete.obs", method = "spearman")
ord <- c("GSE114007","GSE57218","GSE117999","GSE169077","GSE35958","GSE56815","GSE7158")
CM <- CM[ord, ord]
cdt <- as.data.table(as.table(CM)); setnames(cdt, c("A","B","rho"))
cdt[, `:=`(A = factor(A, ord), B = factor(B, rev(ord)))]
ann <- data.table(g = ord,
                  cls = c(rep("OA cartilage",4), "OP BM-MSC", rep("OP monocyte",2)))
f1c <- ggplot(cdt, aes(A, B, fill = rho)) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 1.75,
            colour = ifelse(abs(cdt$rho) > .55, "white", "grey15")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-.4, 1), name = expression(rho[Spearman])) +
  labs(x = NULL, y = NULL, title = "Concordance of pathway-level disease effects",
       subtitle = "Spearman correlation of per-cohort pathway z-scores") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.grid = element_blank())

FIG1 <- (f1a | f1b) / f1c +
  plot_layout(heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "a") & tag_th
save_fig(FIG1, "Fig1_design_and_landscape", 7.2, 6.2)

###############################################################################
## FIGURE 2  基因层面
###############################################################################
volc <- function(d, ttl, col) {
  d <- copy(d)[is.finite(meta_Z)]
  d[, sig := fifelse(meta_FDR < .05 & abs(mean_logFC) > .5, "sig", "ns")]
  top <- rbind(head(d[order(-meta_Z)], 8), head(d[order(meta_Z)], 8))
  ggplot(d, aes(mean_logFC, -log10(meta_P))) +
    geom_point(aes(colour = sig), size = .35, alpha = .55) +
    scale_colour_manual(values = c(sig = col, ns = "grey78"), guide = "none") +
    geom_hline(yintercept = -log10(max(d[meta_FDR < .05, meta_P], na.rm = TRUE)),
               linetype = 2, linewidth = .3, colour = "grey35") +
    ggrepel::geom_text_repel(data = top, aes(label = gene), size = 1.7,
                             max.overlaps = 30, segment.size = .2,
                             min.segment.length = 0, box.padding = .18) +
    labs(x = expression(mean~log[2]~FC), y = expression(-log[10]~italic(P)[meta]),
         title = ttl)
}
f2a <- volc(deOA, "OA cartilage (4 cohorts)", COL_OA)
f2b <- volc(deOP, "OP (3 cohorts)", COL_OP)

mg <- merge(deOA[, .(gene, zOA = meta_Z, fOA = meta_FDR)],
            deOP[, .(gene, zOP = meta_Z, fOP = meta_FDR)], by = "gene")
rr <- cor.test(mg$zOA, mg$zOP, method = "spearman")
nb <- mg[fOA < .05 & fOP < .05, .N]
f2c <- ggplot(mg, aes(zOA, zOP)) +
  geom_hline(yintercept = 0, linewidth = .25, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = .25, colour = "grey70") +
  geom_point(size = .3, alpha = .28, colour = "grey35") +
  geom_point(data = mg[fOA < .05 & fOP < .05], size = .9, colour = "#8E44AD") +
  ggrepel::geom_text_repel(data = mg[fOA < .05 & fOP < .05], aes(label = gene),
                           size = 1.7, colour = "#8E44AD", max.overlaps = 20,
                           segment.size = .2, min.segment.length = 0) +
  annotate("text", x = -Inf, y = Inf, hjust = -.08, vjust = 1.6, size = 2.05,
           label = sprintf("rho == %.3f", rr$estimate), parse = TRUE) +
  annotate("text", x = -Inf, y = Inf, hjust = -.06, vjust = 3.4, size = 2.05,
           label = sprintf("%d genes FDR<0.05 in both", nb)) +
  labs(x = expression(OA~meta~italic(Z)), y = expression(OP~meta~italic(Z)),
       title = "Gene-level concordance between diseases")

## 2d 重叠计数
ov <- data.table(
  set = c("OA only","OP only","Both"),
  n = c(deOA[meta_FDR < .05, .N] - nb, deOP[meta_FDR < .05, .N] - nb, nb))
ov[, set := factor(set, levels = c("OA only","Both","OP only"))]
f2d <- ggplot(ov, aes(set, n, fill = set)) +
  geom_col(width = .6, colour = "grey20", linewidth = .25) +
  geom_text(aes(label = comma(n)), vjust = -.4, size = 2.1) +
  scale_fill_manual(values = c("OA only" = COL_OA, "Both" = "#8E44AD",
                               "OP only" = COL_OP), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .18)), labels = comma) +
  labs(x = NULL, y = "Genes (FDR < 0.05)", title = "Overlap of differential genes")

## 2e WGCNA 模块对应
wov[, `:=`(OA = factor(OA), OP = factor(OP))]
f2e <- ggplot(wov, aes(OP, OA, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = n_shared), size = 1.8, colour = "grey15") +
  scale_fill_gradient(low = "white", high = "#7D3C98", name = "Jaccard",
                      limits = c(0, max(wov$jaccard))) +
  labs(x = "OP module", y = "OA module",
       title = "WGCNA cross-disease module correspondence",
       subtitle = sprintf("all FDR >= %.2f; max Jaccard = %.3f",
                          min(wov$FDR), max(wov$jaccard))) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank())

FIG2 <- (f2a | f2b) / (f2c | f2d) / f2e +
  plot_layout(heights = c(1, 1, 1.05)) +
  plot_annotation(tag_levels = "a") & tag_th
save_fig(FIG2, "Fig2_gene_level_no_overlap", 7.2, 8.4)

###############################################################################
## FIGURE 3  通路层面分区特异收敛
###############################################################################
cc <- copy(ctrl)
cc[, short := c("OA vs OP\n(pooled)",
                "OA RNA-seq vs\nOA array",
                "GSE114007 vs\nGSE57218",
                "GSE57218 vs\nGSE169077",
                "OP monocyte vs\nOP BM-MSC",
                "OA cartilage vs\nOP BM-MSC",
                "OA cartilage vs\nOP monocyte")]
cc[, kind := c("Main","Positive ctrl","Positive ctrl","Positive ctrl",
               "Within-disease\ncross-tissue","Mesenchymal pair","Cross-compartment")]
cc[, short := factor(short, levels = short[order(-enrichment)])]
f3a <- ggplot(cc, aes(short, enrichment, fill = kind)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = .35, colour = "grey40") +
  geom_col(width = .62, colour = "grey20", linewidth = .25) +
  geom_text(aes(label = ifelse(perm_P < .001, "P<0.001", sprintf("P=%.3f", perm_P))),
            vjust = -.45, size = 1.85) +
  scale_fill_manual(values = c("Main" = COL_NEU, "Positive ctrl" = "#5D6D7E",
                               "Within-disease\ncross-tissue" = "#AAB7B8",
                               "Mesenchymal pair" = COL_MES,
                               "Cross-compartment" = COL_HEM), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .2))) +
  labs(x = NULL, y = "Convergence enrichment\n(observed / permuted)",
       title = "Pathway-level convergence with controls") +
  theme(axis.text.x = element_text(angle = 32, hjust = 1, size = BASE - 2, lineheight = .85),
        legend.position = "right")

## 3b/3c 散点
scat <- function(mA, mB, labB, col, ttl) {
  d <- merge(mA[, .(pathway, zA = z, fA = FDR)], mB[, .(pathway, zB = z, fB = FDR)],
             by = "pathway")
  rho <- cor(d$zA, d$zB, method = "spearman")
  hit <- d[fA < .05 & fB < .25 & sign(zA) == sign(zB)]
  ggplot(d, aes(zA, zB)) +
    geom_hline(yintercept = 0, linewidth = .25, colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = .25, colour = "grey70") +
    geom_point(size = .28, alpha = .2, colour = "grey45") +
    geom_point(data = hit, size = .85, colour = col, alpha = .9) +
    geom_smooth(method = "lm", se = FALSE, linewidth = .4, colour = col, formula = y ~ x) +
    annotate("text", x = -Inf, y = Inf, hjust = -.09, vjust = 1.7, size = 2.05,
             label = sprintf("rho == %.3f", rho), parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, hjust = -.07, vjust = 3.5, size = 2.05,
             label = sprintf("n = %d concordant", nrow(hit))) +
    labs(x = expression(OA~cartilage~italic(Z)), y = labB, title = ttl)
}
f3b <- scat(m_oa, m_msc, expression(OP~BM-MSC~italic(Z)), COL_MES,
            "Mesenchymal pair (same compartment)")
f3c <- scat(m_oa, m_mono, expression(OP~monocyte~italic(Z)), COL_HEM,
            "Cross-compartment pair")

## 3d 家族层面
fm <- copy(fam)
fm[, short := c("OA array vs\nOA RNA-seq","OA cartilage vs\nOP BM-MSC",
                "OA cartilage vs\nOP monocyte","OA vs OP\n(pooled)")]
fm[, short := factor(short, levels = short[order(-fold)])]
f3d <- ggplot(fm, aes(short, fold,
                      fill = c("Positive ctrl","Mesenchymal pair",
                               "Cross-compartment","Main"))) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = .35, colour = "grey40") +
  geom_col(width = .6, colour = "grey20", linewidth = .25) +
  geom_text(aes(label = sprintf("P=%.4g\n%d/%.1f", P_perm, observed, expected)),
            vjust = -.28, size = 1.75, lineheight = .9) +
  scale_fill_manual(values = c("Main" = COL_NEU, "Positive ctrl" = "#5D6D7E",
                               "Mesenchymal pair" = COL_MES,
                               "Cross-compartment" = COL_HEM), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .28))) +
  labs(x = NULL, y = "Family-level enrichment",
       title = "After collapsing redundant pathways into families",
       subtitle = sprintf("%d pathways -> %d gene-set families (Jaccard clustering)",
                          nrow(fread('03_results/gsva/pathway_families.csv')),
                          max(fread('03_results/gsva/pathway_families.csv')$family))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = BASE - 2, lineheight = .85))

FIG3 <- (f3a) / (f3b | f3c) / f3d +
  plot_layout(heights = c(1.05, 1, 1.05)) +
  plot_annotation(tag_levels = "a") & tag_th
save_fig(FIG3, "Fig3_compartment_convergence", 7.2, 8.2)

saveRDS(list(m_oa = m_oa, m_msc = m_msc, m_mono = m_mono), "03_results/intermediate/pathway_meta_sets.rds")
cat("FIG 1-3 done\n")
