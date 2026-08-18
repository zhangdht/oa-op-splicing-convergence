###############################################################################
# 05f_splicing_axis.R
# 正式检验：RNA 剪接程序是否为 OP-OA 间充质共享轴的特征性成分
#
# 观察到 A 轴（OA 软骨 ∩ OP BM-MSC 同向）87 条通路中 13 条与 RNA 剪接/mRNA
# 代谢相关，且全部下调。本脚本用正式统计检验该富集，并与对照集比较：
#   (1) 超几何富集检验：A 轴 vs 全部被检通路背景
#   (2) 特异性对照：同样检验 OA 独有集、OP-monocyte 集、B 拮抗轴
#   (3) 样本水平剪接程序评分：在 7 个队列逐一检验（含效应量）
#   (4) SON 在该框架下的定位（回应预设假说）
#   (5) 剪接评分 vs 成骨/软骨分化评分的关系
#
# 输出：03_results/splicing/*.csv, 05_tables/Table3_splicing_axis.csv
###############################################################################

suppressPackageStartupMessages({ library(data.table); library(limma) })
ROOT <- "D:/bioinfo05"; setwd(ROOT)
dir.create("03_results/splicing", recursive = TRUE, showWarnings = FALSE)
set.seed(20260801)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

de_all  <- fread("03_results/gsva/pathway_DE_percohort.csv")
cohorts <- readRDS("03_results/intermediate/cohorts.rds")

meta_set <- function(dt) {
  dt <- dt[is.finite(P) & P > 0 & P < 1]
  dt[, z := qnorm(P / 2, lower.tail = FALSE) * sign(delta)]
  dt[, w := sqrt(n)]
  out <- dt[, .(z = sum(z * w) / sqrt(sum(w^2))), by = pathway]
  out[, FDR := p.adjust(2 * pnorm(-abs(z)), "BH")][]
}
OA <- c("GSE114007","GSE57218","GSE117999","GSE169077")
MSC <- "GSE35958"; MONO <- c("GSE56815","GSE7158")
m_oa <- meta_set(de_all[cohort %in% OA])
m_msc <- meta_set(de_all[cohort %in% MSC])
m_mono <- meta_set(de_all[cohort %in% MONO])

## ------------------------------------- 定义"剪接/mRNA 代谢"通路类别 --------
SPLICE_RX <- "SPLIC|SPLICEOSOME|MRNA_PROCESS|MRNA_METABOL|MRNA_CATABOL|MRNA_TRANSPORT|NUCLEAR_SPECK|RNA_SPLIC|SNRNP|EXON_JUNCTION"
universe <- intersect(intersect(m_oa$pathway, m_msc$pathway), m_mono$pathway)
is_spl <- grepl(SPLICE_RX, universe)
msg("背景通路 ", length(universe), " 条，其中剪接/mRNA 类 ", sum(is_spl), " 条 (",
    round(100*mean(is_spl),1), "%)")

## ------------------------------------------------- (1)(2) 富集检验 --------
A <- merge(m_oa[, .(pathway, z_OA = z, f_OA = FDR)],
           m_msc[, .(pathway, z_MSC = z, f_MSC = FDR)], by = "pathway")
B <- merge(m_oa[, .(pathway, z_OA = z, f_OA = FDR)],
           m_mono[, .(pathway, z_MO = z, f_MO = FDR)], by = "pathway")

sets <- list(
  `A: mesenchymal shared (OA+MSC, same dir)` =
    A[sign(z_OA)==sign(z_MSC) & f_OA<0.05 & f_MSC<0.25, pathway],
  `B: haematopoietic antagonistic (OA vs MONO)` =
    B[sign(z_OA)!=sign(z_MO) & f_OA<0.05 & f_MO<0.05, pathway],
  `OA significant only` = setdiff(m_oa[FDR<0.05, pathway], m_msc[FDR<0.25, pathway]),
  `OP BM-MSC significant only` = setdiff(m_msc[FDR<0.05, pathway], m_oa[FDR<0.05, pathway]),
  `OP monocyte significant` = m_mono[FDR<0.05, pathway]
)

enr <- rbindlist(lapply(names(sets), function(nm) {
  s <- intersect(sets[[nm]], universe)
  k <- sum(grepl(SPLICE_RX, s)); n <- length(s)
  K <- sum(is_spl); N <- length(universe)
  p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  data.table(set = nm, n_set = n, n_splice = k,
             pct_splice = round(100 * k / max(n,1), 1),
             expected = round(n * K / N, 1),
             fold = round((k / max(n,1)) / (K / N), 2),
             P_hyper = signif(p, 3))
}))
enr[, FDR := signif(p.adjust(P_hyper, "BH"), 3)]
fwrite(enr, "03_results/splicing/splicing_enrichment_by_set.csv")
msg("--- 剪接通路富集检验 ---"); print(enr)

# A 轴中剪接通路的方向一致性
Aspl <- A[sign(z_OA)==sign(z_MSC) & f_OA<0.05 & f_MSC<0.25 & grepl(SPLICE_RX, pathway)]
msg("A 轴剪接通路 ", nrow(Aspl), " 条，其中下调(两病均降) ",
    sum(Aspl$z_OA < 0 & Aspl$z_MSC < 0), " 条")
fwrite(Aspl[order(z_OA)], "03_results/splicing/A_axis_splicing_pathways.csv")

## ------------------------------------ (3) 样本水平剪接程序评分逐队列检验 ----
gmt <- function(f) { ln <- readLines(f, warn=FALSE); ln <- ln[nzchar(ln)]
  sp <- strsplit(ln,"\t"); setNames(lapply(sp, function(x) unique(x[-c(1,2)])), sapply(sp,`[`,1)) }
gs <- c(gmt("01_data/genesets/c5.go.bp.v7.5.1.symbols.gmt"),
        gmt("01_data/genesets/c2.cp.reactome.v7.5.1.symbols.gmt"),
        gmt("01_data/genesets/c2.cp.kegg.v7.5.1.symbols.gmt"))
spl_genes <- unique(unlist(gs[grepl(SPLICE_RX, names(gs))]))
msg("剪接程序基因集合并后 ", length(spl_genes), " 个基因")

# 对照基因集：成骨分化、软骨发育、核糖体（技术对照）
ctrl_sets <- list(
  splicing   = spl_genes,
  osteoblast = unique(unlist(gs[grepl("OSSIFICATION|OSTEOBLAST_DIFFERENTIATION|BONE_MINERALIZATION", names(gs))])),
  chondro    = unique(unlist(gs[grepl("CARTILAGE_DEVELOPMENT|CHONDROCYTE_DIFFERENTIATION", names(gs))])),
  ribosome   = unique(unlist(gs[grepl("^KEGG_RIBOSOME$|CYTOPLASMIC_TRANSLATION", names(gs))]))
)

score_one <- function(e, genes) {
  g <- intersect(genes, rownames(e)); if (length(g) < 10) return(NULL)
  z <- t(scale(t(e[g, , drop = FALSE])))          # 基因内 z 标准化
  colMeans(z, na.rm = TRUE)
}

sc_res <- rbindlist(lapply(names(cohorts), function(cn) {
  co <- cohorts[[cn]]; e <- as.matrix(co$expr)
  ph <- co$pheno[match(colnames(e), co$pheno$sample), ]
  rbindlist(lapply(names(ctrl_sets), function(sn) {
    s <- score_one(e, ctrl_sets[[sn]]); if (is.null(s)) return(NULL)
    g <- ph$group
    tt <- t.test(s[g=="Case"], s[g=="Control"])
    # Hedges' g
    n1 <- sum(g=="Case"); n2 <- sum(g=="Control")
    sp <- sqrt(((n1-1)*var(s[g=="Case"]) + (n2-1)*var(s[g=="Control"]))/(n1+n2-2))
    d  <- (mean(s[g=="Case"]) - mean(s[g=="Control"]))/sp
    J  <- 1 - 3/(4*(n1+n2)-9)
    data.table(cohort = cn, disease = co$disease, tissue = co$tissue,
               programme = sn, n_gene = length(intersect(ctrl_sets[[sn]], rownames(e))),
               delta = unname(diff(rev(tt$estimate))), t = unname(tt$statistic),
               P = tt$p.value, hedges_g = d*J)
  }))
}))
sc_res[, FDR := p.adjust(P, "BH"), by = programme]
fwrite(sc_res, "03_results/splicing/programme_scores_percohort.csv")
msg("--- 各队列程序评分（Case vs Control）---")
print(sc_res[programme == "splicing",
             .(cohort, disease, tissue, g = round(hedges_g,2),
               P = signif(P,3), FDR = signif(FDR,3))])

# 病种内 meta（Stouffer）
meta_prog <- sc_res[, {
  z <- qnorm(P/2, lower.tail = FALSE) * sign(hedges_g); w <- sqrt(n_gene*0+1)
  zz <- sum(z)/sqrt(.N)
  .(k = .N, z_meta = zz, P_meta = 2*pnorm(-abs(zz)), g_mean = mean(hedges_g))
}, by = .(programme, disease)]
fwrite(meta_prog, "03_results/splicing/programme_meta_by_disease.csv")
msg("--- 程序评分病种内 meta ---"); print(meta_prog[order(programme, disease)])

## ------------------------------------------------ (4) SON 在此框架的定位 ---
son <- rbindlist(lapply(names(cohorts), function(cn) {
  co <- cohorts[[cn]]; e <- as.matrix(co$expr)
  if (!"SON" %in% rownames(e)) return(NULL)
  ph <- co$pheno[match(colnames(e), co$pheno$sample), ]
  s <- score_one(e, spl_genes)
  v <- e["SON", ]
  tt <- t.test(v[ph$group=="Case"], v[ph$group=="Control"])
  ct <- cor.test(v, s, method = "spearman")
  data.table(cohort = cn, disease = co$disease, tissue = co$tissue,
             SON_logFC = unname(diff(rev(tt$estimate))), SON_P = tt$p.value,
             SON_vs_splicing_rho = unname(ct$estimate), rho_P = ct$p.value)
}))
fwrite(son, "03_results/splicing/SON_in_splicing_axis.csv")
msg("--- SON 与剪接程序 ---"); print(son[, .(cohort, disease,
    SON_logFC = round(SON_logFC,3), SON_P = signif(SON_P,3),
    rho = round(SON_vs_splicing_rho,3), rho_P = signif(rho_P,3))])

## ------------------------------------------------------------ Table 3 -----
tb3 <- Aspl[order(z_OA)][, .(
  Pathway = gsub("_"," ", sub("^(GOBP|KEGG|REACTOME|HALLMARK)_","",pathway)),
  Source  = sub("_.*","",pathway),
  Z_OA_cartilage = round(z_OA,2), FDR_OA = signif(f_OA,3),
  Z_OP_BMMSC = round(z_MSC,2), FDR_OP = signif(f_MSC,3),
  Direction = ifelse(z_OA<0,"Down in both","Up in both"))]
fwrite(tb3, "05_tables/Table3_splicing_axis.csv")

saveRDS(list(enrichment = enr, A_splicing = Aspl, scores = sc_res,
             meta_prog = meta_prog, son = son, spl_genes = spl_genes),
        "03_results/intermediate/splicing_axis.rds")
msg("05f 完成")
