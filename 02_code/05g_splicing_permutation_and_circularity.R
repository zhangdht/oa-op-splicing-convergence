###############################################################################
# 05g_splicing_permutation_and_circularity.R
#
# 两项严格性修正：
#
# [I] 剪接通路富集的经验零分布
#     05f 用超几何检验得到 17.1 倍富集 (P=1.25e-10)。但该检验假设通路彼此独立，
#     而剪接类通路共享大量基因，10 条通路实际上远不足 10 个独立观测，故超几何
#     P 值被严重高估。此处改用与 05d 相同的队列内 case/control 标签置换：
#     每次置换后重新构建 A 轴集合，统计其中剪接通路数，得到经验零分布。
#
# [II] 去除 SON-剪接程序相关性的循环论证
#     05f 中 SON 本身属于剪接基因集，因此 SON 与"剪接程序评分"的相关性天然为正。
#     此处将 SON 从基因集中剔除后重算，并与随机抽取的同规模非剪接基因作对照，
#     判断 SON 与剪接程序的关联是否强于随机基因的背景相关。
#
# 输出：03_results/splicing/splicing_enrichment_permutation.csv
#       03_results/splicing/SON_correlation_deconfounded.csv
###############################################################################

suppressPackageStartupMessages({ library(limma); library(data.table) })
ROOT <- "D:/bioinfo05"; setwd(ROOT)
set.seed(20260801)
NPERM <- 500
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

cohorts <- readRDS("03_results/intermediate/cohorts.rds")
cn_all  <- names(cohorts)
gsva_scores <- lapply(cn_all, function(cn)
  readRDS(file.path("03_results/gsva", paste0("gsva_", cn, ".rds"))))
names(gsva_scores) <- cn_all
grp_list <- lapply(cn_all, function(cn) {
  ph <- cohorts[[cn]]$pheno
  ph <- ph[match(colnames(gsva_scores[[cn]]), ph$sample), ]
  factor(ph$group, levels = c("Control", "Case"))
}); names(grp_list) <- cn_all

SPLICE_RX <- "SPLIC|SPLICEOSOME|MRNA_PROCESS|MRNA_METABOL|MRNA_CATABOL|MRNA_TRANSPORT|NUCLEAR_SPECK|RNA_SPLIC|SNRNP|EXON_JUNCTION"
OA <- c("GSE114007","GSE57218","GSE117999","GSE169077"); MSC <- "GSE35958"

run_limma <- function(cn, grp) {
  sc <- gsva_scores[[cn]]; des <- model.matrix(~ grp)
  fit <- eBayes(lmFit(sc, des))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.table(pathway = rownames(tt), delta = tt$logFC, P = tt$P.Value,
             n = ncol(sc), cohort = cn)
}
meta_set <- function(dt) {
  dt <- dt[is.finite(P) & P > 0 & P < 1]
  dt[, z := qnorm(P/2, lower.tail = FALSE) * sign(delta)][, w := sqrt(n)]
  out <- dt[, .(z = sum(z*w)/sqrt(sum(w^2))), by = pathway]
  out[, FDR := p.adjust(2*pnorm(-abs(z)), "BH")][]
}
# A 轴集合与其中的剪接通路数
axisA_stats <- function(de) {
  a <- meta_set(de[cohort %in% OA]); b <- meta_set(de[cohort %in% MSC])
  m <- merge(a[, .(pathway, zA = z, fA = FDR)], b[, .(pathway, zB = z, fB = FDR)], by = "pathway")
  s <- m[sign(zA) == sign(zB) & fA < 0.05 & fB < 0.25, pathway]
  c(n_set = length(s), n_spl = sum(grepl(SPLICE_RX, s)),
    pct = if (length(s)) 100*sum(grepl(SPLICE_RX, s))/length(s) else 0)
}

## ---------------------------------------------------------------- [I] -----
de_obs <- rbindlist(lapply(cn_all, function(cn) run_limma(cn, grp_list[[cn]])))
obs <- axisA_stats(de_obs)
msg("观测 A 轴: 通路 ", obs["n_set"], " 条, 剪接 ", obs["n_spl"],
    " 条 (", round(obs["pct"],1), "%)")

msg("开始 ", NPERM, " 次置换 ...")
nullm <- matrix(NA_real_, NPERM, 3, dimnames = list(NULL, c("n_set","n_spl","pct")))
for (i in seq_len(NPERM)) {
  de_p <- rbindlist(lapply(cn_all, function(cn)
    run_limma(cn, factor(sample(as.character(grp_list[[cn]])), levels = c("Control","Case")))))
  nullm[i, ] <- axisA_stats(de_p)
  if (i %% 100 == 0) msg("  ", i, "/", NPERM)
}

p_nspl <- (sum(nullm[,"n_spl"] >= obs["n_spl"]) + 1) / (NPERM + 1)
# 条件化检验：仅在集合规模可比的置换中比较比例，避免"集合越大越易命中"的偏倚
sel <- nullm[,"n_set"] >= 20
p_pct <- if (sum(sel) >= 20)
  (sum(nullm[sel,"pct"] >= obs["pct"]) + 1) / (sum(sel) + 1) else NA_real_

res1 <- data.table(
  statistic = c("A-axis size", "splicing pathways in A axis", "% splicing in A axis"),
  observed  = c(obs["n_set"], obs["n_spl"], round(obs["pct"],2)),
  null_mean = round(colMeans(nullm), 2),
  null_q95  = round(apply(nullm, 2, quantile, .95), 2),
  null_max  = round(apply(nullm, 2, max), 2),
  emp_P     = c(NA, signif(p_nspl, 4), signif(p_pct, 4)))
fwrite(res1, "03_results/splicing/splicing_enrichment_permutation.csv")
msg("--- [I] 剪接富集的经验零分布 ---"); print(res1)
msg("规模>=20 的置换次数: ", sum(sel), " / ", NPERM)

## --------------------------------------------------------------- [II] -----
gmt <- function(f) { ln <- readLines(f, warn=FALSE); ln <- ln[nzchar(ln)]
  sp <- strsplit(ln,"\t"); setNames(lapply(sp, function(x) unique(x[-c(1,2)])), sapply(sp,`[`,1)) }
gs <- c(gmt("01_data/genesets/c5.go.bp.v7.5.1.symbols.gmt"),
        gmt("01_data/genesets/c2.cp.reactome.v7.5.1.symbols.gmt"),
        gmt("01_data/genesets/c2.cp.kegg.v7.5.1.symbols.gmt"))
spl_genes <- unique(unlist(gs[grepl(SPLICE_RX, names(gs))]))
spl_noSON <- setdiff(spl_genes, "SON")
msg("剪接基因 ", length(spl_genes), " -> 剔除 SON 后 ", length(spl_noSON))

score_of <- function(e, genes) {
  g <- intersect(genes, rownames(e)); if (length(g) < 10) return(NULL)
  colMeans(t(scale(t(e[g, , drop = FALSE]))), na.rm = TRUE)
}

res2 <- rbindlist(lapply(cn_all, function(cn) {
  e <- as.matrix(cohorts[[cn]]$expr)
  if (!"SON" %in% rownames(e)) return(NULL)
  s_all  <- score_of(e, spl_genes)
  s_noS  <- score_of(e, spl_noSON)
  v <- e["SON", ]
  rho_all <- cor(v, s_all,  method = "spearman")
  rho_no  <- cor(v, s_noS,  method = "spearman")
  # 背景：随机抽取同规模的非剪接基因集，重复 500 次，看随机基因与之相关的分布
  pool <- setdiff(rownames(e), spl_genes)
  ng <- length(intersect(spl_noSON, rownames(e)))
  bg <- replicate(500, {
    rs <- score_of(e, sample(pool, ng)); if (is.null(rs)) NA_real_ else cor(v, rs, method="spearman")
  })
  # SON 与随机基因集的相关背景（SON 是否特别贴近剪接程序）
  p_emp <- (sum(abs(bg) >= abs(rho_no), na.rm = TRUE) + 1) / (sum(!is.na(bg)) + 1)
  data.table(cohort = cn, disease = cohorts[[cn]]$disease, tissue = cohorts[[cn]]$tissue,
             rho_with_SON_in_set = round(rho_all, 3),
             rho_SON_excluded    = round(rho_no, 3),
             bg_rho_mean = round(mean(bg, na.rm=TRUE), 3),
             bg_rho_q975 = round(quantile(abs(bg), .975, na.rm=TRUE), 3),
             emp_P = signif(p_emp, 4))
}))
fwrite(res2, "03_results/splicing/SON_correlation_deconfounded.csv")
msg("--- [II] SON 与剪接程序（剔除 SON 后 + 随机基因集背景）---")
print(res2)
msg("剔除 SON 后仍显著强于随机背景的队列数: ", sum(res2$emp_P < 0.05), " / ", nrow(res2))

saveRDS(list(perm = nullm, obs = obs, enrich = res1, son = res2),
        "03_results/intermediate/splicing_permutation.rds")
msg("05g 完成")
