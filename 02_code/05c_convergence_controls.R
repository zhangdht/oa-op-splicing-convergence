###############################################################################
# 05c_convergence_controls.R
# 通路收敛分析的阳性 / 阴性对照 —— 回答"无重叠是不是组织与平台差异造成的假阴性"
#
# 背景：05b 显示 OP 与 OA 之间通路收敛数(6) 低于置换零期望(9.5)，P=0.923。
# 但 OP 队列来自单核细胞 / BM-MSC，OA 队列来自关节软骨，组织与平台均不同。
# 若流程本身对跨平台、跨组织比较不敏感，则"无收敛"可能只是检验乏力的假象。
#
# 因此构建三类对照，全部使用与 05b 完全相同的统计流程：
#   [P1] 病种内跨队列跨平台（OA 内部两两）—— 应显示强收敛（阳性对照）
#   [P2] 病种内跨组织（OP: monocyte GSE56815/GSE7158 vs BM-MSC GSE35958）
#        —— 检验流程能否跨组织侦测同病收敛（关键阳性对照）
#   [N1] 随机标签置换（打乱病种归属）—— 阴性对照
#
# 输出：03_results/gsva/convergence_controls.csv
#       03_results/gsva/convergence_controls_summary.txt
###############################################################################

suppressPackageStartupMessages({ library(data.table); library(limma) })

ROOT <- "D:/bioinfo05"; setwd(ROOT)
set.seed(20260801)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

cohorts <- readRDS("03_results/intermediate/cohorts.rds")
de_all  <- fread("03_results/gsva/pathway_DE_percohort.csv")

## 与 05b 完全一致的 Stouffer meta 与收敛统计 -------------------------------
meta_set <- function(dt) {
  dt <- dt[is.finite(P) & P > 0 & P < 1]
  dt[, z := qnorm(P / 2, lower.tail = FALSE) * sign(delta)]
  dt[, w := sqrt(n)]
  out <- dt[, .(z_meta = sum(z * w) / sqrt(sum(w^2))), by = pathway]
  out[, P_meta := 2 * pnorm(-abs(z_meta))]
  out[, FDR_meta := p.adjust(P_meta, "BH")]
  out[]
}

# 给定两组队列名，返回收敛通路数、置换期望、置换 P
convergence <- function(setA, setB, label, nperm = 2000) {
  a <- meta_set(de_all[cohort %in% setA])
  b <- meta_set(de_all[cohort %in% setB])
  m <- merge(a[, .(pathway, zA = z_meta, fA = FDR_meta)],
             b[, .(pathway, zB = z_meta, fB = FDR_meta)], by = "pathway")
  if (!nrow(m)) return(NULL)
  obs <- sum(sign(m$zA) == sign(m$zB) & m$fA < 0.05 & m$fB < 0.05)
  perm <- vapply(seq_len(nperm), function(i) {
    idx <- sample.int(nrow(m))
    sum(sign(m$zA) == sign(m$zB[idx]) & m$fA < 0.05 & m$fB[idx] < 0.05)
  }, numeric(1))
  pp <- (sum(perm >= obs) + 1) / (nperm + 1)
  # 另用全通路谱的 Spearman 相关，作为不依赖阈值的连续型收敛度量
  rho <- suppressWarnings(cor(m$zA, m$zB, method = "spearman"))
  rp  <- suppressWarnings(cor.test(m$zA, m$zB, method = "spearman"))$p.value
  data.table(comparison = label,
             setA = paste(setA, collapse = "+"), setB = paste(setB, collapse = "+"),
             n_pathway = nrow(m),
             nsig_A = sum(m$fA < 0.05), nsig_B = sum(m$fB < 0.05),
             observed = obs, perm_mean = mean(perm), perm_sd = sd(perm),
             perm_P = pp, spearman_rho = rho, spearman_P = rp,
             enrichment = ifelse(mean(perm) > 0, obs / mean(perm), NA_real_))
}

res <- list()

## ------------------------------------------- 主分析（与 05b 一致，复核） ----
OA <- c("GSE114007", "GSE57218", "GSE117999", "GSE169077")
OP <- c("GSE56815", "GSE7158", "GSE35958")
res$main <- convergence(OA, OP, "MAIN: OA(all) vs OP(all)")

## ------------------------------------- [P1] 阳性对照：OA 内部跨队列跨平台 ----
# GSE114007 (RNA-seq) vs 其余三个芯片队列 —— 同病、同组织、不同平台
res$P1a <- convergence("GSE114007", c("GSE57218", "GSE117999", "GSE169077"),
                       "P1: OA RNA-seq vs OA microarray (same disease/tissue, diff platform)")
# 两两拆分
res$P1b <- convergence("GSE114007", "GSE57218", "P1b: GSE114007 vs GSE57218 (OA vs OA)")
res$P1c <- convergence("GSE57218", "GSE169077", "P1c: GSE57218 vs GSE169077 (OA vs OA)")

## ------------------------- [P2] 关键阳性对照：OP 内部跨组织（血 vs 骨髓） ----
res$P2 <- convergence(c("GSE56815", "GSE7158"), "GSE35958",
                      "P2: OP monocyte vs OP BM-MSC (same disease, DIFFERENT tissue)")

## ---------------------- [P3] 跨组织但同病的另一验证：OA 软骨 vs OP 单核 ------
# 已在 main 中；这里补充 OA vs 单独 BM-MSC，排除单一队列驱动
res$P3 <- convergence(OA, "GSE35958", "P3: OA(all) vs OP BM-MSC only")
res$P4 <- convergence(OA, c("GSE56815", "GSE7158"), "P4: OA(all) vs OP monocyte only")

## ------------------------------ [N1] 阴性对照：随机重排病种归属后的收敛 ------
allc <- c(OA, OP)
neg <- replicate(20, {
  sh <- sample(allc)
  r <- convergence(sh[1:4], sh[5:7], "neg", nperm = 300)
  if (is.null(r)) NA_real_ else r$observed / max(r$perm_mean, .Machine$double.eps)
})
msg("阴性对照（随机分组）enrichment: 均值 ", round(mean(neg, na.rm = TRUE), 3),
    " (sd ", round(sd(neg, na.rm = TRUE), 3), ")")

out <- rbindlist(Filter(Negate(is.null), res))
setcolorder(out, c("comparison", "n_pathway", "nsig_A", "nsig_B", "observed",
                   "perm_mean", "enrichment", "perm_P", "spearman_rho", "spearman_P"))
fwrite(out, "03_results/gsva/convergence_controls.csv")
print(out[, .(comparison, observed, perm_mean, enrichment = round(enrichment, 2),
              perm_P = signif(perm_P, 3), rho = round(spearman_rho, 3),
              rho_P = signif(spearman_P, 3))])

sink("03_results/gsva/convergence_controls_summary.txt")
cat("Convergence controls for OP-OA pathway-level analysis\n")
cat("Generated:", format(Sys.time()), "\n\n")
print(as.data.frame(out))
cat("\nNegative control (random disease relabelling), enrichment mean = ",
    round(mean(neg, na.rm = TRUE), 3), " sd = ", round(sd(neg, na.rm = TRUE), 3), "\n")
sink()

saveRDS(out, "03_results/intermediate/convergence_controls.rds")
msg("05c 完成")
