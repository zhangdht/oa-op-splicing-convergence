###############################################################################
# 05d_label_permutation_null.R
# 严格零分布：队列内 case/control 标签置换
#
# 05c 中使用的"通路标签重排"零分布会破坏通路间的基因共享相关结构，因此偏保守
# 且校准不佳（随机分组 enrichment 均值 1.82）。本脚本改用金标准零分布：
#   在【每个队列内部】随机置换 case/control 标签 → 重跑 limma → 重算 meta 与
#   跨组收敛统计量。此过程完整保留了通路间相关结构、队列规模、平台噪声与
#   组织特异的表达协方差，仅破坏"疾病状态"这一唯一信息。
#
# 因 GSVA 打分已缓存，仅需重跑 limma，代价可接受。
#
# 检验对象（与 05c 相同的三个关键对比）：
#   MAIN : OA(all)  vs OP(all)
#   P3   : OA(all)  vs OP BM-MSC (间充质分区)
#   P4   : OA(all)  vs OP monocyte (造血分区)
#   P1   : OA RNA-seq vs OA microarray (阳性对照)
#
# 输出：03_results/gsva/label_permutation_null.csv
#       03_results/gsva/label_permutation_null.txt
###############################################################################

suppressPackageStartupMessages({ library(limma); library(data.table) })

ROOT <- "D:/bioinfo05"; setwd(ROOT)
set.seed(20260801)
NPERM <- 500
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

cohorts <- readRDS("03_results/intermediate/cohorts.rds")
cn_all  <- names(cohorts)

## 载入缓存的 GSVA 打分 -----------------------------------------------------
gsva_scores <- lapply(cn_all, function(cn)
  readRDS(file.path("03_results/gsva", paste0("gsva_", cn, ".rds"))))
names(gsva_scores) <- cn_all
msg("载入 GSVA 缓存 ", length(gsva_scores), " 个队列")

# 预取表型与配对信息
meta_info <- lapply(cn_all, function(cn) {
  ph <- cohorts[[cn]]$pheno
  ph <- ph[match(colnames(gsva_scores[[cn]]), ph$sample), , drop = FALSE]
  list(grp = factor(ph$group, levels = c("Control", "Case")),
       donor = if (!is.null(ph$donor)) as.character(ph$donor) else NULL,
       n = ncol(gsva_scores[[cn]]))
})
names(meta_info) <- cn_all

## 单队列 limma（可传入置换后的分组）---------------------------------------
run_limma <- function(cn, grp) {
  sc <- gsva_scores[[cn]]
  des <- model.matrix(~ grp)
  fit <- limma::lmFit(sc, des)            # 置换下不再使用 block，保持零分布一致性
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.table(pathway = rownames(tt), delta = tt$logFC, P = tt$P.Value,
             n = ncol(sc), cohort = cn)
}

meta_set <- function(dt) {
  dt <- dt[is.finite(P) & P > 0 & P < 1]
  dt[, z := qnorm(P / 2, lower.tail = FALSE) * sign(delta)]
  dt[, w := sqrt(n)]
  out <- dt[, .(z_meta = sum(z * w) / sqrt(sum(w^2))), by = pathway]
  out[, FDR_meta := p.adjust(2 * pnorm(-abs(z_meta)), "BH")]
  out[]
}

conv_stat <- function(de, setA, setB) {
  a <- meta_set(de[cohort %in% setA]); b <- meta_set(de[cohort %in% setB])
  m <- merge(a[, .(pathway, zA = z_meta, fA = FDR_meta)],
             b[, .(pathway, zB = z_meta, fB = FDR_meta)], by = "pathway")
  c(nconv = sum(sign(m$zA) == sign(m$zB) & m$fA < 0.05 & m$fB < 0.05),
    rho   = suppressWarnings(cor(m$zA, m$zB, method = "spearman")))
}

OA  <- c("GSE114007", "GSE57218", "GSE117999", "GSE169077")
OP  <- c("GSE56815", "GSE7158", "GSE35958")
MSC <- "GSE35958"; MONO <- c("GSE56815", "GSE7158")
OA_RNA <- "GSE114007"; OA_ARR <- c("GSE57218", "GSE117999", "GSE169077")

comparisons <- list(
  MAIN = list(OA, OP),  P3 = list(OA, MSC),
  P4   = list(OA, MONO), P1 = list(OA_RNA, OA_ARR),
  P2   = list(MONO, MSC)
)

## 观测值 -------------------------------------------------------------------
de_obs <- rbindlist(lapply(cn_all, function(cn) run_limma(cn, meta_info[[cn]]$grp)))
obs <- t(vapply(comparisons, function(x) conv_stat(de_obs, x[[1]], x[[2]]), numeric(2)))
msg("观测值：")
print(round(obs, 4))

## 置换 ---------------------------------------------------------------------
msg("开始 ", NPERM, " 次队列内标签置换 ...")
null_n <- matrix(NA_real_, NPERM, length(comparisons),
                 dimnames = list(NULL, names(comparisons)))
null_r <- null_n
for (i in seq_len(NPERM)) {
  de_p <- rbindlist(lapply(cn_all, function(cn) {
    g <- meta_info[[cn]]$grp
    run_limma(cn, factor(sample(as.character(g)), levels = c("Control", "Case")))
  }))
  for (k in names(comparisons)) {
    s <- conv_stat(de_p, comparisons[[k]][[1]], comparisons[[k]][[2]])
    null_n[i, k] <- s["nconv"]; null_r[i, k] <- s["rho"]
  }
  if (i %% 50 == 0) msg("  ", i, "/", NPERM)
}

## 汇总 ---------------------------------------------------------------------
res <- rbindlist(lapply(names(comparisons), function(k) {
  data.table(
    comparison   = k,
    obs_nconv    = obs[k, "nconv"],
    null_nconv_mean = mean(null_n[, k]), null_nconv_q95 = quantile(null_n[, k], .95),
    P_nconv      = (sum(null_n[, k] >= obs[k, "nconv"]) + 1) / (NPERM + 1),
    obs_rho      = obs[k, "rho"],
    null_rho_mean = mean(null_r[, k]), null_rho_sd = sd(null_r[, k]),
    P_rho        = (sum(abs(null_r[, k]) >= abs(obs[k, "rho"])) + 1) / (NPERM + 1),
    Z_rho        = (obs[k, "rho"] - mean(null_r[, k])) / sd(null_r[, k]))
}))
fwrite(res, "03_results/gsva/label_permutation_null.csv")
print(res)

sink("03_results/gsva/label_permutation_null.txt")
cat("Within-cohort case/control label permutation null (", NPERM, " permutations)\n", sep = "")
cat("Generated:", format(Sys.time()), "\n\n")
print(as.data.frame(res))
cat("\nInterpretation key:\n")
cat("  MAIN = OA(all) vs OP(all)\n  P3 = OA vs OP BM-MSC (mesenchymal)\n")
cat("  P4 = OA vs OP monocyte (haematopoietic)\n  P1 = OA RNA-seq vs OA microarray (positive control)\n")
cat("  P2 = OP monocyte vs OP BM-MSC\n")
sink()

saveRDS(list(obs = obs, null_n = null_n, null_r = null_r, res = res),
        "03_results/intermediate/label_permutation_null.rds")
msg("05d 完成")
