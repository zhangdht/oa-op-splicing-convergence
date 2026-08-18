###############################################################################
# 07_ml_compartment_transfer.R
#
# 科学问题（由分区特异收敛假说导出的可证伪预测）：
#   若 OA 软骨与 OP 骨髓间充质干细胞（BM-MSC）共享同一条 RNA 剪接程序，
#   则「在 OA 软骨上训练的剪接程序分类器」应当能迁移到 OP BM-MSC，
#   却不能（甚至反向）迁移到 OP 循环单核细胞（造血分区）。
#
# 这是本研究最关键的一个 out-of-sample 检验：它不依赖于同一批通路统计量，
# 而是用完全独立的建模框架重新问同一个问题。
#
# 设计：
#   (1) 特征集：剪接程序基因、WGCNA 跨病交集基因、随机同规模基因集（经验零分布）
#   (2) 模型：LASSO logistic / Random Forest / SVM-linear
#   (3) 队列内：留一队列交叉验证（LOCO）
#   (4) 跨分区迁移：OA 软骨 -> OP BM-MSC（同分区）vs OA 软骨 -> OP 单核（跨分区）
#   (5) 反向：OP BM-MSC -> OA 软骨
#   (6) 每个基因在每个队列内独立 z 标准化，消除平台/组织基线差异
#
# 输出：03_results/ml/*.csv
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(randomForest)
  library(pROC); library(e1071)
})
ROOT <- "D:/bioinfo05"; setwd(ROOT)
dir.create("03_results/ml", recursive = TRUE, showWarnings = FALSE)
set.seed(20260801)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

cohorts <- readRDS("03_results/intermediate/cohorts.rds")
OA   <- c("GSE114007","GSE57218","GSE117999","GSE169077")
MSC  <- "GSE35958"
MONO <- c("GSE56815","GSE7158")

## ---------------------------------------------------------- 基因集 --------
splice_genes <- readLines("03_results/splice/core_spliceosome_genes.txt")
inter_genes  <- readLines("03_results/wgcna/shared_module_genes.txt")
msg("剪接程序基因 ", length(splice_genes), " 个；WGCNA 交集基因 ", length(inter_genes), " 个")

## 共同基因宇宙（7 队列都测到）
universe <- Reduce(intersect, lapply(cohorts, function(x) rownames(x$expr)))
msg("7 队列共同基因宇宙：", length(universe))

fs_splice <- intersect(splice_genes, universe)
fs_inter  <- intersect(inter_genes,  universe)
msg("可用剪接特征 ", length(fs_splice), " 个；可用交集特征 ", length(fs_inter), " 个")

## ------------------------------------------- 构建标准化后的队列矩阵 --------
# 每个基因在每个队列内 z 标准化：这样跨平台/跨组织的基线差异被完全消除，
# 分类器只能利用「病例-对照的相对差异模式」，这正是我们想检验的东西。
prep <- function(gse, genes) {
  x  <- cohorts[[gse]]
  e  <- x$expr[genes, , drop = FALSE]
  e  <- t(scale(t(e)))                      # 基因内 z（队列内）
  e[!is.finite(e)] <- 0
  list(X = t(e),                            # 样本 x 基因
       y = factor(ifelse(x$pheno$group == "Case", 1, 0), levels = c(0, 1)),
       gse = gse)
}

## ------------------------------------------------------- 建模函数 --------
fit_predict <- function(Xtr, ytr, Xte, method) {
  ok <- apply(Xtr, 2, function(v) sd(v) > 0)
  Xtr <- Xtr[, ok, drop = FALSE]; Xte <- Xte[, colnames(Xtr), drop = FALSE]
  if (ncol(Xtr) < 2) return(rep(NA_real_, nrow(Xte)))
  out <- tryCatch({
    if (method == "LASSO") {
      nf <- min(5, min(table(ytr)))
      if (nf < 3) {
        f <- glmnet(Xtr, ytr, family = "binomial", alpha = 1)
        s <- f$lambda[max(1, floor(length(f$lambda) * 0.6))]
        as.numeric(predict(f, Xte, s = s, type = "response"))
      } else {
        cvf <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, nfolds = nf)
        as.numeric(predict(cvf, Xte, s = "lambda.min", type = "response"))
      }
    } else if (method == "RF") {
      f <- randomForest(x = Xtr, y = ytr, ntree = 1000)
      as.numeric(predict(f, Xte, type = "prob")[, "1"])
    } else if (method == "SVM") {
      f <- svm(x = Xtr, y = ytr, kernel = "linear", probability = TRUE, scale = FALSE)
      p <- predict(f, Xte, probability = TRUE)
      as.numeric(attr(p, "probabilities")[, "1"])
    }
  }, error = function(e) rep(NA_real_, nrow(Xte)))
  out
}

safe_auc <- function(y, p) {
  if (all(is.na(p)) || length(unique(y)) < 2) return(c(NA, NA, NA))
  r <- tryCatch(pROC::roc(y, p, quiet = TRUE, direction = "<"), error = function(e) NULL)
  if (is.null(r)) return(c(NA, NA, NA))
  ci <- tryCatch(as.numeric(pROC::ci.auc(r)), error = function(e) c(NA, NA, NA))
  c(as.numeric(r$auc), ci[1], ci[3])
}

## -------------------------------------------------- 主实验：迁移矩阵 --------
run_scenario <- function(train_gses, test_gses, genes, label, methods = c("LASSO","RF","SVM")) {
  tr <- lapply(train_gses, prep, genes = genes)
  Xtr <- do.call(rbind, lapply(tr, `[[`, "X"))
  # 注意：unlist() 作用于 factor 列表会返回整数编码，必须先转字符再转整数
  ytr <- factor(unlist(lapply(tr, function(z) as.integer(as.character(z$y)))),
                levels = c(0, 1))
  res <- rbindlist(lapply(test_gses, function(g) {
    te <- prep(g, genes)
    rbindlist(lapply(methods, function(m) {
      p <- fit_predict(Xtr, ytr, te$X, m)
      a <- safe_auc(te$y, p)
      data.table(scenario = label, train = paste(train_gses, collapse = "+"),
                 test = g, method = m, n_test = nrow(te$X),
                 AUC = a[1], CI_lo = a[2], CI_hi = a[3])
    }))
  }))
  res
}

FEATURE_SETS <- list(Splicing = fs_splice, WGCNA_intersection = fs_inter)

all_res <- list()
for (fsname in names(FEATURE_SETS)) {
  g <- FEATURE_SETS[[fsname]]
  msg("=== 特征集：", fsname, " (", length(g), " 基因) ===")

  ## (a) OA 内部：留一队列交叉验证
  loco <- rbindlist(lapply(OA, function(ho) {
    run_scenario(setdiff(OA, ho), ho, g, "OA internal (LOCO)")
  }))

  ## (b) OA -> OP BM-MSC（同为间充质分区）
  tr_msc <- run_scenario(OA, MSC, g, "OA cartilage -> OP BM-MSC (same compartment)")

  ## (c) OA -> OP 单核（造血分区）
  tr_mono <- run_scenario(OA, MONO, g, "OA cartilage -> OP monocyte (cross compartment)")

  ## (d) 反向：OP BM-MSC -> OA
  rev_msc <- run_scenario(MSC, OA, g, "OP BM-MSC -> OA cartilage (reverse)")

  ## (e) OP 单核内部（阳性对照）
  mono_int <- run_scenario("GSE56815", "GSE7158", g, "OP monocyte internal (positive ctrl)")

  r <- rbindlist(list(loco, tr_msc, tr_mono, rev_msc, mono_int))
  r[, feature_set := fsname]
  all_res[[fsname]] <- r
}
ML <- rbindlist(all_res)
fwrite(ML, "03_results/ml/transfer_auc_full.csv")

## 汇总（按场景 x 特征集 x 模型取跨测试队列的样本量加权均值）
summ <- ML[!is.na(AUC), .(mean_AUC = sum(AUC * n_test) / sum(n_test),
                          n_tests = .N), by = .(feature_set, scenario, method)]
fwrite(summ, "03_results/ml/transfer_auc_summary.csv")
print(summ)

## ------------------------------------------ 随机基因集经验零分布 --------
# 关键：剪接特征的迁移 AUC 是否高于同规模随机基因集？
NPERM <- 200
msg("随机基因集零分布，", NPERM, " 次（LASSO，OA -> MSC / OA -> MONO）...")
nullres <- rbindlist(lapply(seq_len(NPERM), function(i) {
  gg <- sample(universe, length(fs_splice))
  a1 <- run_scenario(OA, MSC,  gg, "null", methods = "LASSO")
  a2 <- run_scenario(OA, MONO, gg, "null", methods = "LASSO")
  data.table(iter = i,
             AUC_MSC  = a1$AUC[1],
             AUC_MONO = sum(a2$AUC * a2$n_test, na.rm = TRUE) / sum(a2$n_test[!is.na(a2$AUC)]))
}), fill = TRUE)
fwrite(nullres, "03_results/ml/random_geneset_null.csv")

obs_msc  <- ML[feature_set == "Splicing" & method == "LASSO" &
                 scenario == "OA cartilage -> OP BM-MSC (same compartment)", AUC][1]
obs_mono <- ML[feature_set == "Splicing" & method == "LASSO" &
                 scenario == "OA cartilage -> OP monocyte (cross compartment)"]
obs_mono <- sum(obs_mono$AUC * obs_mono$n_test, na.rm = TRUE) /
  sum(obs_mono$n_test[!is.na(obs_mono$AUC)])

p_msc  <- (1 + sum(nullres$AUC_MSC  >= obs_msc,  na.rm = TRUE)) / (1 + sum(!is.na(nullres$AUC_MSC)))
p_mono <- (1 + sum(nullres$AUC_MONO >= obs_mono, na.rm = TRUE)) / (1 + sum(!is.na(nullres$AUC_MONO)))

txt <- c(
  "=== 剪接程序分类器的跨分区迁移 vs 随机基因集零分布 (LASSO) ===",
  sprintf("随机基因集规模 = %d，置换次数 = %d", length(fs_splice), NPERM),
  sprintf("OA -> OP BM-MSC   观测 AUC = %.3f；零分布均值 %.3f (sd %.3f)；经验 P = %.4f",
          obs_msc, mean(nullres$AUC_MSC, na.rm = TRUE), sd(nullres$AUC_MSC, na.rm = TRUE), p_msc),
  sprintf("OA -> OP monocyte 观测 AUC = %.3f；零分布均值 %.3f (sd %.3f)；经验 P = %.4f",
          obs_mono, mean(nullres$AUC_MONO, na.rm = TRUE), sd(nullres$AUC_MONO, na.rm = TRUE), p_mono),
  "",
  "解释：若同分区迁移显著优于随机、跨分区迁移不优于随机（或反向），",
  "则支持『分区特异共享』而非『疾病特异共享』。"
)
writeLines(txt, "03_results/ml/transfer_vs_null.txt")
cat(paste(txt, collapse = "\n"), "\n")

saveRDS(list(ML = ML, summ = summ, null = nullres,
             obs = c(MSC = obs_msc, MONO = obs_mono),
             p = c(MSC = p_msc, MONO = p_mono),
             fs_splice = fs_splice, fs_inter = fs_inter),
        "03_results/intermediate/ml_transfer.rds")
msg("完成 07_ml_compartment_transfer.R")
