###############################################################################
# 05b_gsva_convergence.R
# 通路层面的 OP-OA 共病收敛分析
#
# 逻辑：03_DEG.R 已证明两病在【单基因】层面几乎无重叠（仅 4 个基因 FDR<0.05
# 同时显著，且方向全部相反）。但共病的生物学基础未必体现在同一基因上——
# 更可能体现在【同一通路被不同基因扰动】。本脚本正式检验这一假设：
#   (1) 每个队列 → GSVA 通路活性矩阵（Hallmark + KEGG + Reactome + GO:BP）
#   (2) 每个队列 → limma 差异通路活性（配对队列用 duplicateCorrelation）
#   (3) 病种内 → Stouffer 加权 Z meta（权重 = sqrt(n)）
#   (4) 跨病种 → 一致性收敛检验（同向 & 双侧 FDR）
#
# 输出：
#   03_results/gsva/gsva_<cohort>.rds
#   03_results/gsva/pathway_meta_OA.csv / pathway_meta_OP.csv
#   03_results/gsva/pathway_convergence.csv      <- 核心结果
#   05_tables/Table2_pathway_convergence.csv
###############################################################################

suppressPackageStartupMessages({
  library(GSVA); library(limma); library(data.table)
})

ROOT <- "D:/bioinfo05"
setwd(ROOT)
dir.create("03_results/gsva", recursive = TRUE, showWarnings = FALSE)
set.seed(20260801)

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

cohorts <- readRDS("03_results/intermediate/cohorts.rds")
msg("载入队列: ", paste(names(cohorts), collapse = ", "))

## ---------------------------------------------------------------- 基因集 ----
gmt_files <- c(
  HALLMARK = "01_data/genesets/h.all.v7.5.1.symbols.gmt",
  KEGG     = "01_data/genesets/c2.cp.kegg.v7.5.1.symbols.gmt",
  REACTOME = "01_data/genesets/c2.cp.reactome.v7.5.1.symbols.gmt",
  GOBP     = "01_data/genesets/c5.go.bp.v7.5.1.symbols.gmt"
)

read_gmt <- function(f) {
  ln <- readLines(f, warn = FALSE)
  ln <- ln[nzchar(ln)]
  sp <- strsplit(ln, "\t")
  setNames(lapply(sp, function(x) unique(x[-c(1, 2)])), vapply(sp, `[`, "", 1))
}

gs_all <- list()
for (nm in names(gmt_files)) {
  g <- read_gmt(gmt_files[[nm]])
  msg(nm, ": ", length(g), " 条通路")
  gs_all <- c(gs_all, g)
}
# GO:BP 太多，限制基因集规模以控制多重检验负担与噪声
sz <- lengths(gs_all)
gs_all <- gs_all[sz >= 15 & sz <= 300]
msg("过滤后 (15<=size<=300): ", length(gs_all), " 条通路")

## ------------------------------------------------------------ GSVA 打分 ----
gsva_scores <- list()
for (cn in names(cohorts)) {
  fout <- file.path("03_results/gsva", paste0("gsva_", cn, ".rds"))
  if (file.exists(fout)) {
    gsva_scores[[cn]] <- readRDS(fout); msg(cn, " 复用缓存"); next
  }
  e <- cohorts[[cn]]$expr
  e <- as.matrix(e)
  # 去掉全 NA / 零方差
  v <- matrixStats::rowVars(e, na.rm = TRUE)
  e <- e[is.finite(v) & v > 0, , drop = FALSE]
  msg(cn, " GSVA 输入 ", nrow(e), " genes x ", ncol(e), " samples")

  par <- GSVA::gsvaParam(exprData = e, geneSets = gs_all,
                         minSize = 15, maxSize = 300, kcdf = "Gaussian")
  sc <- GSVA::gsva(par, verbose = FALSE)
  saveRDS(sc, fout)
  gsva_scores[[cn]] <- sc
  msg(cn, " -> ", nrow(sc), " 条通路打分完成")
}

## --------------------------------------------------- 每队列差异通路活性 ----
de_path <- list()
for (cn in names(cohorts)) {
  sc <- gsva_scores[[cn]]
  ph <- cohorts[[cn]]$pheno
  ph <- ph[match(colnames(sc), ph$sample), , drop = FALSE]
  grp <- factor(ph$group, levels = c("Control", "Case"))
  des <- model.matrix(~ grp)

  paired <- !is.null(ph$donor) && any(duplicated(na.omit(ph$donor)))
  if (paired) {
    cf <- limma::duplicateCorrelation(sc, des, block = ph$donor)
    fit <- limma::lmFit(sc, des, block = ph$donor, correlation = cf$consensus)
  } else {
    fit <- limma::lmFit(sc, des)
  }
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  de_path[[cn]] <- data.table(pathway = rownames(tt),
                              delta = tt$logFC, t = tt$t,
                              P = tt$P.Value, FDR = tt$adj.P.Val,
                              n = ncol(sc), cohort = cn,
                              disease = cohorts[[cn]]$disease)
  msg(cn, " 差异通路 FDR<0.05: ", sum(tt$adj.P.Val < 0.05))
}
de_all <- rbindlist(de_path)
fwrite(de_all, "03_results/gsva/pathway_DE_percohort.csv")

## ------------------------------------------------- 病种内 Stouffer meta ----
meta_one <- function(dt) {
  dt <- dt[is.finite(P) & P > 0 & P < 1]
  dt[, z := qnorm(P / 2, lower.tail = FALSE) * sign(delta)]
  dt[, w := sqrt(n)]
  out <- dt[, .(k = .N,
                z_meta = sum(z * w) / sqrt(sum(w^2)),
                delta_mean = mean(delta),
                concord = mean(sign(delta) == sign(sum(z * w)))),
            by = pathway]
  out[, P_meta := 2 * pnorm(-abs(z_meta))]
  out[, FDR_meta := p.adjust(P_meta, "BH")]
  out[order(P_meta)]
}

meta_OA <- meta_one(de_all[disease == "OA"])
meta_OP <- meta_one(de_all[disease == "OP"])
fwrite(meta_OA, "03_results/gsva/pathway_meta_OA.csv")
fwrite(meta_OP, "03_results/gsva/pathway_meta_OP.csv")
msg("OA meta FDR<0.05: ", sum(meta_OA$FDR_meta < 0.05),
    " | OP meta FDR<0.05: ", sum(meta_OP$FDR_meta < 0.05))

## ------------------------------------------------------- 跨病种收敛检验 ----
cv <- merge(meta_OA[, .(pathway, z_OA = z_meta, d_OA = delta_mean,
                        P_OA = P_meta, FDR_OA = FDR_meta)],
            meta_OP[, .(pathway, z_OP = z_meta, d_OP = delta_mean,
                        P_OP = P_meta, FDR_OP = FDR_meta)],
            by = "pathway")

# 收敛统计量：Fisher 合并 + 同向要求
cv[, chisq := -2 * (log(P_OA) + log(P_OP))]
cv[, P_joint := pchisq(chisq, df = 4, lower.tail = FALSE)]
cv[, FDR_joint := p.adjust(P_joint, "BH")]
cv[, same_dir := sign(z_OA) == sign(z_OP)]
cv[, direction := ifelse(!same_dir, "Discordant",
                  ifelse(z_OA > 0, "Up in both", "Down in both"))]
# 共病收敛分数：同向时取两病 |z| 的调和均值，异向记 0
cv[, conv_score := ifelse(same_dir,
                          2 / (1 / abs(z_OA) + 1 / abs(z_OP)), 0)]
setorder(cv, -conv_score)
fwrite(cv, "03_results/gsva/pathway_convergence.csv")

sig <- cv[same_dir == TRUE & FDR_OA < 0.05 & FDR_OP < 0.05]
msg("双病 FDR<0.05 且同向的通路: ", nrow(sig))
if (nrow(sig)) print(head(sig[, .(pathway, direction, z_OA, z_OP,
                                  FDR_OA, FDR_OP)], 25))

## ------------------------------------- 与随机期望比较（收敛的显著性） -----
# 置换检验：打乱 OP 的通路标签，看同向且双显著的通路数的零分布
nperm <- 2000
obs <- nrow(sig)
perm <- integer(nperm)
zOA <- cv$z_OA; fOA <- cv$FDR_OA; zOP <- cv$z_OP; fOP <- cv$FDR_OP
for (i in seq_len(nperm)) {
  idx <- sample.int(nrow(cv))
  perm[i] <- sum(sign(zOA) == sign(zOP[idx]) & fOA < 0.05 & fOP[idx] < 0.05)
}
p_perm <- (sum(perm >= obs) + 1) / (nperm + 1)
msg("观测收敛通路数 = ", obs, "; 置换均值 = ", round(mean(perm), 1),
    "; 置换 P = ", signif(p_perm, 3))
writeLines(c(sprintf("observed_convergent_pathways\t%d", obs),
             sprintf("perm_mean\t%.2f", mean(perm)),
             sprintf("perm_sd\t%.2f", sd(perm)),
             sprintf("perm_P\t%.4g", p_perm)),
           "03_results/gsva/convergence_permutation.txt")

## ------------------------------------------------------------- Table 2 ----
tb2 <- head(sig[, .(Pathway = pathway,
                    Direction = direction,
                    Z_OA = round(z_OA, 2), FDR_OA = signif(FDR_OA, 3),
                    Z_OP = round(z_OP, 2), FDR_OP = signif(FDR_OP, 3),
                    Convergence = round(conv_score, 2))], 30)
fwrite(tb2, "05_tables/Table2_pathway_convergence.csv")

saveRDS(list(de_all = de_all, meta_OA = meta_OA, meta_OP = meta_OP,
             convergence = cv, perm_P = p_perm),
        "03_results/intermediate/gsva_convergence.rds")
msg("05b 完成")
