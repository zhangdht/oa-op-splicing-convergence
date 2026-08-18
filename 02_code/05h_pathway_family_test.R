###############################################################################
# 05h_pathway_family_test.R
# 通路家族层面的收敛与剪接富集检验 —— 消除通路间基因共享造成的伪独立
#
# 问题：05c 的通路标签重排检验条件于各病的边际显著集，是正确的收敛检验，
# 但它把 5,059 条通路当作独立单位。事实上 GO/Reactome 通路大量共享基因
# （如 10 条剪接通路彼此 Jaccard 很高），实际独立单位远少于名义条数，
# 因此超几何 / 标签重排 P 值均被高估。
#
# 解决：先按基因集 Jaccard 相似度将通路层次聚类成"通路家族"，
# 每个家族取一个代表（家族内 |z| 最大者），在家族层面重做：
#   (1) A 轴收敛检验（OA 软骨 ∩ OP BM-MSC 同向）
#   (2) 剪接家族是否在 A 轴中显著富集
#   (3) 与阳性对照 P1（OA vs OA）、阴性对照 P4（OA vs 单核）并列比较
#
# 输出：03_results/gsva/pathway_families.csv
#       03_results/gsva/family_level_convergence.csv
###############################################################################

suppressPackageStartupMessages({ library(data.table) })
ROOT <- "D:/bioinfo05"; setwd(ROOT)
set.seed(20260801)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

de_all <- fread("03_results/gsva/pathway_DE_percohort.csv")
meta_set <- function(dt) {
  dt <- dt[is.finite(P) & P > 0 & P < 1]
  dt[, z := qnorm(P/2, lower.tail = FALSE) * sign(delta)][, w := sqrt(n)]
  out <- dt[, .(z = sum(z*w)/sqrt(sum(w^2))), by = pathway]
  out[, FDR := p.adjust(2*pnorm(-abs(z)), "BH")][]
}
OA <- c("GSE114007","GSE57218","GSE117999","GSE169077")
MSC <- "GSE35958"; MONO <- c("GSE56815","GSE7158")
OA_RNA <- "GSE114007"; OA_ARR <- c("GSE57218","GSE117999","GSE169077")

## ------------------------------------------------------- 读取基因集 -------
gmt <- function(f) { ln <- readLines(f, warn=FALSE); ln <- ln[nzchar(ln)]
  sp <- strsplit(ln,"\t"); setNames(lapply(sp, function(x) unique(x[-c(1,2)])), sapply(sp,`[`,1)) }
gs <- c(gmt("01_data/genesets/h.all.v7.5.1.symbols.gmt"),
        gmt("01_data/genesets/c2.cp.kegg.v7.5.1.symbols.gmt"),
        gmt("01_data/genesets/c2.cp.reactome.v7.5.1.symbols.gmt"),
        gmt("01_data/genesets/c5.go.bp.v7.5.1.symbols.gmt"))

tested <- Reduce(intersect, list(meta_set(de_all[cohort %in% OA])$pathway,
                                 meta_set(de_all[cohort %in% MSC])$pathway,
                                 meta_set(de_all[cohort %in% MONO])$pathway))
gs <- gs[names(gs) %in% tested]
msg("参与聚类的通路 ", length(gs), " 条")

## ---------------------------------------- Jaccard 相似度 + 层次聚类 -------
# 用稀疏 0/1 基因-通路矩阵计算 Jaccard，避免 5000x5000 双重循环
allg <- unique(unlist(gs))
idx  <- lapply(gs, function(g) match(g, allg))
library(Matrix)
i <- unlist(idx); j <- rep(seq_along(idx), lengths(idx))
M <- sparseMatrix(i = i, j = j, x = 1, dims = c(length(allg), length(gs)))
inter <- as.matrix(Matrix::crossprod(M))              # 交集大小
sz <- diag(inter)
union_ <- outer(sz, sz, "+") - inter
J <- inter / pmax(union_, 1)
msg("Jaccard 矩阵完成: ", nrow(J), " x ", ncol(J))

hc <- hclust(as.dist(1 - J), method = "average")
CUT <- 0.7                                            # 家族内平均 Jaccard >= 0.3
fam <- cutree(hc, h = CUT)
msg("Jaccard 距离 ", CUT, " 切割 -> ", length(unique(fam)), " 个通路家族")
fam_dt <- data.table(pathway = names(gs), family = fam)
fwrite(fam_dt, "03_results/gsva/pathway_families.csv")
msg("家族规模分布："); print(summary(as.integer(table(fam))))

SPLICE_RX <- "SPLIC|SPLICEOSOME|MRNA_PROCESS|MRNA_METABOL|MRNA_CATABOL|MRNA_TRANSPORT|NUCLEAR_SPECK|RNA_SPLIC|SNRNP|EXON_JUNCTION"
fam_dt[, is_splice := grepl(SPLICE_RX, pathway)]
fam_splice <- fam_dt[, .(any_splice = any(is_splice)), by = family]
msg("剪接家族数 ", sum(fam_splice$any_splice), " / ", nrow(fam_splice))

## -------------------------------- 家族代表：|z| 最大者（按各自比较定义） ---
rep_by <- function(mz) {                # mz: pathway,z,FDR
  d <- merge(mz, fam_dt[, .(pathway, family)], by = "pathway")
  d[order(-abs(z)), .SD[1], by = family]
}

conv_family <- function(setA, setB, fdrA = 0.05, fdrB = 0.05, label, nperm = 5000) {
  a <- meta_set(de_all[cohort %in% setA]); b <- meta_set(de_all[cohort %in% setB])
  m <- merge(a[, .(pathway, zA = z, fA = FDR)], b[, .(pathway, zB = z, fB = FDR)], by = "pathway")
  m <- merge(m, fam_dt[, .(pathway, family, is_splice)], by = "pathway")
  # 家族层面：家族内取 |zA|+|zB| 最大的通路为代表
  fm <- m[order(-(abs(zA) + abs(zB))), .SD[1], by = family]
  sigA <- fm[fA < fdrA, family]; sigB <- fm[fB < fdrB, family]
  conv <- fm[family %in% intersect(sigA, sigB) & sign(zA) == sign(zB), family]
  # 条件化置换：固定两侧显著家族数，随机重排 B 侧家族标签
  N <- nrow(fm); nA <- length(sigA); nB <- length(sigB)
  perm <- replicate(nperm, {
    sb <- sample(fm$family, nB)
    length(intersect(sigA, sb)) * 0.5     # 随机同向概率约 1/2
  })
  p <- (sum(perm >= length(conv)) + 1) / (nperm + 1)
  data.table(comparison = label, n_family = N, sigA = nA, sigB = nB,
             observed = length(conv), expected = round(mean(perm), 2),
             fold = round(length(conv) / max(mean(perm), 1e-9), 2),
             P_perm = signif(p, 4),
             n_splice_conv = sum(fm[family %in% conv, is_splice]),
             conv_families = paste(head(fm[family %in% conv, pathway], 40), collapse = " | "))
}

res <- rbindlist(list(
  conv_family(OA_RNA, OA_ARR, 0.05, 0.05, "P1 positive control: OA RNA-seq vs OA microarray"),
  conv_family(OA, MSC, 0.05, 0.25, "A: OA cartilage vs OP BM-MSC (mesenchymal)"),
  conv_family(OA, MONO, 0.05, 0.05, "P4: OA cartilage vs OP monocyte (haematopoietic)"),
  conv_family(OA, c(MSC, MONO), 0.05, 0.05, "MAIN: OA vs OP pooled")
))
fwrite(res[, -"conv_families"], "03_results/gsva/family_level_convergence.csv")
msg("--- 家族层面收敛检验 ---")
print(res[, .(comparison, n_family, sigA, sigB, observed, expected, fold, P_perm, n_splice_conv)])

## ------------------------------ 剪接家族在 A 轴中的富集（家族层面）--------
a <- meta_set(de_all[cohort %in% OA]); b <- meta_set(de_all[cohort %in% MSC])
m <- merge(a[, .(pathway, zA = z, fA = FDR)], b[, .(pathway, zB = z, fB = FDR)], by = "pathway")
m <- merge(m, fam_dt[, .(pathway, family, is_splice)], by = "pathway")
fm <- m[order(-(abs(zA) + abs(zB))), .SD[1], by = family]
Aset <- fm[fA < 0.05 & fB < 0.25 & sign(zA) == sign(zB)]
K <- sum(fm$is_splice); N <- nrow(fm); n <- nrow(Aset); k <- sum(Aset$is_splice)
p_hyp <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
msg("家族层面剪接富集: A 轴 ", n, " 个家族，其中剪接家族 ", k,
    "；背景 ", K, "/", N, "；fold = ", round((k/n)/(K/N), 2),
    "；hypergeometric P = ", signif(p_hyp, 4))
spl_out <- data.table(level = "family", A_size = n, A_splice = k,
                      bg_splice = K, bg_total = N,
                      fold = round((k/max(n,1))/(K/N), 2), P = signif(p_hyp, 4))
fwrite(spl_out, "03_results/gsva/family_level_splicing_enrichment.csv")
if (k) { msg("A 轴中的剪接家族代表通路："); print(Aset[is_splice == TRUE,
      .(pathway, zA = round(zA,2), zB = round(zB,2), fA = signif(fA,3), fB = signif(fB,3))]) }

saveRDS(list(families = fam_dt, convergence = res, splicing = spl_out, Aset = Aset),
        "03_results/intermediate/pathway_family_test.rds")
msg("05h 完成")
