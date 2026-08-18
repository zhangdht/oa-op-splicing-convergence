###############################################################################
# 05e_compartment_convergence.R
# 分区特异收敛的具体通路解析
#
# 05c/05d 确立：OA 软骨 与 OP 骨髓间充质(BM-MSC) 收敛；与 OP 循环单核细胞反向。
# 本脚本回答"具体是哪些通路"，并给出三类通路清单：
#   [A] 间充质共享轴 (OA cartilage ∩ OP BM-MSC, 同向)      -> 共病的真实分子基础
#   [B] 造血拮抗轴  (OA cartilage vs OP monocyte, 反向)     -> 骨-软骨悖论的分子对应
#   [C] 疾病特异轴  (仅一病显著)                            -> 对照
#
# 输出：03_results/gsva/compartment_pathways_A_mesenchymal.csv
#       03_results/gsva/compartment_pathways_B_antagonistic.csv
#       05_tables/Table2_compartment_convergence.csv
#       03_results/intermediate/compartment.rds
###############################################################################

suppressPackageStartupMessages({ library(data.table) })
ROOT <- "D:/bioinfo05"; setwd(ROOT)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

de_all <- fread("03_results/gsva/pathway_DE_percohort.csv")

meta_set <- function(dt) {
  dt <- dt[is.finite(P) & P > 0 & P < 1]
  dt[, z := qnorm(P / 2, lower.tail = FALSE) * sign(delta)]
  dt[, w := sqrt(n)]
  out <- dt[, .(k = .N, z = sum(z * w) / sqrt(sum(w^2)), d = mean(delta)), by = pathway]
  out[, P := 2 * pnorm(-abs(z))][, FDR := p.adjust(P, "BH")]
  out[]
}

OA   <- c("GSE114007", "GSE57218", "GSE117999", "GSE169077")
MSC  <- "GSE35958"
MONO <- c("GSE56815", "GSE7158")

m_oa   <- meta_set(de_all[cohort %in% OA])
m_msc  <- meta_set(de_all[cohort %in% MSC])
m_mono <- meta_set(de_all[cohort %in% MONO])

msg("OA FDR<0.05: ", sum(m_oa$FDR < .05),
    " | BM-MSC: ", sum(m_msc$FDR < .05),
    " | monocyte: ", sum(m_mono$FDR < .05))

## ---------------------------------------------- [A] 间充质共享轴 -----------
A <- merge(m_oa[,   .(pathway, z_OA = z, FDR_OA = FDR)],
           m_msc[,  .(pathway, z_MSC = z, FDR_MSC = FDR)], by = "pathway")
A[, same_dir := sign(z_OA) == sign(z_MSC)]
A[, score := ifelse(same_dir, 2 / (1/abs(z_OA) + 1/abs(z_MSC)), 0)]
A_sig <- A[same_dir & FDR_OA < 0.05 & FDR_MSC < 0.05][order(-score)]
A_rel <- A[same_dir & FDR_OA < 0.05 & FDR_MSC < 0.25][order(-score)]   # 放宽 MSC（n=9）
msg("[A] 间充质共享轴 严格 FDR<0.05 双显著同向: ", nrow(A_sig),
    " | 放宽 MSC FDR<0.25: ", nrow(A_rel))
fwrite(A[order(-score)], "03_results/gsva/compartment_pathways_A_mesenchymal.csv")

## ---------------------------------------------- [B] 造血拮抗轴 -------------
B <- merge(m_oa[,    .(pathway, z_OA = z, FDR_OA = FDR)],
           m_mono[,  .(pathway, z_MONO = z, FDR_MONO = FDR)], by = "pathway")
B[, opposite := sign(z_OA) != sign(z_MONO)]
B[, score := ifelse(opposite, 2 / (1/abs(z_OA) + 1/abs(z_MONO)), 0)]
B_sig <- B[opposite & FDR_OA < 0.05 & FDR_MONO < 0.05][order(-score)]
msg("[B] 造血拮抗轴 双显著反向: ", nrow(B_sig))
fwrite(B[order(-score)], "03_results/gsva/compartment_pathways_B_antagonistic.csv")

## ------------------------------------------------- 通路主题归纳 ------------
themes <- list(
  `ECM & collagen`        = "COLLAGEN|EXTRACELLULAR_MATRIX|ECM_|MATRIX_ORGAN|PROTEOGLYCAN|ELASTIC_FIBRE|FIBRIL",
  `Ossification & skeletal` = "OSSIFICATION|BONE_|SKELETAL_SYSTEM|CARTILAGE|CHONDROCYTE|OSTEOBLAST|OSTEOCLAST",
  `WNT / BMP / TGFB`      = "WNT|BMP|TGF_BETA|TGFB|SMAD|HEDGEHOG|NOTCH",
  `Innate immune / myeloid` = "MYELOID|MONOCYTE|MACROPHAGE|NEUTROPHIL|GRANULOCYTE|INFLAMMAT|INTERLEUKIN|TOLL_LIKE|COMPLEMENT|CHEMOKINE|TNF",
  `Oxidative & mitochondrial` = "OXIDATIVE_PHOSPHOR|MITOCHONDRI|RESPIRATORY_ELECTRON|REACTIVE_OXYGEN|TCA|ATP_SYNTH",
  `Cell cycle & senescence` = "CELL_CYCLE|MITOTIC|DNA_REPLICATION|SENESCEN|APOPTO|P53",
  `RNA processing & splicing` = "SPLIC|SPLICEOSOME|MRNA_PROCESS|RNA_METABOL|NUCLEAR_SPECK",
  `Lipid & metabolism`    = "LIPID|FATTY_ACID|CHOLESTEROL|GLYCOLY|STEROL"
)
tag_theme <- function(p) {
  hits <- names(themes)[vapply(themes, function(rx) grepl(rx, p), logical(1))]
  if (!length(hits)) "Other" else paste(hits, collapse = "; ")
}

summarise_theme <- function(dt, label) {
  if (!nrow(dt)) return(data.table(Axis = label, Theme = NA_character_, n = 0L))
  th <- vapply(dt$pathway, tag_theme, "")
  tb <- as.data.table(table(unlist(strsplit(th, "; "))))[order(-N)]
  setnames(tb, c("Theme", "n")); tb[, Axis := label][]
}

th_A <- summarise_theme(A_rel, "A: Mesenchymal shared (OA cartilage + OP BM-MSC)")
th_B <- summarise_theme(B_sig, "B: Haematopoietic antagonistic (OA cartilage vs OP monocyte)")
msg("--- 主题分布 [A] ---"); print(th_A)
msg("--- 主题分布 [B] ---"); print(th_B)

## ------------------------------------------------------- Table 2 ----------
mk <- function(dt, cols, axis, n = 20) {
  if (!nrow(dt)) return(NULL)
  d <- head(dt, n)
  data.table(Axis = axis,
             Pathway = gsub("_", " ", sub("^(GOBP|KEGG|REACTOME|HALLMARK)_", "", d$pathway)),
             Source = sub("_.*", "", d$pathway),
             Z_OA = round(d$z_OA, 2), FDR_OA = signif(d$FDR_OA, 3),
             Z_other = round(d[[cols[1]]], 2), FDR_other = signif(d[[cols[2]]], 3),
             Theme = vapply(d$pathway, tag_theme, ""))
}
tb2 <- rbindlist(list(
  mk(A_rel, c("z_MSC", "FDR_MSC"),  "A. Mesenchymal shared"),
  mk(B_sig, c("z_MONO", "FDR_MONO"), "B. Haematopoietic antagonistic")), fill = TRUE)
fwrite(tb2, "05_tables/Table2_compartment_convergence.csv")
fwrite(rbindlist(list(th_A, th_B), fill = TRUE),
       "03_results/gsva/compartment_theme_summary.csv")

saveRDS(list(m_oa = m_oa, m_msc = m_msc, m_mono = m_mono,
             A = A, B = B, A_sig = A_sig, A_rel = A_rel, B_sig = B_sig,
             theme_A = th_A, theme_B = th_B),
        "03_results/intermediate/compartment.rds")

msg("--- [A] Top 15 间充质共享通路 ---")
if (nrow(A_rel)) print(head(A_rel[, .(pathway, z_OA = round(z_OA,2), z_MSC = round(z_MSC,2),
                                      FDR_OA = signif(FDR_OA,2), FDR_MSC = signif(FDR_MSC,2))], 15))
msg("--- [B] Top 15 造血拮抗通路 ---")
if (nrow(B_sig)) print(head(B_sig[, .(pathway, z_OA = round(z_OA,2), z_MONO = round(z_MONO,2),
                                      FDR_OA = signif(FDR_OA,2), FDR_MONO = signif(FDR_MONO,2))], 15))
msg("05e 完成")
