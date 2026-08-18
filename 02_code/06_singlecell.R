###############################################################################
# 06_singlecell.R
# 单细胞层面定位剪接程序：GSE324993（人软骨，HC x3 vs OA x3，44,348 细胞，8 亚型）
#
# 目的：bulk 层面发现"剪接程序在间充质分区一致下调"。但 bulk 无法区分
#   (a) 是所有软骨细胞普遍下调，还是
#   (b) 由某一特定亚群驱动，抑或
#   (c) 只是 OA 中细胞组成改变（亚群比例漂移）造成的假象。
# 单细胞可直接区分这三种可能——这是该结论能否成立的关键检验。
#
# 分析：
#   1. 构建 Seurat 对象，沿用作者提供的 8 类软骨细胞注释
#   2. 细胞组成分析（HC vs OA 各亚型比例，含 per-sample 变异）
#   3. AddModuleScore 计算剪接程序评分（对照：成骨/软骨/核糖体）
#   4. 亚型内 HC vs OA 比较（以样本为单位做 pseudobulk，避免伪重复）
#   5. SON 表达在亚型与病例状态间的分布
#
# 输出：03_results/singlecell/*.csv, 04_figures/Fig4_*.pdf/png
###############################################################################

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
})
ROOT <- "D:/bioinfo05"; setwd(ROOT)
dir.create("03_results/singlecell", recursive = TRUE, showWarnings = FALSE)
dir.create("04_figures", showWarnings = FALSE)
set.seed(20260801)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ------------------------------------------------------------- 读入数据 ----
msg("读取表达矩阵 ...")
cnt <- fread("01_data/sc/GSE324993_counts_matrix.csv.gz", data.table = FALSE)
rownames(cnt) <- cnt[[1]]; cnt <- cnt[, -1, drop = FALSE]
cnt <- as.matrix(cnt)
msg("矩阵: ", nrow(cnt), " genes x ", ncol(cnt), " cells")

ann <- fread("01_data/sc/GSE324993_cell_annotation.txt.gz", skip = 4,
             header = TRUE, data.table = FALSE)
colnames(ann)[1:2] <- c("cell", "celltype")
ann <- ann[ann$cell %in% colnames(cnt), ]
msg("注释细胞 ", nrow(ann), " 个；类型 ", length(unique(ann$celltype)))

cnt <- cnt[, ann$cell, drop = FALSE]
meta <- data.frame(row.names = ann$cell,
                   celltype = ann$celltype,
                   sample   = sub("_.*", "", ann$cell))
meta$group <- ifelse(grepl("^HC", meta$sample), "Healthy", "OA")
meta$celltype_short <- sub("\\(.*", "", meta$celltype)
print(table(meta$group, meta$celltype_short))

## --------------------------------------------------------- Seurat 流程 ----
so <- CreateSeuratObject(counts = cnt, meta.data = meta, min.cells = 3)
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- RunUMAP(so, dims = 1:20, verbose = FALSE)
msg("Seurat 处理完成: ", ncol(so), " cells x ", nrow(so), " genes")

## ------------------------------------------------------ 程序基因集打分 ----
gmt <- function(f) { ln <- readLines(f, warn=FALSE); ln <- ln[nzchar(ln)]
  sp <- strsplit(ln,"\t"); setNames(lapply(sp, function(x) unique(x[-c(1,2)])), sapply(sp,`[`,1)) }
gs <- c(gmt("01_data/genesets/c5.go.bp.v7.5.1.symbols.gmt"),
        gmt("01_data/genesets/c2.cp.reactome.v7.5.1.symbols.gmt"),
        gmt("01_data/genesets/c2.cp.kegg.v7.5.1.symbols.gmt"))
SPLICE_RX <- "SPLIC|SPLICEOSOME|MRNA_PROCESS|MRNA_METABOL|MRNA_CATABOL|MRNA_TRANSPORT|NUCLEAR_SPECK|RNA_SPLIC|SNRNP|EXON_JUNCTION"

progs <- list(
  Splicing   = unique(unlist(gs[grepl(SPLICE_RX, names(gs))])),
  Osteogenic = unique(unlist(gs[grepl("OSSIFICATION|OSTEOBLAST_DIFFERENTIATION|BONE_MINERALIZATION", names(gs))])),
  Chondro    = unique(unlist(gs[grepl("CARTILAGE_DEVELOPMENT|CHONDROCYTE_DIFFERENTIATION", names(gs))])),
  Ribosome   = unique(unlist(gs[grepl("^KEGG_RIBOSOME$|CYTOPLASMIC_TRANSLATION", names(gs))])),
  Hypertrophy= unique(unlist(gs[grepl("CHONDROCYTE_HYPERTROPHY|ENDOCHONDRAL", names(gs))]))
)
cov <- sapply(progs, function(g) length(intersect(g, rownames(so))))
msg("各程序在本数据中的基因覆盖数：")
print(cov)
progs <- progs[cov >= 15]
msg("保留程序: ", paste(names(progs), collapse = ", "))

so <- AddModuleScore(so, features = lapply(progs, function(g) intersect(g, rownames(so))),
                     name = "PROG", seed = 1, verbose = FALSE)
sc_cols <- paste0("PROG", seq_along(progs))
names(sc_cols) <- names(progs)
for (i in seq_along(progs)) so[[names(progs)[i]]] <- so[[sc_cols[i]]]

## --------------------------------------------- 1. 细胞组成（比例漂移） ----
comp <- as.data.table(so@meta.data)[, .N, by = .(sample, group, celltype_short)]
comp[, frac := N / sum(N), by = sample]
comp_wide <- dcast(comp, celltype_short ~ sample, value.var = "frac", fill = 0)
fwrite(comp_wide, "03_results/singlecell/celltype_composition.csv")

comp_test <- comp[, {
  h <- frac[group == "Healthy"]; o <- frac[group == "OA"]
  tt <- if (length(h) > 1 && length(o) > 1) t.test(o, h) else list(p.value = NA, estimate = c(NA, NA))
  .(mean_HC = mean(h), mean_OA = mean(o), diff = mean(o) - mean(h), P = tt$p.value)
}, by = celltype_short][order(-abs(diff))]
comp_test[, FDR := p.adjust(P, "BH")]
fwrite(comp_test, "03_results/singlecell/celltype_composition_test.csv")
msg("--- 细胞组成 HC vs OA ---"); print(comp_test)

## -------------------------- 2. 亚型内程序评分：pseudobulk（样本为单位） ----
md <- as.data.table(so@meta.data, keep.rownames = "cell")
pb <- md[, lapply(.SD, mean), by = .(sample, group, celltype_short),
         .SDcols = names(progs)]
fwrite(pb, "03_results/singlecell/programme_pseudobulk.csv")

prog_test <- rbindlist(lapply(names(progs), function(p) {
  d <- pb[, .(sample, group, celltype_short, val = get(p))]
  d[, {
    h <- val[group == "Healthy"]; o <- val[group == "OA"]
    if (length(h) < 2 || length(o) < 2) return(.(mean_HC = mean(h), mean_OA = mean(o),
                                                 delta = NA_real_, P = NA_real_))
    tt <- t.test(o, h)
    sp <- sqrt(((length(o)-1)*var(o) + (length(h)-1)*var(h))/(length(o)+length(h)-2))
    .(mean_HC = mean(h), mean_OA = mean(o), delta = mean(o) - mean(h),
      hedges_g = (mean(o)-mean(h))/sp, P = tt$p.value)
  }, by = celltype_short][, programme := p][]
}))
prog_test[, FDR := p.adjust(P, "BH"), by = programme]
fwrite(prog_test, "03_results/singlecell/programme_by_celltype_test.csv")
msg("--- 剪接程序 亚型内 HC vs OA (pseudobulk) ---")
print(prog_test[programme == "Splicing"][order(delta)][
  , .(celltype_short, delta = round(delta,4), g = round(hedges_g,2),
      P = signif(P,3), FDR = signif(FDR,3))])
msg("--- 对照程序 meta ---")
print(prog_test[, .(mean_delta = round(mean(delta, na.rm=TRUE),4),
                    n_nominal = sum(P < 0.05, na.rm=TRUE)), by = programme])

## ----------------------------------------- 3. 组成漂移 vs 细胞内在变化 ----
# 若剪接程序下调仅由亚型比例改变造成，则"组成校正后"效应应消失。
# 用等权重（各亚型均权）重算全局评分，剔除比例影响。
glob_raw <- md[, .(val = mean(Splicing)), by = .(sample, group)]
glob_adj <- pb[, .(val = mean(Splicing)), by = .(sample, group)]   # 亚型等权
tt_raw <- t.test(glob_raw[group=="OA", val], glob_raw[group=="Healthy", val])
tt_adj <- t.test(glob_adj[group=="OA", val], glob_adj[group=="Healthy", val])
comp_adj <- data.table(
  measure = c("Raw (composition-confounded)", "Composition-adjusted (celltype-equal weight)"),
  mean_HC = c(mean(glob_raw[group=="Healthy", val]), mean(glob_adj[group=="Healthy", val])),
  mean_OA = c(mean(glob_raw[group=="OA", val]),      mean(glob_adj[group=="OA", val])),
  delta   = c(diff(rev(tt_raw$estimate)), diff(rev(tt_adj$estimate))),
  P       = c(tt_raw$p.value, tt_adj$p.value))
fwrite(comp_adj, "03_results/singlecell/composition_adjusted_splicing.csv")
msg("--- 组成校正前后 ---"); print(comp_adj)

## ----------------------------------------------------------- 4. SON ------
if ("SON" %in% rownames(so)) {
  md[, SON := FetchData(so, "SON")[, 1]]
  son_pb <- md[, .(SON = mean(SON)), by = .(sample, group, celltype_short)]
  son_test <- son_pb[, {
    h <- SON[group=="Healthy"]; o <- SON[group=="OA"]
    if (length(h) < 2 || length(o) < 2) .(delta = NA_real_, P = NA_real_)
    else .(mean_HC = mean(h), mean_OA = mean(o), delta = mean(o)-mean(h), P = t.test(o,h)$p.value)
  }, by = celltype_short][order(delta)]
  son_test[, FDR := p.adjust(P, "BH")]
  fwrite(son_test, "03_results/singlecell/SON_by_celltype.csv")
  msg("--- SON 亚型内 HC vs OA ---"); print(son_test)
  # SON 与剪接评分的单细胞相关
  r <- cor(md$SON, md$Splicing, method = "spearman")
  msg("单细胞水平 SON vs 剪接程序 Spearman rho = ", round(r, 3))
  writeLines(sprintf("single_cell_SON_vs_splicing_spearman_rho\t%.4f", r),
             "03_results/singlecell/SON_splicing_correlation.txt")
} else {
  msg("注意：SON 不在该过滤矩阵的 ", nrow(so), " 个基因中，跳过 SON 单细胞分析")
  writeLines("SON not present in GSE324993 filtered gene matrix",
             "03_results/singlecell/SON_absent_note.txt")
}

saveRDS(so, "03_results/intermediate/seurat_GSE324993.rds")
saveRDS(list(comp = comp_test, prog = prog_test, comp_adj = comp_adj),
        "03_results/intermediate/singlecell_results.rds")
msg("06 完成")
