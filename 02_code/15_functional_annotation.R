###############################################################################
# 15_functional_annotation.R
# 目的：为两类基因集做正式的功能注释，补齐论文 Results 中缺失的机制解释
#   (A) WGCNA 跨病种 21 个共享基因 -> 检验它们是否真为髓系/炎症来源
#   (B) 间充质轴剪接程序核心基因 -> 确认剪接注释并给出可解释的 GO 结构
# 输出：03_results/functional/
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
})
ROOT <- "D:/bioinfo05"; setwd(ROOT)
dir.create("03_results/functional", showWarnings = FALSE, recursive = TRUE)

set.seed(123)

# ---------------------------------------------------------------------------
# 背景基因集：用 7 队列共同检测到的基因作为 universe，避免用全基因组造成的
# 富集 P 值虚高（这是 GO 分析最常见的方法学错误）
# ---------------------------------------------------------------------------
cohorts  <- readRDS("03_results/intermediate/cohorts.rds")
universe <- Reduce(intersect, lapply(cohorts, function(x) rownames(x$expr)))
cat("Universe genes (detected in all 7 cohorts):", length(universe), "\n")

sym2eg <- function(x) {
  suppressWarnings(suppressMessages(
    bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ))
}
uni_eg <- sym2eg(universe)

run_go <- function(genes, label, ont = "BP") {
  g <- intersect(genes, universe)
  if (length(g) < 3) {
    cat(sprintf("[%s] too few genes mapped (%d), skipped\n", label, length(g)))
    return(NULL)
  }
  eg <- sym2eg(g)
  res <- tryCatch(
    enrichGO(gene = eg$ENTREZID, universe = uni_eg$ENTREZID, OrgDb = org.Hs.eg.db,
             ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05,
             qvalueCutoff = 0.25, readable = TRUE),
    error = function(e) { cat("[", label, "] enrichGO error:", conditionMessage(e), "\n"); NULL })
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    cat(sprintf("[%s] n=%d mapped=%d -> no GO term passed FDR\n", label, length(g), nrow(eg)))
    return(data.table(set = label, n_input = length(g), n_mapped = nrow(eg),
                      ID = NA_character_, Description = "no significant term",
                      GeneRatio = NA_character_, pvalue = NA_real_,
                      p.adjust = NA_real_, geneID = NA_character_))
  }
  dt <- as.data.table(as.data.frame(res))
  dt[, `:=`(set = label, n_input = length(g), n_mapped = nrow(eg))]
  cat(sprintf("[%s] n=%d mapped=%d -> %d GO BP terms at FDR<0.05\n",
              label, length(g), nrow(eg), nrow(dt)))
  dt[order(p.adjust)][1:min(.N, 25)]
}

out <- list()

# ---------------------------------------------------------------------------
# (A) WGCNA 跨病种共享基因
# ---------------------------------------------------------------------------
shared_file <- "03_results/wgcna/shared_module_genes.txt"
if (file.exists(shared_file)) {
  shared <- unique(trimws(readLines(shared_file)))
  shared <- shared[nchar(shared) > 0 & !grepl("^#", shared)]
  cat("\n== (A) WGCNA cross-disease shared genes:", length(shared), "==\n")
  cat(paste(shared, collapse = ", "), "\n")
  out[["WGCNA_shared"]] <- run_go(shared, "WGCNA cross-disease shared genes")
}

# ---------------------------------------------------------------------------
# (B) 剪接程序核心基因
# ---------------------------------------------------------------------------
sp_file <- "03_results/splice/core_spliceosome_genes.txt"
if (file.exists(sp_file)) {
  sp <- unique(trimws(readLines(sp_file)))
  sp <- sp[nchar(sp) > 0]
  cat("\n== (B) Core splicing programme genes:", length(sp), "==\n")
  out[["splicing"]] <- run_go(sp, "Core splicing programme")
}

# ---------------------------------------------------------------------------
# (C) 间充质轴一致下调通路的前导基因（leading-edge 近似）
#     取 A 轴剪接通路成员 ∩ 两病种 meta z 均为负的基因
# ---------------------------------------------------------------------------
# 关键：间充质轴 = OA 软骨 meta vs OP BM-MSC 单队列（GSE35958），
# 不能用 OP 三队列合并 meta，否则混入造血分区（单核）信号
meta_oa_f <- "03_results/deg/META_OA.csv"
msc_f     <- "03_results/deg/DEG_GSE35958.csv"
mono_f    <- "03_results/deg/DEG_GSE56815.csv"
if (all(file.exists(meta_oa_f, msc_f)) && file.exists(sp_file)) {
  moa  <- fread(meta_oa_f)
  msc  <- fread(msc_f)
  mm <- merge(moa[, .(gene, z_oa = meta_Z)],
              msc[, .(gene, t_msc = t, p_msc = P.Value)], by = "gene")
  if (file.exists(mono_f)) {
    mono <- fread(mono_f)
    mm <- merge(mm, mono[, .(gene, t_mono = t)], by = "gene", all.x = TRUE)
  }
  lead <- mm[gene %in% sp & z_oa < 0 & t_msc < 0][order(z_oa + t_msc)]
  fwrite(lead, "03_results/functional/mesenchymal_leading_edge_genes.csv")
  cat("\n== (C) Splicing genes down in BOTH OA cartilage and OP BM-MSC:", nrow(lead), "==\n")
  if (nrow(lead) > 0) cat(paste(head(lead$gene, 30), collapse = ", "), "\n")

  # 分区特异性检验：这些基因在造血分区（单核）中是否也下调？
  if ("t_mono" %in% names(lead) && nrow(lead) > 3) {
    n_down_mono <- sum(lead$t_mono < 0, na.rm = TRUE)
    n_tot <- sum(!is.na(lead$t_mono))
    bt <- binom.test(n_down_mono, n_tot, p = 0.5)
    cat(sprintf("   of these, %d/%d (%.0f%%) are also down in OP monocytes (binom P = %.3g)\n",
                n_down_mono, n_tot, 100 * n_down_mono / n_tot, bt$p.value))
    writeLines(c(
      sprintf("Splicing genes concordantly down in OA cartilage and OP BM-MSC: %d", nrow(lead)),
      sprintf("Of these, also down in OP monocytes (GSE56815): %d/%d (%.1f%%)",
              n_down_mono, n_tot, 100 * n_down_mono / n_tot),
      sprintf("Binomial test vs 50%%: P = %.4g", bt$p.value),
      "Interpretation: a proportion at or below 50% indicates the mesenchymal splicing",
      "down-regulation is not mirrored in the haematopoietic compartment."
    ), "03_results/functional/compartment_specificity_of_leading_edge.txt")
  }
  if (nrow(lead) >= 5) out[["leading_edge"]] <- run_go(lead$gene, "Mesenchymal splicing leading edge")
}

res_all <- rbindlist(Filter(Negate(is.null), out), fill = TRUE)
if (nrow(res_all) > 0) {
  keep <- intersect(c("set", "n_input", "n_mapped", "ID", "Description",
                      "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count", "geneID"),
                    names(res_all))
  fwrite(res_all[, ..keep], "03_results/functional/GO_enrichment_all.csv")
  cat("\nWritten: 03_results/functional/GO_enrichment_all.csv  (", nrow(res_all), "rows )\n")
}

sink("03_results/functional/summary.txt")
cat("Functional annotation summary\n")
cat("Universe:", length(universe), "genes detected in all seven cohorts\n\n")
if (nrow(res_all) > 0) {
  for (s in unique(res_all$set)) {
    sub <- res_all[set == s]
    cat("###", s, " (input", sub$n_input[1], "genes, mapped", sub$n_mapped[1], ")\n")
    if (all(is.na(sub$ID))) {
      cat("  no GO BP term passed FDR < 0.05\n\n")
    } else {
      for (i in seq_len(min(nrow(sub), 10)))
        cat(sprintf("  %-52s FDR=%.3g  n=%s\n", substr(sub$Description[i], 1, 52),
                    sub$p.adjust[i], sub$GeneRatio[i]))
      cat("\n")
    }
  }
}
sink()
cat("\nDONE\n")
