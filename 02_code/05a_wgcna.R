## =====================================================================
## bioinfo05 | Step 05 : Unbiased discovery of the shared OP-OA programme
##   5.1  WGCNA in the two best-powered discovery cohorts
##   5.2  module-trait association
##   5.3  cross-disease module correspondence
## =====================================================================
suppressPackageStartupMessages({
  library(WGCNA); library(data.table); library(limma)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)

ROOT  <- "D:/bioinfo05"
DIR_I <- file.path(ROOT, "03_results/intermediate")
DIR_W <- file.path(ROOT, "03_results/wgcna")
DIR_E <- file.path(ROOT, "03_results/enrich")
DIR_L <- file.path(ROOT, "07_logs")
for (d in c(DIR_W, DIR_E)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(DIR_L, "05_unbiased_discovery.log"), open = "wt")
lg <- function(...) { m <- paste0(format(Sys.time(), "[%H:%M:%S] "), ...)
                      cat(m, "\n"); writeLines(m, log_con); flush(log_con) }

cohorts <- readRDS(file.path(DIR_I, "cohorts.rds"))

## =====================================================================
## 5.1 / 5.2  WGCNA
## =====================================================================
run_wgcna <- function(nm, n_top = 5000) {
  lg("=== WGCNA : ", nm, " ===")
  co <- cohorts[[nm]]; e <- co$expr; ph <- co$pheno
  mad_v <- apply(e, 1, mad, na.rm = TRUE)
  e <- e[order(mad_v, decreasing = TRUE)[seq_len(min(n_top, nrow(e)))], , drop = FALSE]
  datExpr <- t(e)
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  lg("  matrix: ", paste(dim(datExpr), collapse = " x "))

  powers <- c(1:10, seq(12, 30, 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",
                           verbose = 0)
  pw <- sft$powerEstimate
  if (is.na(pw)) {
    ## fall back to the smallest power reaching R^2 >= 0.8, else 12
    idx <- which(-sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq >= 0.80)
    pw  <- if (length(idx)) sft$fitIndices$Power[idx[1]] else 12
  }
  lg("  soft-threshold power = ", pw)

  net <- blockwiseModules(datExpr, power = pw, networkType = "signed",
                          TOMType = "signed", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          maxBlockSize = 6000, verbose = 0, randomSeed = 42)
  mc <- labels2colors(net$colors)
  lg("  modules: ", length(unique(mc)), " | sizes: ",
     paste(head(sort(table(mc), decreasing = TRUE), 12), collapse = ","))

  MEs <- moduleEigengenes(datExpr, mc)$eigengenes
  MEs <- orderMEs(MEs)
  trait <- as.numeric(ph$group[match(rownames(datExpr), ph$sample)] == "Case")
  mtc <- cor(MEs, trait, use = "p")
  mtp <- corPvalueStudent(mtc, nrow(datExpr))
  res <- data.frame(module = sub("^ME", "", rownames(mtc)),
                    cor = mtc[, 1], p = mtp[, 1],
                    size = as.integer(table(mc)[sub("^ME", "", rownames(mtc))]))
  res <- res[order(res$p), ]
  lg("  top module-trait associations:")
  for (i in seq_len(min(5, nrow(res))))
    lg(sprintf("    %-14s n=%4d  r=%+.3f  P=%.3g",
               res$module[i], res$size[i], res$cor[i], res$p[i]))

  ## gene significance & module membership
  gs <- cor(datExpr, trait, use = "p")[, 1]
  list(cohort = nm, power = pw, colors = setNames(mc, colnames(datExpr)),
       MEs = MEs, moduleTrait = res, geneSignificance = gs,
       genes = colnames(datExpr), trait = trait, datExpr = datExpr)
}

w_oa <- run_wgcna("GSE114007")   # OA discovery, cartilage RNA-seq
w_op <- run_wgcna("GSE56815")    # OP discovery, monocytes n=80
saveRDS(list(OA = w_oa, OP = w_op), file.path(DIR_I, "wgcna.rds"))

fwrite(w_oa$moduleTrait, file.path(DIR_W, "moduleTrait_OA_GSE114007.csv"))
fwrite(w_op$moduleTrait, file.path(DIR_W, "moduleTrait_OP_GSE56815.csv"))
fwrite(data.frame(gene = names(w_oa$colors), module = w_oa$colors,
                  GS = w_oa$geneSignificance[names(w_oa$colors)]),
       file.path(DIR_W, "genes_modules_OA.csv"))
fwrite(data.frame(gene = names(w_op$colors), module = w_op$colors,
                  GS = w_op$geneSignificance[names(w_op$colors)]),
       file.path(DIR_W, "genes_modules_OP.csv"))

## =====================================================================
## 5.3  cross-disease module correspondence
## =====================================================================
lg("=== cross-disease module overlap ===")
sig_oa <- subset(w_oa$moduleTrait, p < 0.05 & module != "grey")
sig_op <- subset(w_op$moduleTrait, p < 0.05 & module != "grey")
lg("  OA trait-associated modules: ", paste(sig_oa$module, collapse = ", "))
lg("  OP trait-associated modules: ", paste(sig_op$module, collapse = ", "))

univ <- intersect(names(w_oa$colors), names(w_op$colors))
lg("  shared gene universe: ", length(univ))
ov <- expand.grid(OA = sig_oa$module, OP = sig_op$module, stringsAsFactors = FALSE)
ov$n_OA <- ov$n_OP <- ov$n_shared <- NA_integer_; ov$P <- NA_real_
for (i in seq_len(nrow(ov))) {
  a <- intersect(names(w_oa$colors)[w_oa$colors == ov$OA[i]], univ)
  b <- intersect(names(w_op$colors)[w_op$colors == ov$OP[i]], univ)
  s <- intersect(a, b)
  ov$n_OA[i] <- length(a); ov$n_OP[i] <- length(b); ov$n_shared[i] <- length(s)
  ov$P[i] <- phyper(length(s) - 1, length(a), length(univ) - length(a), length(b),
                    lower.tail = FALSE)
}
ov$FDR <- p.adjust(ov$P, "BH")
ov$jaccard <- ov$n_shared / (ov$n_OA + ov$n_OP - ov$n_shared)
ov <- ov[order(ov$P), ]
print(head(ov, 12)); fwrite(ov, file.path(DIR_W, "cross_disease_module_overlap.csv"))

## the consensus comorbidity gene set = genes in trait-associated modules of BOTH
oa_genes <- intersect(names(w_oa$colors)[w_oa$colors %in% sig_oa$module], univ)
op_genes <- intersect(names(w_op$colors)[w_op$colors %in% sig_op$module], univ)
shared_module_genes <- intersect(oa_genes, op_genes)
lg("  OA trait-module genes: ", length(oa_genes),
   " | OP trait-module genes: ", length(op_genes),
   " | intersection: ", length(shared_module_genes))
writeLines(shared_module_genes, file.path(DIR_W, "shared_module_genes.txt"))

## =====================================================================

lg("WGCNA + module correspondence finished.")
close(log_con)
