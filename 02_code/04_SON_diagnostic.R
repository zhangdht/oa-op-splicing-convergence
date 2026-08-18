## =====================================================================
## bioinfo05 | Step 04 : Pre-specified hypothesis test for SON
##  (a) is SON differentially expressed?                    -> no (step 03)
##  (b) is the SON-anchored co-expression network rewired?
##  (c) is the spliceosome / nuclear-speckle programme moved?
## =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(limma)
})

ROOT  <- "D:/bioinfo05"
DIR_I <- file.path(ROOT, "03_results/intermediate")
DIR_S <- file.path(ROOT, "03_results/splice")
DIR_L <- file.path(ROOT, "07_logs")
dir.create(DIR_S, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(DIR_L, "04_SON_diagnostic.log"), open = "wt")
lg <- function(...) { m <- paste0(format(Sys.time(), "[%H:%M:%S] "), ...)
                      cat(m, "\n"); writeLines(m, log_con); flush(log_con) }

cohorts <- readRDS(file.path(DIR_I, "cohorts.rds"))
deg     <- readRDS(file.path(DIR_I, "deg_all.rds"))

## ---------------------------------------------------------------------
## (a) per-cohort SON statistics
## ---------------------------------------------------------------------
lg("=== (a) SON differential expression per cohort ===")
son_tab <- do.call(rbind, lapply(names(deg), function(nm) {
  d <- deg[[nm]]; s <- d[d$gene == "SON", ]
  if (!nrow(s)) return(NULL)
  data.frame(cohort = nm, disease = cohorts[[nm]]$disease,
             n = ncol(cohorts[[nm]]$expr),
             logFC = s$logFC, P = s$P.Value, FDR = s$adj.P.Val,
             rank = which(d$gene == "SON"), n_genes = nrow(d))
}))
print(son_tab); fwrite(son_tab, file.path(DIR_S, "SON_per_cohort_DE.csv"))

## ---------------------------------------------------------------------
## (b) differential co-expression: does SON's partner set change?
##     For every cohort, compute Spearman rho(SON, gene) separately in
##     cases and controls, then test the shift with Fisher's z.
## ---------------------------------------------------------------------
lg("=== (b) SON differential co-expression ===")
dcx <- list()
for (nm in names(cohorts)) {
  co <- cohorts[[nm]]; e <- co$expr; ph <- co$pheno
  if (!"SON" %in% rownames(e)) next
  nA <- sum(ph$group == "Control"); nB <- sum(ph$group == "Case")
  if (min(nA, nB) < 5) { lg("  ", nm, " skipped (n too small)"); next }
  ## keep reasonably expressed / variable genes (SON always retained)
  v <- apply(e, 1, sd, na.rm = TRUE)
  keep <- (is.finite(v) & v > quantile(v, .25, na.rm = TRUE)) | rownames(e) == "SON"
  e <- e[keep, , drop = FALSE]
  son <- e["SON", ]
  rC <- suppressWarnings(cor(t(e[, ph$group == "Control", drop = FALSE]),
                             son[ph$group == "Control"], method = "spearman"))[, 1]
  rD <- suppressWarnings(cor(t(e[, ph$group == "Case", drop = FALSE]),
                             son[ph$group == "Case"], method = "spearman"))[, 1]
  zC <- atanh(pmax(pmin(rC, .999), -.999)); zD <- atanh(pmax(pmin(rD, .999), -.999))
  se <- sqrt(1 / (nA - 3) + 1 / (nB - 3))
  zz <- (zD - zC) / se
  p  <- 2 * pnorm(-abs(zz))
  d  <- data.frame(gene = rownames(e), r_control = rC, r_case = rD,
                   delta_z = zz, P = p, FDR = p.adjust(p, "BH"),
                   cohort = nm, disease = co$disease, stringsAsFactors = FALSE)
  d <- d[d$gene != "SON", ]
  dcx[[nm]] <- d
  lg(sprintf("  %-10s  genes=%d  FDR<0.05 rewired partners = %d  (median |r| ctrl %.2f -> case %.2f)",
             nm, nrow(d), sum(d$FDR < 0.05, na.rm = TRUE),
             median(abs(d$r_control), na.rm = TRUE), median(abs(d$r_case), na.rm = TRUE)))
}
saveRDS(dcx, file.path(DIR_I, "SON_diffcoexp.rds"))

## ---------------------------------------------------------------------
## (c) spliceosome / nuclear-speckle programme
## ---------------------------------------------------------------------
lg("=== (c) splicing-machinery programme ===")
gmt <- function(f) {
  ln <- readLines(f); s <- strsplit(ln, "\t")
  setNames(lapply(s, function(x) unique(x[-c(1, 2)])), sapply(s, `[`, 1))
}
go  <- gmt(file.path(ROOT, "01_data/genesets/c5.go.bp.v7.5.1.symbols.gmt"))
kg  <- gmt(file.path(ROOT, "01_data/genesets/c2.cp.kegg.v7.5.1.symbols.gmt"))
rc  <- gmt(file.path(ROOT, "01_data/genesets/c2.cp.reactome.v7.5.1.symbols.gmt"))
all_sets <- c(go, kg, rc)

pick <- grep("SPLIC|SPECKLE|MRNA_PROCESSING|SPLICEOSOM", names(all_sets), value = TRUE)
lg("  candidate splicing gene sets: ", length(pick))
core <- unique(unlist(all_sets[c("KEGG_SPLICEOSOME",
                                 grep("^GOBP_MRNA_SPLICING_VIA_SPLICEOSOME$", names(all_sets), value = TRUE),
                                 grep("^GOBP_REGULATION_OF_MRNA_SPLICING_VIA_SPLICEOSOME$", names(all_sets), value = TRUE))]))
core <- core[!is.na(core)]
lg("  core spliceosome gene set size: ", length(core))
writeLines(core, file.path(DIR_S, "core_spliceosome_genes.txt"))

## simple single-sample score = mean z of the gene set
ss_score <- function(e, gs) {
  gs <- intersect(gs, rownames(e))
  z  <- t(scale(t(e[gs, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

sp_res <- do.call(rbind, lapply(names(cohorts), function(nm) {
  co <- cohorts[[nm]]; e <- co$expr; ph <- co$pheno
  sc <- ss_score(e, core)
  tt <- t.test(sc ~ ph$group)
  data.frame(cohort = nm, disease = co$disease,
             n_genes_used = length(intersect(core, rownames(e))),
             mean_control = mean(sc[ph$group == "Control"]),
             mean_case = mean(sc[ph$group == "Case"]),
             delta = mean(sc[ph$group == "Case"]) - mean(sc[ph$group == "Control"]),
             P = tt$p.value, stringsAsFactors = FALSE)
}))
sp_res$FDR <- p.adjust(sp_res$P, "BH")
print(sp_res); fwrite(sp_res, file.path(DIR_S, "spliceosome_programme_score.csv"))

## correlation between SON and the spliceosome programme
lg("=== SON vs spliceosome programme correlation ===")
cr <- do.call(rbind, lapply(names(cohorts), function(nm) {
  co <- cohorts[[nm]]; e <- co$expr
  sc <- ss_score(e, setdiff(core, "SON"))
  ct <- cor.test(e["SON", ], sc, method = "spearman")
  data.frame(cohort = nm, disease = co$disease,
             rho = unname(ct$estimate), P = ct$p.value)
}))
print(cr); fwrite(cr, file.path(DIR_S, "SON_vs_spliceosome_corr.csv"))

close(log_con)
