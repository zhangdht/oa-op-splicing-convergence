## =====================================================================
## bioinfo05 | Step 03 : Differential expression in every cohort
##            + cross-cohort meta-analysis within each disease
## =====================================================================
suppressPackageStartupMessages({
  library(limma); library(data.table)
})

ROOT  <- "D:/bioinfo05"
DIR_R <- file.path(ROOT, "03_results/deg")
DIR_I <- file.path(ROOT, "03_results/intermediate")
DIR_L <- file.path(ROOT, "07_logs")
dir.create(DIR_R, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(DIR_L, "03_DEG.log"), open = "wt")
lg <- function(...) { m <- paste0(format(Sys.time(), "[%H:%M:%S] "), ...)
                      cat(m, "\n"); writeLines(m, log_con); flush(log_con) }

cohorts <- readRDS(file.path(DIR_I, "cohorts.rds"))

## ---------------------------------------------------------------------
## limma differential expression; uses a within-donor blocked design when
## the cohort provides paired samples (duplicateCorrelation).
## ---------------------------------------------------------------------
run_limma <- function(co, name) {
  e  <- co$expr
  ph <- co$pheno
  ## drop genes with near-zero variance
  v  <- apply(e, 1, var, na.rm = TRUE)
  e  <- e[is.finite(v) & v > 1e-8, , drop = FALSE]

  design <- model.matrix(~ group, data = ph)
  paired <- !all(is.na(ph$donor)) && anyDuplicated(ph$donor[!is.na(ph$donor)]) > 0

  if (paired) {
    dc  <- duplicateCorrelation(e, design, block = ph$donor)
    fit <- lmFit(e, design, block = ph$donor, correlation = dc$consensus.correlation)
    lg("    paired design, consensus correlation = ", round(dc$consensus.correlation, 3))
  } else {
    fit <- lmFit(e, design)
  }
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  tt  <- topTable(fit, coef = "groupCase", number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  tt$cohort <- name; tt$disease <- co$disease
  tt
}

deg <- list()
for (nm in names(cohorts)) {
  lg("== DEG ", nm, " (", cohorts[[nm]]$disease, ") ==")
  tt <- run_limma(cohorts[[nm]], nm)
  deg[[nm]] <- tt
  n_sig <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  n_sig_fc <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
  lg("    genes tested: ", nrow(tt),
     " | FDR<0.05: ", n_sig, " | FDR<0.05 & |logFC|>0.5: ", n_sig_fc)
  son <- tt[tt$gene == "SON", ]
  if (nrow(son)) lg(sprintf("    SON: logFC=%.3f  P=%.3g  FDR=%.3g",
                            son$logFC, son$P.Value, son$adj.P.Val))
  fwrite(tt, file.path(DIR_R, paste0("DEG_", nm, ".csv")))
}
saveRDS(deg, file.path(DIR_I, "deg_all.rds"))

## =====================================================================
## Cross-cohort meta-analysis (random-effects on the limma t-derived
## standardised effect; Stouffer weighted-Z as the primary summary)
## =====================================================================
meta_z <- function(deg_list, cohort_names, cohorts) {
  ## weighted Z with weights = sqrt(n)
  common <- Reduce(intersect, lapply(cohort_names, function(x) deg_list[[x]]$gene))
  lg("    common genes: ", length(common))
  Z <- W <- FC <- matrix(NA_real_, nrow = length(common), ncol = length(cohort_names),
                         dimnames = list(common, cohort_names))
  for (cn in cohort_names) {
    d <- deg_list[[cn]]; rownames(d) <- d$gene; d <- d[common, ]
    z <- qnorm(pmax(pmin(1 - d$P.Value / 2, 1 - 1e-16), 1e-16)) * sign(d$logFC)
    Z[, cn]  <- z
    W[, cn]  <- sqrt(ncol(cohorts[[cn]]$expr))
    FC[, cn] <- d$logFC
  }
  zc <- rowSums(Z * W, na.rm = TRUE) / sqrt(rowSums(W^2, na.rm = TRUE))
  p  <- 2 * pnorm(-abs(zc))
  data.frame(gene = common, meta_Z = zc, meta_P = p,
             meta_FDR = p.adjust(p, "BH"),
             mean_logFC = rowMeans(FC, na.rm = TRUE),
             n_consistent = apply(sign(FC), 1, function(x) max(table(x[!is.na(x)]))),
             stringsAsFactors = FALSE)
}

lg("== meta-analysis : OA ==")
oa_names <- names(cohorts)[sapply(cohorts, function(x) x$disease) == "OA"]
meta_oa  <- meta_z(deg, oa_names, cohorts)
meta_oa  <- meta_oa[order(meta_oa$meta_P), ]
fwrite(meta_oa, file.path(DIR_R, "META_OA.csv"))
lg("    FDR<0.05: ", sum(meta_oa$meta_FDR < 0.05))
lg("    SON: ", paste(capture.output(print(meta_oa[meta_oa$gene == "SON", ])), collapse = " | "))

lg("== meta-analysis : OP ==")
op_names <- names(cohorts)[sapply(cohorts, function(x) x$disease) == "OP"]
meta_op  <- meta_z(deg, op_names, cohorts)
meta_op  <- meta_op[order(meta_op$meta_P), ]
fwrite(meta_op, file.path(DIR_R, "META_OP.csv"))
lg("    FDR<0.05: ", sum(meta_op$meta_FDR < 0.05))
lg("    SON: ", paste(capture.output(print(meta_op[meta_op$gene == "SON", ])), collapse = " | "))

saveRDS(list(OA = meta_oa, OP = meta_op), file.path(DIR_I, "meta_deg.rds"))

## =====================================================================
## Shared differential signal (the comorbidity candidate pool)
## =====================================================================
lg("== shared OP/OA differential genes ==")
m <- merge(meta_oa[, c("gene", "meta_Z", "meta_P", "meta_FDR", "mean_logFC")],
           meta_op[, c("gene", "meta_Z", "meta_P", "meta_FDR", "mean_logFC")],
           by = "gene", suffixes = c("_OA", "_OP"))
m$same_direction <- sign(m$meta_Z_OA) == sign(m$meta_Z_OP)
m$shared_score   <- -log10(m$meta_P_OA) - log10(m$meta_P_OP)
m <- m[order(-m$shared_score), ]
fwrite(m, file.path(DIR_R, "shared_OP_OA_meta.csv"))

sig <- subset(m, meta_FDR_OA < 0.05 & meta_FDR_OP < 0.05)
lg("    genes FDR<0.05 in BOTH diseases: ", nrow(sig))
lg("    of which concordant direction : ", sum(sig$same_direction))
fwrite(sig, file.path(DIR_R, "shared_OP_OA_FDR05.csv"))
lg("    SON rank in shared score: ", which(m$gene == "SON"), " / ", nrow(m))
print(head(m[, c("gene", "meta_Z_OA", "meta_FDR_OA", "meta_Z_OP", "meta_FDR_OP",
                 "same_direction", "shared_score")], 25))
cat("\n--- SON ---\n"); print(m[m$gene == "SON", ])

close(log_con)
