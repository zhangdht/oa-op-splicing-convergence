## =====================================================================
## bioinfo05 | Step 01 : Preprocessing of osteoporosis (OP) cohorts
## Project : SON as a shared splicing-regulatory hub in OP-OA comorbidity
## Author  : analysis pipeline
## =====================================================================
suppressPackageStartupMessages({
  library(limma); library(data.table)
  library(hgu133plus2.db); library(hgu133a.db); library(AnnotationDbi)
})

ROOT <- "D:/bioinfo05"
DIR_GEO  <- file.path(ROOT, "01_data/geo")
DIR_OUT  <- file.path(ROOT, "01_data/processed")
DIR_LOG  <- file.path(ROOT, "07_logs")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(DIR_LOG, "01_preprocess_OP.log"), open = "wt")
lg <- function(...) { msg <- paste0(format(Sys.time(), "[%H:%M:%S] "), ...)
                      cat(msg, "\n"); writeLines(msg, log_con); flush(log_con) }

## ---------------------------------------------------------------------
## helper : read a GEO series matrix file into expression + phenotype
## ---------------------------------------------------------------------
read_series_matrix <- function(path) {
  lines <- readLines(gzfile(path))
  beg <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end",   lines)
  tab <- fread(text = paste(lines[(beg + 1):(end - 1)], collapse = "\n"),
               sep = "\t", header = TRUE, data.table = FALSE)
  rownames(tab) <- tab[[1]]; tab <- tab[, -1, drop = FALSE]
  expr <- as.matrix(tab)
  storage.mode(expr) <- "numeric"

  gsm  <- strsplit(gsub('"', '', lines[grep("^!Sample_geo_accession", lines)[1]]), "\t")[[1]][-1]
  chr  <- lines[grep("^!Sample_characteristics_ch1", lines)]
  chl  <- lapply(chr, function(x) strsplit(gsub('"', '', x), "\t")[[1]][-1])
  pdat <- data.frame(gsm = gsm, stringsAsFactors = FALSE)
  for (i in seq_along(chl)) pdat[[paste0("ch", i)]] <- chl[[i]]
  ttl  <- lines[grep("^!Sample_title", lines)][1]
  if (!is.na(ttl)) pdat$title <- strsplit(gsub('"', '', ttl), "\t")[[1]][-1]
  list(expr = expr, pheno = pdat)
}

## ---------------------------------------------------------------------
## helper : probe -> gene symbol collapse (keep max mean expression)
## ---------------------------------------------------------------------
collapse_to_symbol <- function(expr, db) {
  sym <- AnnotationDbi::mapIds(db, keys = rownames(expr), column = "SYMBOL",
                               keytype = "PROBEID", multiVals = "first")
  keep <- !is.na(sym)
  expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  ord  <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
  expr <- expr[ord, , drop = FALSE]; sym <- sym[ord]
  expr <- expr[!duplicated(sym), , drop = FALSE]
  rownames(expr) <- sym[!duplicated(sym)]
  expr
}

## ---------------------------------------------------------------------
## helper : ensure log2 scale
## ---------------------------------------------------------------------
ensure_log2 <- function(expr) {
  qx <- as.numeric(quantile(expr, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
  needs_log <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (needs_log) { expr[expr < 0] <- 0; expr <- log2(expr + 1)
                   lg("    log2 transform applied") }
  expr
}

## =====================================================================
## GSE56815  |  GPL96  |  monocytes, 40 high-BMD vs 40 low-BMD
##            -> primary OP discovery cohort (best powered)
## =====================================================================
lg("== GSE56815 ==")
s <- read_series_matrix(file.path(DIR_GEO, "GSE56815_series_matrix.txt.gz"))
lg("  raw probes x samples: ", paste(dim(s$expr), collapse = " x "))

bmd_col <- which(sapply(s$pheno, function(x) any(grepl("bone mineral density", x))))
men_col <- which(sapply(s$pheno, function(x) any(grepl("^state:", x))))
grp <- ifelse(grepl("low BMD", s$pheno[[bmd_col]]), "OP", "Control")
ph  <- data.frame(sample = s$pheno$gsm,
                  group  = factor(grp, levels = c("Control", "OP")),
                  menopause = sub("state: ", "", s$pheno[[men_col]]),
                  series = "GSE56815", tissue = "monocyte",
                  stringsAsFactors = FALSE)
lg("  groups: ", paste(names(table(ph$group)), table(ph$group), collapse = " / "))

e <- ensure_log2(s$expr)
e <- collapse_to_symbol(e, hgu133a.db)
e <- normalizeBetweenArrays(e, method = "quantile")
lg("  final genes x samples: ", paste(dim(e), collapse = " x "))
stopifnot(identical(colnames(e), ph$sample))
fwrite(data.frame(symbol = rownames(e), e), file.path(DIR_OUT, "GSE56815_expr_symbol.tsv.gz"), sep = "\t")
fwrite(ph, file.path(DIR_OUT, "GSE56815_pheno.tsv"), sep = "\t")

## =====================================================================
## GSE7158   |  GPL570 |  monocytes, high vs low peak bone mass
##            -> OP validation cohort
## =====================================================================
lg("== GSE7158 ==")
s <- read_series_matrix(file.path(DIR_GEO, "GSE7158_series_matrix.txt.gz"))
lg("  raw probes x samples: ", paste(dim(s$expr), collapse = " x "))
lab <- s$pheno$ch1
grp <- ifelse(grepl("Low", lab, ignore.case = TRUE), "OP", "Control")
ph  <- data.frame(sample = s$pheno$gsm,
                  group  = factor(grp, levels = c("Control", "OP")),
                  series = "GSE7158", tissue = "monocyte", stringsAsFactors = FALSE)
lg("  groups: ", paste(names(table(ph$group)), table(ph$group), collapse = " / "))

e <- ensure_log2(s$expr)
e <- collapse_to_symbol(e, hgu133plus2.db)
e <- normalizeBetweenArrays(e, method = "quantile")
lg("  final genes x samples: ", paste(dim(e), collapse = " x "))
fwrite(data.frame(symbol = rownames(e), e), file.path(DIR_OUT, "GSE7158_expr_symbol.tsv.gz"), sep = "\t")
fwrite(ph, file.path(DIR_OUT, "GSE7158_pheno.tsv"), sep = "\t")

## =====================================================================
## GSE35958  |  GPL570 |  bone-marrow MSC, primary OP vs non-OP elderly
##            -> OP bone-compartment cohort
## =====================================================================
lg("== GSE35958 ==")
s <- read_series_matrix(file.path(DIR_GEO, "GSE35958_series_matrix.txt.gz"))
lg("  raw probes x samples: ", paste(dim(s$expr), collapse = " x "))
is_op <- apply(s$pheno[, grep("^ch", names(s$pheno))], 1,
               function(r) any(grepl("primary osteoporosis", r)))
ph <- data.frame(sample = s$pheno$gsm,
                 group  = factor(ifelse(is_op, "OP", "Control"), levels = c("Control", "OP")),
                 series = "GSE35958", tissue = "BM-MSC", stringsAsFactors = FALSE)
lg("  groups: ", paste(names(table(ph$group)), table(ph$group), collapse = " / "))

e <- ensure_log2(s$expr)
e <- collapse_to_symbol(e, hgu133plus2.db)
e <- normalizeBetweenArrays(e, method = "quantile")
lg("  final genes x samples: ", paste(dim(e), collapse = " x "))
fwrite(data.frame(symbol = rownames(e), e), file.path(DIR_OUT, "GSE35958_expr_symbol.tsv.gz"), sep = "\t")
fwrite(ph, file.path(DIR_OUT, "GSE35958_pheno.tsv"), sep = "\t")

lg("OP preprocessing finished.")
close(log_con)
