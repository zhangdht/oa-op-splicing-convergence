## =====================================================================
## bioinfo05 | Step 02 : Harmonisation of OA cohorts + global QC
## =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(ggplot2)
})

ROOT <- "D:/bioinfo05"
DIR_P <- file.path(ROOT, "01_data/processed")
DIR_R <- file.path(ROOT, "03_results/intermediate")
DIR_F <- file.path(ROOT, "04_figures")
DIR_L <- file.path(ROOT, "07_logs")
dir.create(DIR_R, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(DIR_L, "02_harmonise_QC.log"), open = "wt")
lg <- function(...) { m <- paste0(format(Sys.time(), "[%H:%M:%S] "), ...)
                      cat(m, "\n"); writeLines(m, log_con); flush(log_con) }

## ---------------------------------------------------------------------
## Excel date-corrupted gene symbols  ->  official HGNC symbols
## (a well documented artefact of spreadsheet handling of gene lists)
## ---------------------------------------------------------------------
fix_date_symbols <- function(x) {
  d <- grepl("^(19|20)[0-9]{2}-[0-9]{2}-[0-9]{2}", x)
  if (!any(d)) return(x)
  mo <- c("01" = "JAN", "02" = "FEB", "03" = "MAR", "04" = "APR",
          "05" = "MAY", "06" = "JUN", "07" = "JUL", "08" = "AUG",
          "09" = "SEP", "10" = "OCT", "11" = "NOV", "12" = "DEC")
  fam <- c(JAN = "JADE", FEB = "FEB", MAR = "MARCHF", APR = "APR",
           MAY = "MAY",  JUN = "JUN", JUL = "JUL",    AUG = "AUG",
           SEP = "SEPTIN", OCT = "OCT", NOV = "NOV",  DEC = "DEC")
  for (i in which(d)) {
    p <- strsplit(x[i], "-")[[1]]
    mm <- mo[p[2]]; dd <- as.integer(p[3])
    if (!is.na(mm) && mm %in% names(fam)) x[i] <- paste0(fam[[mm]], dd)
  }
  x
}

read_expr <- function(gse) {
  e <- fread(file.path(DIR_P, paste0(gse, "_expr_symbol.tsv.gz")), data.table = FALSE)
  sym <- fix_date_symbols(as.character(e[[1]]))
  e <- as.matrix(e[, -1, drop = FALSE]); storage.mode(e) <- "numeric"
  ## collapse duplicated symbols by maximum mean expression
  ord <- order(rowMeans(e, na.rm = TRUE), decreasing = TRUE)
  e <- e[ord, , drop = FALSE]; sym <- sym[ord]
  keep <- !duplicated(sym) & !is.na(sym) & sym != ""
  e <- e[keep, , drop = FALSE]; rownames(e) <- sym[keep]
  e
}

## =====================================================================
## Load every cohort, attach a harmonised two-level phenotype
## =====================================================================
cohorts <- list()

## ---- OA discovery : GSE114007 (knee cartilage, RNA-seq, log2) --------
lg("== GSE114007 ==")
e  <- read_expr("GSE114007")
ph <- fread(file.path(DIR_P, "GSE114007_pheno.tsv"), data.table = FALSE)
ph$sample <- ph$sample; rownames(ph) <- ph$sample
ph <- ph[colnames(e), ]
cohorts$GSE114007 <- list(
  expr = e, disease = "OA", tissue = "knee cartilage", platform = "RNA-seq",
  pheno = data.frame(sample = ph$sample,
                     group = factor(ifelse(ph$group == "OA", "Case", "Control"),
                                    levels = c("Control", "Case")),
                     age = suppressWarnings(as.numeric(ph$age)), sex = ph$sex,
                     donor = NA_character_, stringsAsFactors = FALSE))
lg("  ", paste(dim(e), collapse = " x "), " | ",
   paste(names(table(cohorts$GSE114007$pheno$group)),
         table(cohorts$GSE114007$pheno$group), collapse = " / "))

## ---- OA validation : GSE57218 (RAAK, paired preserved / OA) ----------
lg("== GSE57218 ==")
e  <- read_expr("GSE57218")
ph <- fread(file.path(DIR_P, "GSE57218_pheno.tsv"), data.table = FALSE)
rownames(ph) <- ph$sample; ph <- ph[colnames(e), ]
keep <- ph$group %in% c("OA", "Preserved")          # within-donor contrast
e <- e[, keep, drop = FALSE]; ph <- ph[keep, ]
cohorts$GSE57218 <- list(
  expr = e, disease = "OA", tissue = "knee/hip cartilage", platform = "GPL570",
  pheno = data.frame(sample = ph$sample,
                     group = factor(ifelse(ph$group == "OA", "Case", "Control"),
                                    levels = c("Control", "Case")),
                     age = suppressWarnings(as.numeric(ph$age)), sex = ph$sex,
                     donor = as.character(ph$donor_id), stringsAsFactors = FALSE))
lg("  ", paste(dim(e), collapse = " x "), " | ",
   paste(names(table(cohorts$GSE57218$pheno$group)),
         table(cohorts$GSE57218$pheno$group), collapse = " / "))

## ---- OA validation : GSE117999 (paired, knee cartilage) --------------
lg("== GSE117999 ==")
e  <- read_expr("GSE117999")
ph <- fread(file.path(DIR_P, "GSE117999_pheno.tsv"), data.table = FALSE)
rownames(ph) <- ph$sample; ph <- ph[colnames(e), ]
cohorts$GSE117999 <- list(
  expr = e, disease = "OA", tissue = "knee cartilage", platform = "GPL10295",
  pheno = data.frame(sample = ph$sample,
                     group = factor(ifelse(ph$group == "OA", "Case", "Control"),
                                    levels = c("Control", "Case")),
                     age = suppressWarnings(as.numeric(ph$age)), sex = ph$sex,
                     donor = as.character(ph$donor_id), stringsAsFactors = FALSE))
lg("  ", paste(dim(e), collapse = " x "), " | ",
   paste(names(table(cohorts$GSE117999$pheno$group)),
         table(cohorts$GSE117999$pheno$group), collapse = " / "))

## ---- OA validation : GSE169077 --------------------------------------
lg("== GSE169077 ==")
e  <- read_expr("GSE169077")
ph <- fread(file.path(DIR_P, "GSE169077_pheno.tsv"), data.table = FALSE)
rownames(ph) <- ph$sample; ph <- ph[colnames(e), ]
cohorts$GSE169077 <- list(
  expr = e, disease = "OA", tissue = "knee cartilage", platform = "GPL96",
  pheno = data.frame(sample = ph$sample,
                     group = factor(ifelse(ph$group == "OA", "Case", "Control"),
                                    levels = c("Control", "Case")),
                     age = NA_real_, sex = NA_character_, donor = NA_character_,
                     stringsAsFactors = FALSE))
lg("  ", paste(dim(e), collapse = " x "), " | ",
   paste(names(table(cohorts$GSE169077$pheno$group)),
         table(cohorts$GSE169077$pheno$group), collapse = " / "))

## ---- OP cohorts ------------------------------------------------------
for (g in c("GSE56815", "GSE7158", "GSE35958")) {
  lg("== ", g, " ==")
  e  <- read_expr(g)
  ph <- fread(file.path(DIR_P, paste0(g, "_pheno.tsv")), data.table = FALSE)
  rownames(ph) <- ph$sample; ph <- ph[colnames(e), ]
  cohorts[[g]] <- list(
    expr = e, disease = "OP", tissue = ph$tissue[1],
    platform = ifelse(g == "GSE56815", "GPL96", "GPL570"),
    pheno = data.frame(sample = ph$sample,
                       group = factor(ifelse(ph$group == "OP", "Case", "Control"),
                                      levels = c("Control", "Case")),
                       age = NA_real_, sex = NA_character_, donor = NA_character_,
                       stringsAsFactors = FALSE))
  lg("  ", paste(dim(e), collapse = " x "), " | ",
     paste(names(table(cohorts[[g]]$pheno$group)),
           table(cohorts[[g]]$pheno$group), collapse = " / "))
}

## =====================================================================
## Sanity checks
## =====================================================================
lg("---- sanity checks ----")
for (nm in names(cohorts)) {
  co <- cohorts[[nm]]
  stopifnot(identical(colnames(co$expr), co$pheno$sample))
  son <- if ("SON" %in% rownames(co$expr)) "yes" else "NO"
  lg(sprintf("  %-10s %s  %5d genes  %3d samples  SON present: %s",
             nm, co$disease, nrow(co$expr), ncol(co$expr), son))
}

saveRDS(cohorts, file.path(DIR_R, "cohorts.rds"))

## cohort description table for the manuscript
desc <- do.call(rbind, lapply(names(cohorts), function(nm) {
  co <- cohorts[[nm]]
  data.frame(Cohort = nm, Disease = co$disease, Tissue = co$tissue,
             Platform = co$platform,
             Control = sum(co$pheno$group == "Control"),
             Case = sum(co$pheno$group == "Case"),
             Genes = nrow(co$expr), stringsAsFactors = FALSE)
}))
fwrite(desc, file.path(ROOT, "05_tables/Table1_cohorts.csv"))
lg("cohort table written")
print(desc)

close(log_con)
