###############################################################################
# 14_combine_paper.R
###############################################################################

suppressPackageStartupMessages({ library(officer); library(magrittr) })
setwd("D:/bioinfo05")

man <- read_docx("06_paper/Manuscript.docx") %>%
  body_add_break() %>%
  body_add_docx(src = "06_paper/Tables.docx")

print(man, target = "06_paper/Full_Paper.docx")
cat("Full_Paper.docx generated\n")
