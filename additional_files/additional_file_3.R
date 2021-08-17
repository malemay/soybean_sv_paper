#!/prg/R/4.0/bin/Rscript

# DEPENDENCY : gene_analysis/bp_under_summary.RData
load("../gene_analysis/bp_under_summary.RData")

# DEPENDENCY : scripts/format_go_csv.R
source("../scripts/format_go_csv.R")

# OUTPUT : additional_file_3.csv
format_go_csv(bp_under_summary, "additional_file_3.csv", go_test = TRUE)

