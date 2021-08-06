#!/prg/R/4.0/bin/Rscript

# DEPENDENCY : gene_analysis/bp_over_summary.RData
load("../gene_analysis/bp_over_summary.RData")

# DEPENDENCY : scripts/format_go_csv.R
source("../scripts/format_go_csv.R")

# OUTPUT : additional_bp_over_file.csv
format_go_csv(bp_over_summary, "additional_bp_over_file.csv", go_test = TRUE)

