#!/prg/R/4.0/bin/Rscript

# DEPENDENCY : gene_analysis/pfam_over_summary.RData
load("../gene_analysis/pfam_over_summary.RData")

# DEPENDENCY : scripts/format_go_csv.R
source("../scripts/format_go_csv.R")

# OUTPUT : additional_file_4.csv
format_go_csv(pfam_over_summary, "additional_file_4.csv", go_test = FALSE)

