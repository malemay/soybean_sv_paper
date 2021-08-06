#!/prg/R/4.0/bin/Rscript

# DEPENDENCY : gene_analysis/pfam_under_summary.RData
load("../gene_analysis/pfam_under_summary.RData")

# DEPENDENCY : scripts/format_go_csv.R
source("../scripts/format_go_csv.R")

# OUTPUT : additional_pfam_under_file.csv
format_go_csv(pfam_under_summary, "additional_pfam_under_file.csv", go_test = FALSE)

