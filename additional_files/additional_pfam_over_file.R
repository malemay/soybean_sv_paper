#!/prg/R/4.0/bin/Rscript

# DEPENDENCY : gene_analysis/pfam_over_summary.RData
load("../gene_analysis/pfam_over_summary.RData")

# DEPENDENCY : scripts/format_go_csv.R
source("../scripts/format_go_csv.R")

# OUTPUT : additional_pfam_over_file.csv
format_go_csv(pfam_over_summary, "additional_pfam_over_file.csv", go_test = FALSE)

