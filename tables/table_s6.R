#!/prg/R/4.0/bin/Rscript

# Preparing the data for table S6 of the manuscript

# DEPENDENCY: gene_analysis/pfam_over_summary.RData
load("../gene_analysis/pfam_over_summary.RData")

# DEPENDENCY : scripts/format_go_table.R
source("../scripts/format_go_table.R")

# OUTPUT: tables/table_s6.csv
write.table(format_go_table(pfam_over_summary), file = "table_s6.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

