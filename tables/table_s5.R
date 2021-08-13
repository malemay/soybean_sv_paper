#!/prg/R/4.0/bin/Rscript

# Preparing the data for table S5 of the manuscript
# DEPENDENCY: gene_analysis/bp_under_summary.RData
load("../gene_analysis/bp_under_summary.RData")

# DEPENDENCY : scripts/format_go_table.R
source("../scripts/format_go_table.R")

# OUTPUT: tables/table_s5.csv
write.table(format_go_table(bp_under_summary), file = "table_s5.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

