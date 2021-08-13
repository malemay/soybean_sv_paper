#!/prg/R/4.0/bin/Rscript

# Preparing the data for table S4 of the manuscript
# DEPENDENCY: gene_analysis/bp_over_summary.RData
load("../gene_analysis/bp_over_summary.RData")

# DEPENDENCY : scripts/format_go_table.R
source("../scripts/format_go_table.R")

# OUTPUT: tables/table_s4.csv
write.table(format_go_table(bp_over_summary), file = "table_s4.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

