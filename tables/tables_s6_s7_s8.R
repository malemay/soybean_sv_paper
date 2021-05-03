#!/prg/R/4.0/bin/Rscript

# Preparing the data for Tables S6 to S8 of the manuscript
# DEPENDENCY: gene_analysis/bp_over_summary.RData
load("../gene_analysis/bp_over_summary.RData")
# DEPENDENCY: gene_analysis/bp_under_summary.RData
load("../gene_analysis/bp_under_summary.RData")
# DEPENDENCY: gene_analysis/pfam_over_summary.RData
load("../gene_analysis/pfam_over_summary.RData")

# Writing a function that will be used to format each of the tables since they are based on the same
format_go_table <- function(x){
	x <- x[x$Pvalue < 0.05, ]
	x$Pvalue <- sprintf("%.4g", x$Pvalue)
	x$OddsRatio <- sprintf("%.2f", x$OddsRatio)
	x$ExpCount <- sprintf("%.2f", x$ExpCount)
	x
}

# OUTPUT: tables/table_s6.csv
write.table(format_go_table(bp_over_summary), file = "table_s6.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
# OUTPUT: tables/table_s7.csv
write.table(format_go_table(bp_under_summary), file = "table_s7.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
# OUTPUT: tables/table_s8.csv
write.table(format_go_table(pfam_over_summary), file = "table_s8.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

