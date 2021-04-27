#!/prg/R/4.0/bin/Rscript

# Loading the required packages
library(ggplot2)

# Reading the data on the overlap bet
# DEPENDENCY : overlap_data.txt
overlap_data <- read.table("../gene_analysis/overlap_data.txt",
			   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gene_table <- as.matrix(table(overlap_data$overlap, overlap_data$svtype))
gene_table <- data.frame(DEL = gene_table[, "DEL"],
			 INS = gene_table[, "INS"])

gene_table[, "DELprop"] <- round(gene_table[, "DEL"] / sum(gene_table[, "DEL"]) * 100, 1)
gene_table[, "INSprop"] <- round(gene_table[, "INS"] / sum(gene_table[, "INS"]) * 100, 1)

# Loading the data on the widths of different features in the genome
# DEPENDENCY : feature_wdiths.RData
load("../gene_analysis/feature_widths.RData")
stopifnot(identical(names(feature_widths), rownames(gene_table)))
gene_table[, "gprop"] <- feature_widths / sum(feature_widths) * 100
gene_table$gprop <- sprintf("%.1f", gene_table$gprop)

# Creating a column holding the name of the feature
gene_table$feature <- rownames(gene_table)

# Saving to file
write.table(gene_table, file = "table_3.csv", row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

