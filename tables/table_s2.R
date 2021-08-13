#!/prg/R/4.0/bin/Rscript

# Reading the data on the overlap between  SVs and genic features
# DEPENDENCY: gene_analysis/overlap_data.txt
overlap_data <- read.table("../gene_analysis/overlap_data.txt",
			   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gene_table <- as.matrix(table(overlap_data$overlap, overlap_data$svtype))
gene_table <- data.frame(DEL = gene_table[, "DEL"],
			 INS = gene_table[, "INS"])

gene_table[, "DELprop"] <- sprintf("%.1f", gene_table[, "DEL"] / sum(gene_table[, "DEL"]) * 100)
gene_table[, "INSprop"] <- sprintf("%.1f", gene_table[, "INS"] / sum(gene_table[, "INS"]) * 100)

# Loading the data on the widths of different features in the genome
# DEPENDENCY: gene_analysis/feature_widths.RData
load("../gene_analysis/feature_widths.RData")
stopifnot(identical(names(feature_widths), rownames(gene_table)))
gene_table[, "gprop"] <- feature_widths / sum(feature_widths) * 100
gene_table$gprop <- sprintf("%.1f", gene_table$gprop)

# Creating a column holding the name of the feature
gene_table$feature <- rownames(gene_table)


# Applying the same approach, but to the results in non-repeated regions
# DEPENDENCY: gene_analysis/overlap_norepeat.txt
overlap_norepeat <- read.table("../gene_analysis/overlap_norepeat.txt",
			   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

norepeat_table <- as.matrix(table(overlap_norepeat$overlap, overlap_norepeat$svtype))
norepeat_table <- data.frame(nrDEL = norepeat_table[, "DEL"],
			     nrINS = norepeat_table[, "INS"])

norepeat_table[, "nrDELprop"] <- sprintf("%.1f", norepeat_table[, "nrDEL"] / sum(norepeat_table[, "nrDEL"]) * 100)
norepeat_table[, "nrINSprop"] <- sprintf("%.1f", norepeat_table[, "nrINS"] / sum(norepeat_table[, "nrINS"]) * 100)

# DEPENDENCY: gene_analysis/norepeat_widths.RData
load("../gene_analysis/norepeat_widths.RData")
stopifnot(identical(names(norepeat_widths), rownames(norepeat_table)))
norepeat_table[, "nrgprop"] <- norepeat_widths / sum(norepeat_widths) * 100
norepeat_table$nrgprop <- sprintf("%.1f", norepeat_table$nrgprop)

# Binding the two tables together before outputting the results
stopifnot(identical(rownames(norepeat_table), gene_table$feature))
gene_table <- cbind(gene_table, norepeat_table)

# Saving to file
# OUTPUT : tables/table_s2.csv
write.table(gene_table, file = "table_s2.csv", row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

