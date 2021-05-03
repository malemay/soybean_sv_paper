#!/prg/R/4.0/bin/R

# Loading the data from the allele frequency randomization analysis
# DEPENDENCY: gene_analysis/allele_frequency_permutations.RData
load("../gene_analysis/allele_frequency_permutations.RData")

# Computing the differences between the values for deletions
for(i in 1:3) {
	for(j in (i + 1):4) {
		i_name <- names(deletion_permutations)[i]
		j_name <- names(deletion_permutations)[j]
		deletion_permutations[[paste0(i_name, "_", j_name)]] <- deletion_permutations[[i]] - deletion_permutations[[j]]
	}
}

# Computing unilateral p-values
deletion_pvalues <- deletion_differences
for(i in names(deletion_pvalues)) {
	if(deletion_differences[i] < 0) {
		deletion_pvalues[i] <- sum(deletion_differences[i] > deletion_permutations[[i]]) / nrow(deletion_permutations)
	} else {
		deletion_pvalues[i] <- sum(deletion_differences[i] < deletion_permutations[[i]]) / nrow(deletion_permutations)
	}
}

# Computing the differences between the values for insertions
for(i in 1:3) {
	for(j in (i + 1):4) {
		i_name <- names(insertion_permutations)[i]
		j_name <- names(insertion_permutations)[j]
		insertion_permutations[[paste0(i_name, "_", j_name)]] <- insertion_permutations[[i]] - insertion_permutations[[j]]
	}
}

# Computing unilateral p-values
insertion_pvalues <- insertion_differences
for(i in names(insertion_pvalues)) {
	if(insertion_differences[i] < 0) {
		insertion_pvalues[i] <- sum(insertion_differences[i] > insertion_permutations[[i]]) / nrow(insertion_permutations)
	} else {
		insertion_pvalues[i] <- sum(insertion_differences[i] < insertion_permutations[[i]]) / nrow(insertion_permutations)
	}
}

# Putting all the results into a table for writing to a csv file and preparing the LaTeX table
output_table <- data.frame(regions = names(deletion_differences),
			   stringsAsFactors = FALSE)

# Adding the results to the columns
output_table$del_difference <- deletion_differences[output_table$regions]
output_table$ins_difference <- insertion_differences[output_table$regions]
output_table$del_pvalue <- deletion_pvalues[output_table$regions]
output_table$ins_pvalue <- insertion_pvalues[output_table$regions]

# Formatting the columns for the final table
output_table$regionfmt <- gsub("_", " - ", output_table$regions)
output_table$deldiff <- sprintf("%.4f", output_table$del_difference)
output_table$insdiff <- sprintf("%.4f", output_table$ins_difference)
output_table$delp <- ifelse(output_table$del_pvalue == 0, "$< 10^{-4}$", sprintf("%.3f", output_table$del_pvalue))
output_table$insp <- ifelse(output_table$ins_pvalue == 0, "$< 10^{-4}$", sprintf("%.3f", output_table$ins_pvalue))

# Outputting the table to file
# OUTPUT: tables/table_s5.csv
write.table(output_table, file = "table_s5.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")

