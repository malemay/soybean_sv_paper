#!/prg/R/4.0/bin/R

# Loading the data from the allele frequency randomization analysis
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

