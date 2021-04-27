# File initially created on April 19, 2021

# The code in this file will test the pairwise differences between the allele
# frequencies of deletions and insertions overlapping various genic features
# using a permutation approach.

# The "overlap" column of the data object will be randomly shuffled 10,000
# times and the mean for each group computed each time. The comparison of the
# observed differences in allele frequency and the random distrbution of those
# differences will allow a p-value to be attributed to each of the 12 pairwise
# differences between groups (6 for deletions, and 6 for insertions). A correction
# will be applied to compensate for multiple testing.

# Loading the data.frame of allele frequencies and genic feature overlap
# DEPENDENCY : overlap_data.txt
overlap_data <- read.table("overlap_data.txt", header = TRUE, stringsAsFactors = FALSE)

# Removing the entries with allele frequences = 1 because they might represent errors in the reference
overlap_data <- overlap_data[overlap_data$af != 1, ]

# Splitting the data into separate data.frames for deletions and insertions
deletions  <- overlap_data[overlap_data$svtype == "DEL", ]
insertions <- overlap_data[overlap_data$svtype == "INS", ]

# Computing the observed differences in mean allele frequencies between overlap groups for deletions
deletion_means <- tapply(deletions$af, deletions$overlap, mean)
deletion_differences <- numeric(0)

for(i in 1:(length(deletion_means) - 1)) {
	for(j in (i + 1):(length(deletion_means))) {
		i_name <- names(deletion_means)[i]
		j_name <- names(deletion_means)[j]
		deletion_differences <- c(deletion_differences, deletion_means[i] - deletion_means[j])
		names(deletion_differences)[length(deletion_differences)] <- paste0(i_name, "_", j_name)
	}
}

# Now computing the same differences for insertions
insertion_means <- tapply(insertions$af, insertions$overlap, mean)
insertion_differences <- numeric(0)

for(i in 1:(length(insertion_means) - 1)) {
	for(j in (i + 1):(length(insertion_means))) {
		i_name <- names(insertion_means)[i]
		j_name <- names(insertion_means)[j]
		insertion_differences <- c(insertion_differences, insertion_means[i] - insertion_means[j])
		names(insertion_differences)[length(insertion_differences)] <- paste0(i_name, "_", j_name)
	}
}

# Computing the mean frequencies over 10000 permutations for both deletions and insertions
deletion_permutations  <- matrix(nrow = 10000, ncol = 4)
insertion_permutations <- matrix(nrow = 10000, ncol = 4)

for(i in 1:10000) {

	if(i %% 100 == 0) message("Iteration ", i)

	# First computing for deletions
	deletions$shuffled_overlap <- sample(deletions$overlap)
	i_deletions <- tapply(deletions$af, deletions$shuffled_overlap, mean)
	stopifnot(identical(names(i_deletions), c("cds", "gene", "intergenic", "upstream5kb")))
	deletion_permutations[i, ] <- as.numeric(i_deletions)

	# Then computing for insertions
	insertions$shuffled_overlap <- sample(insertions$overlap)
	i_insertions <- tapply(insertions$af, insertions$shuffled_overlap, mean)
	stopifnot(identical(names(i_insertions), c("cds", "gene", "intergenic", "upstream5kb")))
	insertion_permutations[i, ] <- as.numeric(i_insertions)
}

deletion_permutations <- as.data.frame(deletion_permutations)
names(deletion_permutations) <- c("cds", "gene", "intergenic", "upstream5kb")
insertion_permutations <- as.data.frame(insertion_permutations)
names(insertion_permutations) <- c("cds", "gene", "intergenic", "upstream5kb")

# Saving some of the objects in a RData file
save(deletion_means, deletion_differences, deletion_permutations, 
     insertion_means, insertion_differences, insertion_permutations,
     file = "allele_frequency_permutations.RData")

