
revcomp <- function(sequence) {

	# Checking that only one sequence is provided
	stopifnot(length(sequence) == 1)

	# The lookup table that will be used for replacement
	rep_table <- c("A" = "T",
		       "T" = "A",
		       "G" = "C",
		       "C" = "G",
		       "N" = "N")

	# Splitting the sequence into its constituent nucleotides
	sequence <- strsplit(sequence, "")[[1]]

	# Replacing the nucleotides by their complement
	sequence <- rep_table[sequence]

	# Returning the inverted sequence
	paste0(rev(sequence), collapse = "")
}

