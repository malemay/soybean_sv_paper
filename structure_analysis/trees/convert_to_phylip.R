#!/prg/R/4.0/bin/Rscript

# A simple script to convert the two distance files to the format required by Phylip
to_phylip <- function(input) {
	# Reading the input and reformatting the samples column
	distance_matrix <- read.table(input)
	samples <- read.table(paste0(input, ".id"))[[1]]
	samples <- substr(samples, 1, 7)
	samples <- paste0(samples, "   ")

	# Combining the data together and outputting to file
	distance_matrix <- cbind(samples, distance_matrix)

	output <- paste0(input, ".phyl")
	cat("102\n", file = output)

	write.table(distance_matrix, file = output, sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

	invisible(NULL)
}

