
# This function takes an alignment file as input and writes a file for which
# only alignments with a high percent identity (e.g. > 0.75) to the reference
# is observed in terminal regions (first and last 500 positions of the alignment)
filter_alignments <- function(alignment_file, window_size, threshold, output_file) {
	# Reading the alignments from file using the scan() function
	alignment <- scan(alignment_file, what = character(), sep = "\n", quiet = TRUE)
	if(grepl("\\*", alignment[1])) alignment <- alignment[-1]

	# Reformatting the input such that it is a list of character vectors with one element per character
	sequences <- alignment[seq(2, length(alignment), by = 2)]
	sequences <- as.list(sequences)
	
	names(sequences) <- alignment[seq(1, length(alignment), by = 2)]
	names(sequences) <- sub("^>", "", names(sequences))
	stopifnot(names(sequences)[1] == "ref")

	sequences <- lapply(sequences, function(x) strsplit(x, "")[[1]])

	# Computing the percent identity of the first and last window_size positions in the alignment 
	ref <- sequences[[1]]
	slen <- length(ref)
	stopifnot(length(unique(unlist(lengths(sequences)))) == 1)

	pass <- lapply(sequences, function(x) {
			       first_window <- sum(ref[1:window_size] == x[1:window_size])
			       last_window  <- sum(ref[(slen - window_size + 1):slen] == x[(slen - window_size + 1):slen])
			       pident <- c(first_window, last_window) / window_size 
			       all(pident > threshold) })

	# Using the output to filter the sequences to write to file
	pass <- unlist(pass)
	sequences <- sequences[pass]

	if(file.exists(output_file)) file.remove(output_file)

	for(i in 1:length(sequences)) {
		cat(">", names(sequences)[i], "\n", sep = "", file = output_file, append = TRUE)
		cat(paste0(sequences[[i]], collapse = ""), "\n", sep = "", file = output_file, append = TRUE)
	}

	invisible(NULL)
}

# Reading the data with the number of ALT allele DNA TEs per alignment
# DEPENDENCY : te_analysis/multiple_alignments/mafft_metadata.txt
files <- read.table("mafft_metadata.txt", stringsAsFactors = FALSE)
names(files) <- c("path", "n_alt")
files <- files[complete.cases(files), ]
files$dir <- sub("/.*", "", files$path)
files$filename <- sub(".*/", "", files$path)

for(i in 1:nrow(files)) {
	output_dir <- paste0("filtered_alignments/", files[i, "dir"])
	if(!dir.exists(output_dir)) dir.create(output_dir)
	output_file <- paste0(output_dir, "/", files[i, "filename"])
	filter_alignments(alignment_file = files[i, "path"],
			  window_size = 500, threshold = 0.75, 
			  output_file = output_file)
}

