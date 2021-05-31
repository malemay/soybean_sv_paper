#/prg/R/4.0/bin/Rscript

library(stringr)

# This code extracts the sequences ± 50 bp from the putative TE insertions
# for each of the multiple alignments that have been kept after filtering.
# The objective is to then use these sequences with GenericRepeatFinder
# to identity the TSD and TIR sequences of each of the insertions.

extract_te_sequences <- function(alignment_file, window_size, output_file) {
	# Reading the alignments from file using the scan() function
	alignment <- scan(alignment_file, what = character(), sep = "\n", quiet = TRUE)
	if(grepl("\\*", alignment[1])) alignment <- alignment[-1]

	# Reformatting the input such that it is a list of character vectors with one element per character
	sequences <- alignment[seq(2, length(alignment), by = 2)]
	sequences <- as.list(sequences)
	
	names(sequences) <- alignment[seq(1, length(alignment), by = 2)]
	names(sequences) <- sub("^>", "", names(sequences))
	stopifnot(names(sequences)[1] == "ref")

	# Get the SVTYPE and length of the SV from the name of the file
	base <- basename(alignment_file)
	svtype <- strsplit(base, "_")[[1]][3]
	svlen  <- as.numeric(strsplit(base, "_")[[1]][4])
	svlen <- ifelse(svtype == "DEL", svlen, 0)
	boundaries <- c(500 - window_size, 500 + svlen + window_size)
	n_pos <- str_locate_all(sequences[["ref"]], "[ATGCN]")[[1]]
	boundaries <- n_pos[boundaries, 1]

	# Extracting the right window for each of the sequences
	sequences <- lapply(sequences, function(x) substr(x, boundaries[1], boundaries[2]))

	# Writing the sequences to file
	if(file.exists(output_file)) file.remove(output_file)

	for(i in 1:length(sequences)) {
		cat(">", names(sequences)[i], "\n", sep = "", file = output_file, append = TRUE)
		cat(sequences[[i]], "\n", sep = "", file = output_file, append = TRUE)
	}

}

# DEPENDENCY : te_analysis/multiple_alignments/mafft_metadata.txt
# DEPENDENCY : filtered multiple alignments
files <- read.table("mafft_metadata.txt", stringsAsFactors = FALSE)
files <- files[complete.cases(files), ]
files <- paste0("filtered_alignments/", files[[1]])

output_files <- sub("\\.txt$", "_cropped.txt", files)

for(i in 1:length(files)) {
	extract_te_sequences(files[i], 60, output_files[i])
}

# Now I want to classify the sequences in those files as TE insertions or deletions
# The number of gaps per sequence should be a good indicator of this. Deletions
# should have > 50% of gaps while insertions should have < 50%

classify_sequences <- function(alignment_file) {
	# Reading the alignments from file using the scan() function
	alignment <- scan(alignment_file, what = character(), sep = "\n", quiet = TRUE)
	if(grepl("\\*", alignment[1])) alignment <- alignment[-1]

	# Reformatting the input such that it is a list of character vectors with one element per character
	sequences <- alignment[seq(2, length(alignment), by = 2)]
	sequences <- as.list(sequences)
	
	names(sequences) <- alignment[seq(1, length(alignment), by = 2)]
	names(sequences) <- sub("^>", "", names(sequences))
	stopifnot(names(sequences)[1] == "ref")

	output <- sapply(sequences, function(x) str_count(x, "-") / nchar(x))
	output <- ifelse(output > 0.6, "DEL", ifelse(output < 0.4, "INS", "unknown"))

	return(output)
}

# Before outputting the insertion sequences for MITE discovery with GRF, I want to
# identify potentially problematic TEs by counting the number of deletions, insertions
# and unknowns for each of the multiple alignments
files <- read.table("mafft_metadata.txt", stringsAsFactors = FALSE)
files <- files[complete.cases(files), ]
files$cropped_file <- sub("\\.txt", "_cropped.txt", paste0("filtered_alignments/", files[[1]]))

files$del <- NA
files$ins <- NA
files$unknown <- NA

for(i in 1:nrow(files)) {
	i_result <- classify_sequences(files[i, "cropped_file"])
	files[i, "del"] <- sum(i_result == "DEL")
	files[i, "ins"] <- sum(i_result == "INS")
	files[i, "unknown"] <- sum(i_result == "unknown")
}

# Let us first have a look at those that have unknown values (proportion of gaps between 40% and 50%)
files[files$unknown > 0, ]
# [1] V1           V2           cropped_file del          ins          unknown     
# <0 rows> (or 0-length row.names)

# Now let us look at the case for which either DEL or INS are 0
files[files$del == 0 | files$ins == 0, ]
#                                                  V1 V2                                                                cropped_file del ins unknown
# 148 DTM_uuu_Gm18_11/Gm12_4854493_DEL_1459_mafft.txt  2 filtered_alignments/DTM_uuu_Gm18_11/Gm12_4854493_DEL_1459_mafft_cropped.txt   0   7       0
# 167 DTM_uuu_Gm20_8/Gm13_27877010_DEL_1473_mafft.txt  5 filtered_alignments/DTM_uuu_Gm20_8/Gm13_27877010_DEL_1473_mafft_cropped.txt   0  11       0

# The polymorphic sequences of these were filtered out because their alignments had too
# many differences relative to the reference, so it's okay. We will remove them from the
# analysis
potential_mites <- files[files$del > 0 & files$ins > 0, ]

# Now a function that identifies the insertion sequences from the alignments and outputs
# them to a fasta file that will be used for input with GenericRepeatFinder
write_te_fasta <- function(alignment_file, output_file) {

	# Reading the sequences but then keeping only those that contain the insertion
	alignment <- scan(alignment_file, what = character(), sep = "\n", quiet = TRUE)
	if(grepl("\\*", alignment[1])) alignment <- alignment[-1]

	# Reformatting the input such that it is a list of character vectors with one element per character
	sequences <- alignment[seq(2, length(alignment), by = 2)]
	sequences <- as.list(sequences)
	
	names(sequences) <- alignment[seq(1, length(alignment), by = 2)]
	names(sequences) <- sub("^>", "", names(sequences))
	stopifnot(names(sequences)[1] == "ref")

	# Getting the classification of the sequences
	classification <- classify_sequences(alignment_file)
	insertions <- names(classification[classification == "INS"])
	insertions <- sequences[insertions]

	# Removing the gaps before writing to file
	insertions <- lapply(insertions, function(x) gsub("-+", "", x))

	# Writing the result as a fasta file
	if(file.exists(output_file)) file.remove(output_file)

	for(i in 1:length(insertions)) {
		cat(">", names(insertions)[i], "\n", sep = "", file = output_file, append = TRUE)
		cat(insertions[[i]], "\n", sep = "", file = output_file, append = TRUE)
	}

	invisible(NULL)
}

for(i in 1:nrow(potential_mites)) {
	input_file  <- potential_mites[i, "cropped_file"]
	output_file <- sub("cropped\\.txt", "insertion.fa", input_file)
	write_te_fasta(input_file, output_file)
}

