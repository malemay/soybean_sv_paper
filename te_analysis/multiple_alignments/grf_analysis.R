#!/prg/R/4.0/bin/Rscript

# This code generates the TSD and TIR sequences for each of the polymorphic insertions
# using the GenericRepeatFinder program and then parses the output to identify the
# most likely candidates. These will be manually inspected to identify which ones
# likely do correspond to the real TSD/TIR. Then the presence of TSDs within cultivars
# showing the insertions and those not showing the insertion will be analyzed and
# the number of differences (mismatches and indels) between the two TIRs of the
# sequences will be summarized.

# This function parses the output of the GRF candidate.fasta file
parse_grf <- function() {
	grf_lines <- scan("candidate.fasta", what = character(), sep = "\n", quiet = TRUE)

	if(!length(grf_lines)) return(
				      data.frame(sample = character(),
						 start = numeric(),
						 end = numeric(),
						 mismatches = character(),
						 TSD = character(),
						 stringsAsFactors = FALSE))

	metadata <- grep("^>", grf_lines, value = TRUE)
	sequences <- grep("^>", grf_lines, value = TRUE, invert = TRUE)

	metadata <- sub("^>", "", metadata)
	metadata <- strsplit(metadata, ":")
	metadata <- as.data.frame(do.call("rbind", metadata))
	colnames(metadata) <- c("sample", "start", "end", "mismatches", "TSD")
	metadata$start <- as.numeric(metadata$start)
	metadata$end <- as.numeric(metadata$end)
	metadata$sequence <- sequences
	metadata
}


# This first function launches GenericRepeatFinder in MITE mode on the given file
launch_grf <- function(filepath, grf_path, expected_pos = 60) {
	# Extracting path components
	fasta_dir  <- dirname(filepath)
	fasta_file <- basename(filepath)

	# Changing to the directory
	odir <- setwd(fasta_dir)
	on.exit(setwd(odir), add = TRUE)

	# Preparing the GRF command
	svlen <- as.numeric(strsplit(fasta_file, "_")[[1]][4])
	te_type <- substr(strsplit(fasta_dir, "/")[[1]][2], 1, 3) 
	min_tsd <- ifelse(te_type == "DTH", 3, ifelse(te_type == "DTT", 2, 5))
	max_tsd <- ifelse(te_type == "DTH", 3, ifelse(te_type == "DTT", 2, 10))
	command <- paste0(grf_path, " -i ", fasta_file, " -o . -c 1 -p 30 -s 10 --seed_mismatch 8 ",
			  "--min_tr 12 --min_space ", svlen - 60, " --max_space ", svlen + 20, 
			  " --min_tsd ", min_tsd, " --max_tsd ", max_tsd)
	system(command)

	grf_output <- parse_grf()
	if(!nrow(grf_output)) {
		grf_output$size <- numeric()
		grf_output$expected_size <- numeric()
		grf_output$size_diff <- numeric()
		grf_output$pos_diff <- numeric()
		grf_output$combined_diff <- numeric()
		grf_output$tsd_length <- integer()
		return(grf_output)
	}

	# Adding columns for the difference in expected size, difference in expected position, and sum of both
	# This combined difference will be used to filter for the most likely TIR/TSD combination
	grf_output$size <- grf_output$end - grf_output$start + 1
	grf_output$expected_size <- svlen

	grf_output$size_diff <- abs(grf_output$size - grf_output$expected_size)
	grf_output$pos_diff  <- abs(grf_output$start - expected_pos)
	grf_output$combined_diff <- grf_output$size_diff + grf_output$pos_diff

	# Also adding a column for the length of the TSD ; this one will be used to select the most likely DTM
	grf_output$tsd_length <- nchar(grf_output$TSD)

	# Filtering for specific TSD sequences if the te_type is DTH or DTT
	if(te_type == "DTH") grf_output <- grf_output[grf_output$TSD %in% c("TAA", "TTA"), ]
	if(te_type == "DTT") grf_output <- grf_output[grf_output$TSD == "TA", ]

	grf_output <- grf_output[order(grf_output$tsd_length, grf_output$combined_diff, decreasing = c(TRUE, FALSE)), ]
	grf_output <- grf_output[!duplicated(grf_output$sample), ]

	# Now writing this to file if the data.frame is not empty
	if(nrow(grf_output)) {
		output_file <- sub("insertion.fa$", "TIR_TSD.txt", fasta_file)
		write.table(grf_output, file = output_file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
	}

	return(invisible(NULL))
}

files <- system("ls filtered_alignments/*/*insertion.fa", intern = TRUE)

# Getting the path to the GRF executable from the command line
grf <- commandArgs(trailingOnly = TRUE)[1]

for(i in files) {
	launch_grf(i, grf_path = grf, expected_pos = 60)
}

# Then I will loop over all the TIR_TSD.txt files that were created in order to annotate the
# alignments with the positions of the TSDs

annotate_alignment <- function(filepath) {
	# Extracting path components
	tsd_dir  <- dirname(filepath)
	tsd_file <- basename(filepath)

	# Changing to the directory
	odir <- setwd(tsd_dir)
	on.exit(setwd(odir), add = TRUE)

	# Reading the file
	tsd <- read.table(tsd_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

	# Now I want to read the cropped alignement file that corresponds
	alignment_file <- sub("TIR_TSD", "cropped", tsd_file)
	alignments <- scan(alignment_file, what = character(), sep = "\n", quiet = TRUE)

	# I want to loop over all the lines in the tsd table so I can add asterisks to annotate every sequence
	for(i in 1:nrow(tsd)) {
		# First we get the positions of the TSDs in the continuous sequence
		tsd_length <- tsd[i, "tsd_length"]

		tsd_end1   <- tsd[i, "start"] - 1
		tsd_start1 <- tsd_end1 - tsd_length + 1

		tsd_start2 <- tsd[i, "end"] + 1
		tsd_end2   <- tsd_start2 + tsd_length - 1

		# We then need to translate this as a position in the gapped sequence
		# We get the line number that corresponds to that entry
		line_num <- grep(tsd[i, "sample"], alignments)

		letters_table <- stringr::str_locate_all(alignments[line_num + 1], "[ATGCN]")[[1]]
		tsd_start1 <- letters_table[tsd_start1, 1]
		tsd_end1   <- letters_table[tsd_end1, 1]
		tsd_start2 <- letters_table[tsd_start2, 1]
		tsd_end2   <- letters_table[tsd_end2, 1]
		
		# We add the asterisks at the necessary positions on that line
		# But first we need to add at least as much padding as the length of the sequence
		alignments[line_num] <- paste0(alignments[line_num], paste0(rep(" ", nchar(alignments[line_num + 1])), collapse = ""))
		substr(alignments[line_num], tsd_start1, tsd_end1) <- paste0(rep("*", tsd_end1 - tsd_start1 + 1), collapse = "")
		substr(alignments[line_num], tsd_start2, tsd_end2) <- paste0(rep("*", tsd_end2 - tsd_start2 + 1), collapse = "")

	}

	output_file <- sub("TIR_TSD", "cropped_annotated", tsd_file)
	cat(alignments, file = output_file, sep = "\n")
}

files <- system("ls filtered_alignments/*/*TIR_TSD.txt", intern = TRUE)

for(i in files) {
	annotate_alignment(i)
}

