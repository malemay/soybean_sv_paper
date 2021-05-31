#!/prg/R/4.0/bin/Rscript

# This code file takes the list of sequences identified as actual TIR sequences
# for various polymorphic TE insertions in the file filtered_alignments/correct_tsd_tir_sequences.txt
# and analyzes the divergence between both TIR sequences per polymorphism and TE family.
# The idea is to have a look at whether some families are less divergent between the
# two repeats, as this might indicate recent or ongoing activity.

# Since the reference genome is less error-prone than the assemblies from the long reads,
# the analysis will be done with and without the TIR sequences as derived from the reference
# and results will be compared before either of the two approaches is chosen for presenting
# the final results.

# Loading the ggplot2 library
library(ggplot2)
library(dplyr)

# Setting the working directory
setwd("filtered_alignments")

# Reading the files and samples with TIR sequences from file
# DEPENDENCY : te_analysis/multiple_alignments/filtered_alignments/correct_tsd_tir_sequences.txt
tir_samples <- scan("correct_tsd_tir_sequences.txt", what = character(), sep = "\n")

# Initalizing a data.frame
tir_df <- data.frame(file = character(),
		     sample = character(),
		     stringsAsFactors = FALSE)

# Looping over all the files
for(i in tir_samples) {
	sstring <- strsplit(i, " ")[[1]]
	# Then looping over all the samples
	for(j in sstring[-1]) {
		ij_row <- data.frame(file = sstring[1],
				     sample = j,
				     stringsAsFactors = FALSE)
		tir_df <- rbind(tir_df, ij_row)
	}
}

# Just a sanity check to see that the sample names were properly typed
table(tir_df$sample)
# 
#        AC2001          ALTA    MAPLE_ISLE  MAPLE_PRESTO            NA    OAC_09_35C   OAC_DRAYTON     OAC_EMBRO  OAC_LAKEVIEW     OAC_MADOC 
#            15            14            12             5             9             5             5             5            12             6 
#    OAC_OXFORD  OAC_PETREL_2  OAC_PRUDENCE OAC_STRATFORD       OT09-03    QS5091.50J           ref        ROLAND 
#             7             6            13            11             3             4            18             7 

# We need to transform NAs explicitly and remove them
tir_df[tir_df$sample == "NA", ] <- NA
sum(!complete.cases(tir_df))
tir_df <- tir_df[complete.cases(tir_df), ]

# Formatting some of the columns in the tir_df data.frame
tir_df$family <- dirname(tir_df$file)
tir_df$superfamily <- substr(tir_df$family, 1, 3)
tir_df$tir_file <- sub("cropped_annotated", "TIR_TSD", tir_df$file)

# I want to get the alignment string for each of these TIR for a given sample and polymorphic TE
tir_df$alignment_string <- NA

for(i in 1:nrow(tir_df)) {
	tirs <- read.table(tir_df[i, "tir_file"], header = TRUE, stringsAsFactors = FALSE)
	tir_string <- tirs[tirs$sample == tir_df[i, "sample"], "mismatches"]
	stopifnot(length(tir_string) == 1)
	tir_df[i, "alignment_string"] <- tir_string
}

# Creating a function that parses an alignment string into a named numeric vector with the sum
# of the number of characters corresponding to each letter
parse_tir_string <- function(tir_string) {
	
	# Extracting the letters and their numbers using regular expressions
	numbers <- regmatches(tir_string, gregexpr("\\d+", tir_string))[[1]]
	letters  <- regmatches(tir_string, gregexpr("[a-zA-Z]+", tir_string))[[1]]
	stopifnot(length(numbers) == length(letters))

	# Computing the sum of each letter
	count <- table(rep(letters, times = numbers))
	output <- as.numeric(count)
	names(output) <- names(count)

	output
}

# Using this function to fill a column in the tir_df with the percent similarity
# (sum of matching (m) nucleotides over the total length of the alignment)
tir_df$similarity <- NA

for(i in 1:nrow(tir_df)) {
	parsed_tir <- parse_tir_string(tir_df[i, "alignment_string"])
	tir_df[i, "similarity"] <- parsed_tir["m"] / sum(parsed_tir)
}

# Saving the tir_df to file
# OUTPUT : te_analysis/multiple_alignments/filtered_alignments/tir_df.txt
write.table(tir_df, file = "tir_df.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# It is more relevant to present the mean similarity for a given event (INS/DEL)
# One analysis will be done by excluding the values computed from the reference
# while another analysis will be done with them
all_means <- tir_df %>% group_by(tir_file) %>% summarise(similarity = mean(similarity))
all_means$family <- dirname(all_means$tir_file)
all_means$superfamily <- substr(all_means$family, 1, 3)

noref_means <- tir_df[tir_df$sample != "ref", ] %>% group_by(tir_file) %>% summarise(similarity = mean(similarity))
noref_means$family <- dirname(noref_means$tir_file)
noref_means$superfamily <- substr(noref_means$family, 1, 3)

# We also need to filter out the events that are actually duplicates:
# DTH_uuu_Gm7_61/Gm03_623884_INS_384_mafft_TIR_TSD.txt
# DTH_uuu_Gm7_61/Gm10_2711253_DEL_423_mafft_TIR_TSD.txt
# DTM_uuu_Gm8_29/Gm08_12914621_DEL_228_mafft_TIR_TSD.txt

all_means <- all_means[!all_means$tir_file %in% c("DTH_uuu_Gm7_61/Gm03_623884_INS_384_mafft_TIR_TSD.txt",
						  "DTH_uuu_Gm7_61/Gm10_2711253_DEL_423_mafft_TIR_TSD.txt",
						  "DTM_uuu_Gm8_29/Gm08_12914621_DEL_228_mafft_TIR_TSD.txt"), ]

noref_means <- noref_means[!noref_means$tir_file %in% c("DTH_uuu_Gm7_61/Gm03_623884_INS_384_mafft_TIR_TSD.txt",
							"DTH_uuu_Gm7_61/Gm10_2711253_DEL_423_mafft_TIR_TSD.txt",
							"DTM_uuu_Gm8_29/Gm08_12914621_DEL_228_mafft_TIR_TSD.txt"), ]

# And then the same plot without the values obtained from the reference
ggplot(noref_means, aes(x = family, y = similarity, color = superfamily)) +
	geom_jitter(size = 4, height = 0, width = 0.15) +
	theme_bw()

# Writing the results to files
# OUTPUT : te_analysis/multiple_alignments/filtered_alignments/tir_similarity_all.txt
write.table(all_means, file = "tir_similarity_all.txt",
	    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# OUTPUT : te_analysis/multiple_alignments/filtered_alignments/tir_similarity_noref.txt
write.table(noref_means, file = "tir_similarity_noref.txt",
	    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

