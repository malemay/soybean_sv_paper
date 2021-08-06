#!/prg/R/4.0/bin/Rscript

# Formatting the additional CSV file with the data on the similarity of terminal inverted repeats

# DEPENDENCY : te_analysis/multiple_alignments/filtered_alignments/tir_similarity_noref.txt
tir_data <- read.table("../te_analysis/multiple_alignments/filtered_alignments/tir_similarity_noref.txt",
		       header = TRUE, stringsAsFactors = FALSE)

# Reformatting the first column to keep only useful data
tir_data$tir_file <- gsub("^.*/", "", tir_data$tir_file)
tir_data$tir_file <- gsub("_mafft_TIR_TSD\\.txt", "", tir_data$tir_file)

# Splitting the contents of the concatenated information into different fields
tir_data$sv_chrom <- sapply(strsplit(tir_data$tir_file, "_"), FUN = `[[`, 1)
tir_data$sv_pos <- as.numeric(sapply(strsplit(tir_data$tir_file, "_"), FUN = `[[`, 2))
tir_data$sv_type <- sapply(strsplit(tir_data$tir_file, "_"), FUN = `[[`, 3)
tir_data$sv_length <- as.numeric(sapply(strsplit(tir_data$tir_file, "_"), FUN = `[[`, 4))

# Reordering and renaming the columns as necessary
tir_data <- tir_data[, c("sv_chrom", "sv_pos", "sv_type", "sv_length", "family", "superfamily", "similarity")]
names(tir_data)[5] <- "soytedb_match_id"
names(tir_data)[7] <- "tir_proportion_matching_nucleotides"
tir_data <- tir_data[order(tir_data$tir_proportion_matching_nucleotides, decreasing = TRUE), ]

# Outputting to file
write.table(tir_data, file = "additional_tir_similarity_file.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

