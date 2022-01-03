#!/prg/R/4.0/bin/Rscript

# Combining the Illumina metadata from various sources to prepare a supplementary table

# Reading the files needed for the table
# DEPENDENCY : illumina_data/metadata/sra_metadata.csv
sra_metadata   <- read.table("../illumina_data/metadata/sra_metadata.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
# DEPENDENCY : illumina_data/metadata/illumina_sample_ids.txt
cultivars      <- read.table("../illumina_data/metadata/illumina_sample_ids.txt", header = TRUE, stringsAsFactors = FALSE)
# DEPENDENCY : utilities/read_lengths.txt
nreads         <- read.table("../utilities/read_lengths.txt", header = TRUE, stringsAsFactors = FALSE)
# DEPENDENCY : illumina_data/metadata/core_set.txt
core_set       <- read.table("../illumina_data/metadata/core_set.txt", header = FALSE, stringsAsFactors = FALSE)
# DEPENDENCY : utilities/line_ids.txt
nanopore_lines <- read.table("../utilities/line_ids.txt", header = TRUE, stringsAsFactors = FALSE)

# Reading the output from the samtools coverage command for each sample will be slightly more elaborate
# DEPENDENCY : illumina_data/metadata/ILLUMINA_COVERAGE (the coverage computed for each sample using samtools coverage)
samtools_coverage <- list()

for(i in cultivars$ind) {
	samtools_coverage[[i]] <- read.table(paste0("../illumina_data/metadata/coverage/", i, "_coverage.txt"), header = TRUE, stringsAsFactors = FALSE, comment.char = "")
}

# Formatting the data.frames so they all show the same cultivars in the same order
sra_metadata$LibraryName <- sub("_", "", sra_metadata$LibraryName)
sra_metadata <- sra_metadata[order(sra_metadata$LibraryName), ]
stopifnot(identical(sra_metadata$LibraryName, cultivars$ind))
stopifnot(identical(sra_metadata$LibraryName, nreads$ind))

# Some other sanity checks
stopifnot(all(sra_metadata$spots == nreads$nreads / 2))
stopifnot(all(sra_metadata$avgLength == nreads$readlength * 2))

# Formatting the core_set line names
core_set[[1]] <- sub("_", "", core_set[[1]])
stopifnot(all(core_set[[1]] %in% cultivars$ind))

# Computing the sequencing depth of each sample from the output of the samtools coverage command
compute_mean <- function(x) {
	x <- x[grepl("^Gm[0-9]{2}$", x[[1]]), ]
	weighted.mean(x$meandepth, x$endpos)
}

coverage_means <- sapply(samtools_coverage, compute_mean)
stopifnot(identical(names(coverage_means), cultivars$ind))

# Assembling the metadata that we want to show in the table
illumina_metadata <- data.frame(ID = cultivars$ind,
				name = cultivars$name,
				length = nreads$readlength,
				spots = sra_metadata$spots,
				depth = as.numeric(coverage_means),
				core = cultivars$ind %in% core_set[[1]],
				nanopore = cultivars$ind %in% nanopore_lines$id,
				run = sra_metadata$Run,
				stringsAsFactors = FALSE)

# Formatting the metadata columns before writing to a CSV file
illumina_metadata$name <- gsub("_", " ", illumina_metadata$name)
illumina_metadata$name <- sub("AC([^ ])", "AC \\1", illumina_metadata$name)
illumina_metadata$name <- sub("Maple([^ ])", "Maple \\1", illumina_metadata$name)

illumina_metadata$core <- ifelse(illumina_metadata$core, "X", " ")
illumina_metadata$nanopore <- ifelse(illumina_metadata$nanopore, "X", " ")

illumina_metadata$spots <- prettyNum(illumina_metadata$spots, big.mark = ",")
illumina_metadata$depth <- sprintf("%.1f", illumina_metadata$depth)

# Writing the table to disk
write.table(illumina_metadata, file = "table_s7.csv", sep = ";", quote = FALSE, col.names = TRUE, row.names = FALSE)

