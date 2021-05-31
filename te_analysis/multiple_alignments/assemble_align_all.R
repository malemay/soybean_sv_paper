#!/prg/R/4.0/bin/Rscript

# Loading the required libraries
library(Rsamtools)
library(Biostrings)

# First we need to define some functions that will be used as part of the analysis

# A function that parses a line in a VCF file in order to extract metadata
# about the structural variants it represents
parse_svinfo <- function(vcf_line) {

  stopifnot(length(vcf_line) == 1)

  # Splitting the vcf line into its fields
  split_line <- strsplit(vcf_line, "\\s+")[[1]]

  # Assigning some elements directly
  chr <- split_line[1]
  start <- as.numeric(split_line[2])
  ref <- split_line[4]
  alt <- split_line[5]

  # Some require a bit of processing
  svtype <- sub(".*SVTYPE=([[:upper:]]+).*", "\\1", split_line[8])
  svlen  <- ifelse(svtype == "INS", 
		   nchar(alt) - nchar(ref),
		   nchar(ref) - nchar(alt))

  list(chr = chr, svtype = svtype, start = start, svlen = svlen)
}


# This function accepts arguments from update_breakpoints() in order
# to launch the SV realignment pipeline which involves the following steps:

# - Assembly of the Nanopore reads from the region of interest with wtdbg2
# - Polishing of the resulting assembly using the Nanopore reads themselvs
#   aligned to the assembled contig with minimap2
# - Alignement of the resulting contig to the reference genome using AGE

# This function is merely a wrapper around the bash script age_realign.sh

# The function creates some temporary files that are used by the bash script.
# It also manages those files by removing those that are no longer used at the
# end and returning the names of the alignment file (.age.txt) and contig fasta
# file (.ontpolish.fa) to the function that calls it so it knows where to find them

# This function also creates the command line that will be parsed using getopts
# by the bash script.

call_age <- function(chr, svtype, start, svlen, reference_window, reads_window,
		     age_script, samtools, minimap2, wtdbg2, wtpoa_cns, nanopore_bam,
		     output_dir) {

	# Creating temporary files that will be used by the bash script

	# The file to which the Nanopore reads will be written
	reads_fasta <- tempfile("reads_fasta", tmpdir = getwd(), fileext = ".fa")
	file.create(reads_fasta)
	on.exit(file.remove(reads_fasta), add = TRUE)

	# The file to which the assembled contig will be written
	contig_fasta <- tempfile("contig_fasta", tmpdir = getwd(), fileext = ".fa")
	file.create(contig_fasta)
	on.exit(file.remove(contig_fasta), add = TRUE)

	# The file to which the polished contig will be written
	polished_fasta <- paste0(output_dir, "/", chr, "_", start, "_", svtype, "_", abs(svlen), "_", nanopore_bam, ".fa")
	file.create(polished_fasta)

	# The prefix that fill be used
	prefix <- tempfile("prefix", tmpdir = "")
	prefix <- sub("/", "", prefix)

	command <- paste0(age_script,
			  " -s ", samtools,
			  " -m ", minimap2,
			  " -w ", wtdbg2,
			  " -c ", wtpoa_cns,
			  " -t ", svtype,
			  " -h ", chr,
			  " -e ", start,
			  " -l ", svlen,
			  " -b ", paste0(nanopore_bam, ".bam"),
			  " -n ", reads_fasta,
			  " -f ", contig_fasta,
			  " -p ", polished_fasta,
			  " -y ", prefix,
			  " -u ", reads_window)

	# Launching the command itself
	system(command)

	# Returning the names of the files to be used by update_breakpoints()
	return(invisible(NULL))
}

# A function that assembles contigs from Oxford Nanopore data using
# SVs that are stored in a vcf file
assemble_contigs <- function(vcf_file, samples, reference_window, reads_window,
			     age_script, samtools, minimap2, wtdbg2, wtpoa_cns,
			     assembly_dir, refgenome, pos_window, ginsi_path) {
	vcf_lines <- scan(vcf_file, what = character(), sep = "\n", quiet = TRUE)
	vcf_lines <- grep("^#", vcf_lines, invert = TRUE, value = TRUE)
	vcf_lines

	for(i in vcf_lines) {
		svinfo <- parse_svinfo(i)

		for(sample in samples) {
			message(svinfo$chr, " ", svinfo$start, " ", svinfo$svtype, " ", svinfo$svlen, " ", sample)
			call_age(chr = svinfo$chr,
				 svtype = svinfo$svtype,
				 start = svinfo$start,
				 svlen = svinfo$svlen,
				 reference_window = reference_window,
				 reads_window = reads_window, 
				 age_script = age_script,
				 samtools = samtools,
				 minimap2 = minimap2,
				 wtdbg2 = wtdbg2,
				 wtpoa_cns = wtpoa_cns,
				 nanopore_bam = sample,
				 output_dir = assembly_dir)
		}

		message("Performing multiple alignment ", svinfo$chr, " ", svinfo$start, " ", svinfo$svtype, " ", svinfo$svlen)
		multiple_alignment(chrom = svinfo$chr, pos = svinfo$start,
				   svtype = svinfo$svtype, svlen = svinfo$svlen,
				   assembly_dir = assembly_dir, refgenome = refgenome,
				   pos_window = pos_window, ginsi_path = ginsi_path)

	}

	invisible(NULL)
}


# This function aligns local assemblies against a reference sequence in order to
# extract only the relevant subsequences. These can then be aligned using multiple alignment
extract_query_sequences <- function(query, ref, output_fasta) {

	# Computing the pairwise alignment of each query against the reference
	alignment <- Biostrings::pairwiseAlignment(query, ref, type = "local-global")

	# Then looping to write the aligned part of the query sequences to file
	stopifnot(length(query) == length(alignment))

	# First writing the reference sequence to file, overwriting it at the same time
	cat(">ref\n", as.character(ref), "\n", sep = "", file = output_fasta)

	for(i in 1:length(alignment)) {
		cat(">", names(query)[i], "\n", as.character(aligned(alignment[i], degap = TRUE)), "\n", 
		    sep = "", file = output_fasta, append = TRUE)
	}

	return(invisible(NULL))
}


# This function takes a putative TE position as input and generates the mutliple alignment from them
multiple_alignment <- function(chrom, pos, svtype, svlen, assembly_dir, refgenome, pos_window, ginsi_path) {
	# First generating the pattern to look for in files
	pattern <- paste0(chrom, "_", as.character(pos), "_", svtype, "_", as.character(svlen))

	# Get the names of the files to get the local assemblies from
	files <- dir(assembly_dir, pattern = pattern)
	files <- paste0(assembly_dir, "/", files)

	# Calculating a reference genome sv size for reference extraction purposes
	refsize <- ifelse(svtype == "DEL", svlen, 1)

	# Getting the reference sequence ± pos_window nt from the position of interest
	ref <- Rsamtools::scanFa(refgenome,
				 param = GenomicRanges::GRanges(seqnames = chrom,
								 ranges = IRanges::IRanges(start = pos - pos_window,
											   end   = pos + pos_window + refsize)))

	# Initializing a fasta file to concatenate all the assemblies
	multifasta <- paste0(assembly_dir, "/", pattern, ".fa")

	if(file.exists(multifasta)) {
		file.remove(multifasta)
	}

	# Concatenating them all in a single file
	for(i in 1:length(files)) {
		i_fasta <- Rsamtools::scanFa(files[i])
		line_name <- sub("\\.fa", "", sub(paste0(".*", pattern, "_"), "", files[i]))
		if(length(i_fasta)) {
			cat(">", line_name, "\n", sep = "", file = multifasta, append = TRUE)
			cat(as.character(i_fasta), "\n", sep = "", file = multifasta, append = TRUE)
		}
	}

	# Reading the queries as a DNAStringSet from the file we just wrote to disk
	queries <- scanFa(multifasta)

	# Aligning the query sequences against the reference by pairwise alignment
	# This is done to extract only the relevant sequences for multiple alignment
	sequences_fasta <- paste0(assembly_dir, "/", pattern, "_sequences.fa")
	extract_query_sequences(queries, ref, sequences_fasta)

	# Performing the alignment with mafft
	mafft_fa   <- paste0(assembly_dir, "/", pattern, "_mafft.fa")
	mafft_file <- paste0(assembly_dir, "/", pattern, "_mafft.txt")
	system(paste0(ginsi_path, " --reorder ", sequences_fasta, " > ", mafft_fa))
	reformat_mafft(mafft_fa, mafft_file)

	# I want to add an asterisk indicating where the polymorphic TE should start in the alignment files
	# First with the mafft file
	mafft_ref <- scan(mafft_file, what = character(), sep = "\n", quiet = TRUE, n = 2)[2]
	te_pos <- stringr::str_locate_all(mafft_ref, "[ATGCN]")[[1]][c(pos_window + 1, pos_window + 1 + refsize), 1]
	mafft_lines <- scan(mafft_file, what = character(), sep = "\n", quiet = TRUE)

	starline <- paste0(rep(" ", max(te_pos)), collapse = "")
	substr(starline, te_pos[1], te_pos[1]) <- "*"
	substr(starline, te_pos[2], te_pos[2]) <- "*"
	cat(starline, "\n", sep = "", file = mafft_file)
	cat(mafft_lines, sep = "\n", file = mafft_file, append = TRUE)
}

# A function to reformat the output of MAFFT
reformat_mafft <- function(input_file, output_file) {
	input <- scanFa(input_file)
	#label_lines <- grep("^>", input)
	#start_lines <- label_lines + 1
	#end_lines   <- c(label_lines[-1] - 1, length(input))

	# Making sure the output file is empty before we start writing to it
	if(file.exists(output_file)) {
		file.remove(output_file)
	}

	for(i in 1:length(input)) {
		cat(">", names(input)[i], "\n", file = output_file, sep = "", append = TRUE)
		cat(as.character(input)[i], "\n",
		    sep = "", file = output_file, append = TRUE)
	}

	return(invisible(NULL))
}

# Creating a vector with the samples to be analyzed
# DEPENDENCY : Oxford Nanopore reads aligned with NGMLR
samples <- gsub("\\.bam", "", dir(pattern = ".*bam$"))

# Getting a vector fo the VCF files over which to loop
# DEPENDENCY : VCF files with the polymorphic DNA TEs for which to assemble the SV region
vcf_files <- dir(pattern = "\\.vcf")

# Getting the paths to the executables from the command line
executables <- commandArgs(trailingOnly = TRUE)
samtools  <- executables[1]
minimap2  <- executables[2]
wtdbg2    <- executables[3]
wtpoa_cns <- executables[4]
ginsi     <- executables[5]

for(i in vcf_files) {
	
	dir.create(sub("\\.vcf", "", i))

	assemble_contigs(i, 
			 samples = samples,
			 reference_window = 500,
			 reads_window = 200, 
			 # DEPENDENCY : te_analysis/multiple_alignments/age_realign.sh
			 age_script = "./age_realign.sh",
			 samtools = samtools,
			 minimap2 = minimap2,
			 wtdbg2 = wtdbg2,
			 wtpoa_cns = wtpoa_cns,
			 assembly_dir = sub("\\.vcf", "", i),
			 # DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
			 refgenome = "../../refgenome/Gmax_508_v4.0_mit_chlp.fasta",
			 pos_window = 500,
			 ginsi_path = ginsi)
}

