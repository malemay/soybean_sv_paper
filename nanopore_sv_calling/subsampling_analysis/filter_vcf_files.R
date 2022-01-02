#!/prg/R/4.0/bin/Rscript

# Loading the GenomicRanges package
library(GenomicRanges)

# Loading the function
# DEPENDENCY : scripts/filter_sniffles.R
source("../../scripts/filter_sniffles.R")

# Getting the list of vcf files to process
input_files  <- dir(".", pattern = ".*seed[0-9]+.vcf", recursive = TRUE)

for(index in 1:80) {
	input_file <- input_files[index]
	output_file <- sub("\\.vcf", "_filtered.vcf", input_file)

	# Filerting the input file and writing the result to a file with suffix filtered.vcf
	filter_sniffles(vcf.file = input_file,
			keep_hets = FALSE,
			max_N_del = 0,
			keep_homref = FALSE,
			chrom_pattern = "^Gm[0-9]{2}$",
			remove_N_ins = TRUE,
			N_ins_range = 20,
			# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
			refseq = "../../refgenome/Gmax_508_v4.0_mit_chlp.fasta",
			min_sv_size = 50,
			max_sv_size = 500000,
			output_file = output_file,
			keep.ins.seq = TRUE, 
			keep.ref.seq = TRUE, 
			sample.name = "",
			qual.field = "QUAL",
			other.field = '',
			check.inv = FALSE,
			keep.ids = FALSE,
			nocalls = TRUE,
			out.fmt = "gr",
			min.sv.size = 30)
}

# Updating the nanopore_sv_calling/subsampling_analysis/SUBSAMPLE_SNIFFLES_CALLING file for the Makefile to work
file.create("SUBSAMPLE_SNIFFLES_CALLING")

