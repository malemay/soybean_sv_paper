#!/prg/R/4.0/bin/Rscript

# Benchmarking the variants found by each of the Sniffles discovery steps based on
# subsampled .bam files using the sveval package.

# Here, the truth set consists in the calls made for each sample PRIOR to breakpoint refinement

# Four samples with average sequencing depth > 16 were chosen for this analysis:
# AC2001, ALTA, MAPLE_ISLE, and OAC_LAKEVIEW

# Each sample was subsampled at a fraction of 0.2, 0.4, 0.6, and 0.8 of the original
# coverage 5 times using different seeds for each sampling

# Loading the sveval library (modified by me to allow for argument output.all)
library(sveval)

# Sourcing the extract_rates.R file which contains the fonction to gather precision/sensitivity rates
# DEPENDENCY : scripts/extract_rates.R
source("../../scripts/extract_rates.R")

# And another function that reads a vcf file into a GRanges object and filters it according to various criteria
# DEPENDENCY : scripts/read_filter_vcf.R
source("../../scripts/read_filter_vcf.R")

# Creating a vector of sample names, files to process, and truth files
samples <- c("AC2001", "ALTA", "MAPLE_ISLE", "OAC_LAKEVIEW")

# DEPENDENCY : nanopore_sv_calling/subsampling_analysis/SUBSAMPLE_SNIFFLES_CALLING (calls on subsampled reads)
input_files <- dir(".", pattern = ".*seed[0-9]+_filtered.vcf", recursive = TRUE)
stopifnot(length(input_files) == 80)
input_samples <- sub("/.*", "", input_files)

# DEPENDENCY : nanopore_sv_calling/SV_CALLING (filtered SV calls, not breakpoint-refined)
truth_files <- paste0("../", samples, "_hom70_filtered.vcf")
stopifnot(all(file.exists(truth_files)))
names(truth_files) <- samples

# Extracting the relevant information for this task
for(slurm_index in 1:80) {
	input_file <- input_files[slurm_index]
	input_sample <- input_samples[slurm_index]
	truth_file <- truth_files[input_sample]
	i <- sub(".*/(.*)_filtered.vcf", "\\1", input_file)
	j <- "sniffles"

	# Reading the vcf file for the truth variants (unsampled file)
	truth <- read_filter_vcf(vcf.file = truth_file, 
				 keep_hets = FALSE, max_N_del = 0, 
				 refseq = "../../refgenome/Gmax_508_v4.0_mit_chlp.fasta",
				 keep_homref = FALSE, chrom_pattern = "^Gm[0-9]{2}$", 
				 remove_N_ins = TRUE, keep.ins.seq = TRUE, 
				 keep.ref.seq = TRUE, sample.name = "",
				 qual.field = "QUAL", check.inv = FALSE, 
				 keep.ids = FALSE, nocalls = FALSE, 
				 out.fmt = "gr", min.sv.size = 50)

	# Keeping the user informed
	message("Processing ", i, " --- ", j)

	# Reading the vcf file for the Illumina variants
	calls <- read_filter_vcf(vcf.file = input_file, 
				 keep_hets = FALSE, max_N_del = NULL, 
				 keep_homref = FALSE, chrom_pattern = "^Gm[0-9]{2}$", 
				 remove_N_ins = TRUE, keep.ins.seq = TRUE, 
				 keep.ref.seq = TRUE, sample.name = i,
				 qual.field = "DV", check.inv = FALSE, 
				 keep.ids = FALSE, nocalls = FALSE,
				 out.fmt = "gr", min.sv.size = 50)

	# Calling svevalOl on that file
	sveval_output  <- svevalOl(calls.gr = calls,
				   truth.gr = truth,
				   ins.seq.comp = TRUE,
				   min.size = 50,
				   qual.field = "QUAL",
				   sample.name = "",
				   qual.ths = c(2, 3),
				   check.inv = FALSE,
				   geno.eval = FALSE,
				   stitch.hets = FALSE,
				   merge.hets = FALSE,
				   nb.cores = 1,
				   log.level = "INFO",
				   output.all = TRUE)

	# The output file is saved to an RData file
	save(sveval_output, file = paste0("benchmark_files/", i, "_", j, ".RData"))

	# Processing the results with extract_rates
	ij_rates  <- extract_rates(sveval_output, c("DEL", "INS", "INV", "DUP"),
				   cultivar = i, pipeline = j)	

	saveRDS(ij_rates, file = paste0("benchmark_files/", i, "_", j, "_rates.rds"))
}

# Marking this step as done for the Makefile
file.create("SUBSAMPLE_BENCHMARKS")

