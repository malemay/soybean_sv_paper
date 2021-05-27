#!/usr/bin/Rscript

# File initially created on September 24, 2020

# Benchmarking the structural variants called by Illumina from the variants discovered by Oxford Nanopore
#  on the fourth version of the soybean reference genome assembly with sveval. Three
#  different methods (BayesTyper, vg, Paragraph) will be benchmarked.

# The variants benchmarked here have not had their breakpoints refined using AGE
# A different benchmark will be done for these

# Structural variants will be filtered post-reading with read_sv_filter (a wrapper around readSVvcf)
#  for the following criteria :
#  - Removal of variants that are not insertions or deletions
#  - Optional removal of heterozygous, missing, and homozygous reference genotype calls
#  - Removal of "*" alleles (BayesTyper)
#  - Removal of variants located on unmapped scaffolds or organellar genomes
#  - Removal of insertion alternate alleles containing N nucleotides or unknown sequence (<INS>)
#  - Removal of insertion reference alleles with > 1 N nucleotide (to accomodate for Sniffles'
#     default "N" reference allele for insertions)
#  - Removal of deletions whose sequence contains > 10% of N (the threshold can be tuned)

# Only the benchmark with geno.eval = FALSE is done because heterozygous variants
# were not used for genotyping with Illumina

# Eventually, inversions and duplications will also be evaluated. However, for this
# first benchmark effort we do not need that much detail.

# Loading the sveval library (modified by me to allow for argument output.all)
library(sveval)

# Sourcing the extract_rates.R file which contains the fonction to gather precision/sensitivity rates
# DEPENDENCY : scripts/extract_rates.R
source("../../../scripts/extract_rates.R")

# And another function that reads a vcf file into a GRanges object and filters it according to various criteria
# DEPENDENCY : scripts/read_filter_vcf.R
source("../../../scripts/read_filter_vcf.R")

# This vector will hold the correspondence between Illumina CAD IDs and line names
line_names <- c("CAD1010" = "MAPLE_PRESTO",
		"CAD1052" = "OAC_EMBRO",
		"CAD1064" = "OAC_CARMAN",
		"CAD1070" = "QS5091.50J")

# Initializing a list to store the results of each iteration
rates_list <- list()

# Looping over all 4 samples
for(i in names(line_names)) {

	# Reading the vcf file for the Nanopore variants
	# DEPENDENCY : unrefined (raw) reference SVs called with Sniffles, filtered and normalized
	truth <- read_filter_vcf(vcf.file = paste0("../", line_names[i], "_normalized_ids.vcf"),
				 keep_hets = FALSE, max_N_del = 0, 
				 # DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
				 refseq = "../../../refgenome/Gmax_508_v4.0_mit_chlp.fasta",
				 keep_homref = FALSE, chrom_pattern = "^Gm[0-9]{2}$", 
				 remove_N_ins = TRUE, keep.ins.seq = TRUE, 
				 keep.ref.seq = TRUE, sample.name = "",
				 qual.field = "QUAL", check.inv = FALSE, 
				 keep.ids = FALSE, nocalls = FALSE, 
				 out.fmt = "gr", min.sv.size = 50)

	# Looping over the methods
	for(j in c("paragraph", "bayestyper", "vg")) {

		# Keeping the user informed
		message("Processing ", i, " --- ", j)

		# Setting the name of the input vcf depending on the pipeline
		# DEPENDENCY : raw SVs called by BayesTyper, Paragraph and vg
		if(j == "bayestyper") {
			input_vcf <- "../bayestyper/bayestyper_split.vcf"
		} else if (j == "paragraph") {
			input_vcf <- paste0("../paragraph/", i, "_results/genotypes.vcf.gz")
		} else if (j == "vg") {
			input_vcf <- paste0("../vg/", i, "/", i, "_calls.vcf.gz")
		}

		# Reading the vcf file for the Illumina variants
		calls <- read_filter_vcf(vcf.file = input_vcf, 
					 keep_hets = FALSE, max_N_del = 0, 
					 keep_homref = FALSE, chrom_pattern = "^Gm[0-9]{2}$", 
					 remove_N_ins = TRUE, keep.ins.seq = TRUE, 
					 keep.ref.seq = TRUE, sample.name = i,
					 qual.field = ifelse(j == "paragraph", "DP", "GQ"), 
					 other.field = "ACO", check.inv = FALSE, 
					 keep.ids = FALSE, nocalls = FALSE,
					 out.fmt = "gr", min.sv.size = 50)

		# Calling svevalOl on that file
		sveval_output  <- 
			svevalOl(calls.gr = calls,
				 truth.gr = truth,
				 ins.seq.comp = TRUE,
				 min.size = 50,
				 qual.field = ifelse(j == "paragraph", "DP", "GQ"),
				 sample.name = i,
				 qual.ths = NULL,
				 qual.quantiles = seq(0, 1, 0.05),
				 check.inv = FALSE,
				 geno.eval = FALSE,
				 stitch.hets = FALSE,
				 merge.hets = FALSE,
				 nb.cores = 1,
				 log.level = "INFO",
				 output.all = TRUE)

		# Processing the results with extract_rates
	        ij_rates  <- extract_rates(sveval_output, cultivar = i, pipeline = j)	
		rates_list[[paste0(i, "_", j)]] <- ij_rates
	}
}

# Merging the outputs of all analyses together
deletions  <- do.call("rbind", lapply(rates_list, function(x) `[[`(x, "DEL")))
insertions  <- do.call("rbind", lapply(rates_list, function(x) `[[`(x, "INS")))

# Creating a merged dataset and saving it to file
sveval_nogeno_rates <- list(DEL = deletions, INS = insertions)

# OUTPUT : breakpoint_refinement_analysis/raw_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData
save(sveval_nogeno_rates, file = "nogeno_RData/sveval_nogeno_rates.RData")

