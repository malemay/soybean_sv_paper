#!/usr/bin/Rscript

# File initially created on Wednesday, January 20, 2021

# Benchmarking the structural variants called by Illumina from the variants discovered by Illumina and Oxford Nanopore combined,
#  but with Illumina variants being favoured at the SVmerge step, on the fourth version of the soybean reference genome assembly with sveval. 

# The only difference with the code in benchmark.R is that here we restrict our analysis to the non-repeated regions found in the imported bed file

# Loading the required libraries
library(rtracklayer)
library(sveval) # (modified by me to allow for argument output.all)

# Loading the bed file that contains the non-repeat regions
# DEPENDENCY : refgenome/repeat_regions/non_repeated_regions.bed
non_repeats <- import("../../../refgenome/repeat_regions/non_repeated_regions.bed")

# Sourcing the extract_rates.R file which contains the fonction to gather precision/sensitivity rates
# DEPENDENCY : scripts/extract_rates.R
source("../../../scripts/extract_rates.R")

# And another function that reads a vcf file into a GRanges object and filters it according to various criteria
# DEPENDENCY : scripts/read_filter_vcf.R
source("../../../scripts/read_filter_vcf.R")

# Reading the data.frame giving the correspondence between line names and line ids
# DEPENDENCY : utilities/line_ids.txt
line_ids <- read.table("../../../utilities/line_ids.txt", header = TRUE, stringsAsFactors = FALSE)

# This vector will hold the correspondence between Illumina CAD IDs and line names
line_names <- line_ids$name
names(line_names) <- line_ids$id

# Initializing a list to store the results of each iteration
rates_list <- list()

# Looping over all 17 samples
for(i in names(line_names)) {

	# Reading the vcf file for the Nanopore variants
	# DEPENDENCY : normalized SVs discovered by Sniffles from the Oxford Nanopore data
	truth <- read_filter_vcf(vcf.file = paste0("../../../nanopore_sv_calling/", line_names[i], "_normalized_ids.vcf"),
				 keep_hets = FALSE, max_N_del = 0, 
				 # DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
				 refseq = "../../../refgenome/Gmax_508_v4.0_mit_chlp.fasta",
				 keep_homref = FALSE, chrom_pattern = "^Gm[0-9]{2}$", 
				 remove_N_ins = TRUE, keep.ins.seq = TRUE, 
				 keep.ref.seq = TRUE, sample.name = "",
				 qual.field = "QUAL", check.inv = FALSE, 
				 keep.ids = FALSE, nocalls = FALSE, 
				 out.fmt = "gr", min.sv.size = 50)

	# Setting j to paragraph because it's the only tool we are using here
	j <- "paragraph"

	# Keeping the user informed
	message("Processing ", i, " --- ", j)

	# Setting the name of the input vcf depending on the pipeline
	# DEPENDENCY : Combined Illumina/Oxford Nanopore SVs genotyped with Paragraph
	input_vcf  <- paste0("../", i, "_results/genotypes.vcf.gz")

	# Reading the vcf file for the Illumina variants
	calls <- read_filter_vcf(vcf.file = input_vcf, 
				 keep_hets = FALSE, max_N_del = 0, 
				 keep_homref = FALSE, chrom_pattern = "^Gm[0-9]{2}$", 
				 remove_N_ins = TRUE, keep.ins.seq = TRUE, 
				 keep.ref.seq = TRUE, sample.name = i,
				 qual.field = "DP", check.inv = FALSE, 
				 other.field = "ClusterIDs",
				 keep.ids = FALSE, nocalls = FALSE,
				 out.fmt = "gr", min.sv.size = 50)

	# Calling svevalOl on that file
	sveval_output  <- svevalOl(calls.gr = calls,
				   truth.gr = truth,
				   ins.seq.comp = TRUE,
				   min.size = 50,
				   qual.field = "DP",
				   sample.name = i,
				   qual.ths = NULL,
				   qual.quantiles = seq(0, 1, 0.05),
				   check.inv = FALSE,
				   geno.eval = FALSE,
				   stitch.hets = FALSE,
				   merge.hets = FALSE,
				   nb.cores = 1,
				   # the bed file comes in here
				   bed.regions = non_repeats,
				   bed.regions.ol = 0.8,
				   log.level = "INFO",
				   output.all = TRUE)

	# Processing the results with extract_rates
	ij_rates  <- extract_rates(sveval_output, c("DEL", "INS", "INV", "DUP"),
				   cultivar = i, pipeline = j)	

	rates_list[[paste0(i, "_", j)]] <- ij_rates
}

# Merging the outputs of all analyses together
deletions  <- do.call("rbind", lapply(rates_list, function(x) `[[`(x, "DEL")))
insertions  <- do.call("rbind", lapply(rates_list, function(x) `[[`(x, "INS")))
inversions  <- do.call("rbind", lapply(rates_list, function(x) `[[`(x, "INV")))
duplications  <- do.call("rbind", lapply(rates_list, function(x) `[[`(x, "DUP")))

# Creating a merged dataset and saving it to file
sveval_norepeat_rates <- list(DEL = deletions, INS = insertions,
			    INV = inversions, DUP = duplications)

save(sveval_norepeat_rates, file = "norepeat_RData/sveval_norepeat_rates.RData")

