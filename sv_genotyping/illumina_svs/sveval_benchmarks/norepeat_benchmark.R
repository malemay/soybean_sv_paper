# File initially created on Wednesday, January 20, 2021

# Benchmarking the structural variants called by Illumina from the variants discovered by Illumina
#  on the fourth version of the soybean reference genome assembly with sveval.  The difference with 
#  the code in benchmark.R is that here we only include SVs located in non-repeated regions according 
#  to the following file

# Loading the required libraries
library(rtracklayer)
library(ggplot2)
library(sveval) # (modified by me to allow for argument output.all)

# DEPENDENCY : refgenome/repeat_regions/non_repeated_regions.bed
non_repeats <- import("../../../refgenome/repeat_regions/non_repeated_regions.bed")

# Sourcing the extract_rates.R file which contains the fonction to gather precision/sensitivity rates
# DEPENDENCY : scripts/extract_rates.R
source("../../../scripts/extract_rates.R")

# And another function that reads a vcf file into a GRanges object and filters
#  it according to various criteria
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
# DEPENDENCY : nanopore_sv_calling/AC2001_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/ALTA_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/MAPLE_ISLE_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/MAPLE_PRESTO_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_09_35C_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_CARMAN_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_DRAYTON_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_EMBRO_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_LAKEVIEW_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_MADOC_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_OXFORD_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_PETREL_2_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_PRUDENCE_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_STRATFORD_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/OT09-03_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/QS5091.50J_normalized_ids.vcf
# DEPENDENCY : nanopore_sv_calling/ROLAND_normalized_ids.vcf
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
	truth <- read_filter_vcf(vcf.file = paste0("../../../nanopore_sv_calling/", line_names[i], "_normalized_ids.vcf"),
				 keep_hets = FALSE, max_N_del = 0, 
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
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1077_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1018_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1092_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1010_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1089_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1064_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1056_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1052_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1015_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1065_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1049_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1022_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1002_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1074_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1096_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1070_results/genotypes.vcf.gz
# DEPENDENCY : sv_genotyping/illumina_svs/CAD1087_results/genotypes.vcf.gz
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
				   # here are the two input lines relative to the bed input
				   bed.regions = non_repeats,
				   bed.regions.ol = 0.8,
				   log.level = "INFO",
				   output.all = TRUE)

	# Processing the results with extract_rates
	ij_rates  <- extract_rates(sveval_output, c("DEL", "INS", "INV", "DUP"), cultivar = i, pipeline = j)	
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

# OUTPUT : sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData
save(sveval_norepeat_rates, file = "norepeat_RData/sveval_norepeat_rates.RData")

