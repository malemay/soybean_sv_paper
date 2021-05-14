#!/usr/bin/Rscript

# Loading the function
# DEPENDENCY : scripts/fix_sniffles.R
source("../scripts/fix_sniffles.R")

# This step fixes the coordinates and alt/ref sequences of the filtered
# variants so that they can be processed by bcftools norm and later
# merged using either SVmerge or SURVIVOR

# Looping over the samples
# DEPENDENCY : utilities/line_ids.txt
for(i in read.table("../utilities/line_ids.txt", header = TRUE, stringsAsFactors = FALSE)[[1]]) {

	input_file  <- paste0(i, "_realigned_metainfo.vcf") 
	output_file <- paste0(i, "_fixed.vcf")

	# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
	fix_sniffles(input_file, output_file, refgenome = "../refgenome/Gmax_508_v4.0_mit_chlp.fasta")
}

