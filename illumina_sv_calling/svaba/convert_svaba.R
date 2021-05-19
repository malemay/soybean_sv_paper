#!/usr/bin/Rscript

# Loading the GenomicRanges library
library(GenomicRanges)

# Sourcing the functions used to process svaba output
# DEPENDENCY : scripts/svaba_process.R
source("../../scripts/svaba_process.R")

# Creating a list of IDs over which to iterate
# DEPENDENCY : utilities/all_lines.txt
ids <- scan("../../utilities/all_lines.txt", what = character(), sep = "\n", quiet = TRUE)

# Looping over all individuals to classify and convert the alleles
for(i in ids) {

  cat("Processing ", i, "\n")

  # DEPENDENCY : SvABA SV vcf files
  svaba_classifier(vcf_file = paste0(i, ".svaba.sv.vcf"), 
                   output_file = paste0(i, ".classified.vcf"))

  # DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
  svaba_converter(vcf_file = paste0(i, ".classified.vcf"), 
                  output_file = paste0(i, ".converted.vcf"),
		  refgenome = "../../refgenome/Gmax_508_v4.0_mit_chlp.fasta",
		  contig_pattern = "^Gm[0-9]{2}$",
                  max_span = 500000)
}

