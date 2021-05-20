#!/usr/bin/Rscript

# Using a custom R function to re-select the SVs clustered by SVmerge by preferentially
# choosing realigned variants

# DEPENDENCY : scripts/merge_realigned_variants.R
source("../../scripts/merge_realigned_variants.R")

# Generating a character vector of input vcfs
# DEPENDENCY : sv_genotyping/nanopore_svs/svmerge_files.txt
input_vcfs <- scan("svmerge_files.txt", what = character(), sep = "\n", quiet = TRUE)

# Lauching the command
merge_aligned_variants("svmerged_preliminary.clustered.vcf",
		       input_vcfs,
		       "svmerged.clustered.vcf")

