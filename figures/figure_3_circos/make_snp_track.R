#!/usr/bin/Rscript

# This file prepares the circos data file for plotting SNP density as 
# a heatmap right inside of the gene density heatmap

# Loading the required libraries
library(GenomicRanges)
library(VariantAnnotation)

# Reading the VCF, limiting the input to that needed for overlap checking
# DEPENDENCY: structure_analysis/platypus_filtered_snps.vcf
vcf <- readVcf("../../structure_analysis/platypus_filtered_snps.vcf",
	       param = ScanVcfParam(info = NA, geno = NA))
vcf <- rowRanges(vcf)

chrs <- paste0("Gm", sprintf("%.2d", 1:20))
seqlevels(vcf, pruning.mode = "coarse") <- chrs

# Reading the 3-Mb bins from the .RData file
# DEPENDENCY: figures/figure_3_circos/gmax4_3Mb_bins.RData
load("gmax4_3Mb_bins.RData")

# Let us write a function that takes a GRanges object as input, some bins to count over
# and the name of a file to which the Circos data track will be output
output_track <- function(x, bins, output_file) {
	bins$count <- countOverlaps(bins, x)
	bins <- as.data.frame(bins)
	bins <- bins[, c("seqnames", "start", "end", "count")]
	write.table(bins, file = output_file, col.names = FALSE, row.names = FALSE,
		    quote = FALSE, sep = "\t")
}

# OUTPUT: figures/figure_3_circos/snp_density.txt
output_track(vcf, gmax4_3Mb_bins, "snp_density.txt")

