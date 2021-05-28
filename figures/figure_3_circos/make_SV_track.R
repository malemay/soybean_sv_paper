#!/usr/bin/Rscript

# This script prepares the data for the stacked histogram of insertions and deletions
# that will be plotted as a track between the SNP density track and the transposable
# element histogram tracks

# Loading the required libraries
library(VariantAnnotation)
library(GenomicRanges)

# Reading the filtered SV vcf file
# DEPENDENCY: sv_genotyping/combined_svs/combined_paragraph_filtered.vcf
vcf <- readVcf("../../sv_genotyping/combined_svs/combined_paragraph_filtered.vcf",
	       param = ScanVcfParam(info = "SVTYPE", geno = NA))

deletions  <- vcf[info(vcf)$SVTYPE == "DEL"]
insertions <- vcf[info(vcf)$SVTYPE == "INS"]

# Loading the 3-Mb bins from the RData file
# DEPENDENCY: figures/figure_3_circos/gmax4_3Mb_bins.RData
load("gmax4_3Mb_bins.RData")

# Adding the overlaps with the deletions
gmax4_3Mb_bins$deletions <- countOverlaps(gmax4_3Mb_bins, deletions)
gmax4_3Mb_bins$insertions <- countOverlaps(gmax4_3Mb_bins, insertions)

# Coercing the GRanges object into a data.frame for manipulation and writing to disk
sv_counts <- as.data.frame(gmax4_3Mb_bins)
sv_counts$seqnames <- as.character(sv_counts$seqnames)
sv_counts$value <- paste0(sv_counts$deletions, ",", sv_counts$insertions)
sv_counts <- sv_counts[, c("seqnames", "start", "end", "value")]

# Writing to file for plotting with Circos
# OUTPUT: figures/figure_3_circos/sv_counts.txt
write.table(sv_counts, file = "sv_counts.txt", col.names = FALSE, row.names = FALSE,
	    sep = "\t", quote = FALSE)

