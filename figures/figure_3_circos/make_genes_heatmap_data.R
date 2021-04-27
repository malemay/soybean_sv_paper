#!/usr/bin/Rscript

# This script prepares a data file to use for displaying gene regions
# as a heatmap on top of the ideograms

# Loading the GenomicFeatures package in order to read the gene annotations
library(GenomicFeatures)
library(GenomicRanges)

gmax_genes <- makeTxDbFromGFF("/home/malem420/refgenome/Gmax_v4/Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene.gff3")

# Extracting the genes and formatting for output in the format required by Circos
genes <- genes(gmax_genes)
chrs <- paste0("Gm", sprintf("%.2d", 1:20))
seqlevels(genes, pruning.mode = "coarse") <-  chrs

# Reading in the bins over which to tally the number of genes
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

output_track(genes, gmax4_3Mb_bins, "genes_heatmap.txt")

