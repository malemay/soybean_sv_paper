#!/usr/bin/Rscript

# This script simply creates a GRanges object that represents evenly-spaced
# 3-Mb bins along the reference chromosomes of the soybean v4 genome.

# This will make it easier to bin data and generate input files for plotting
# with Circos

# Loading the required libraries
library(Rsamtools)
library(GenomicRanges)

# Reading the fasta index of the reference genome and pruning the unanchored scaffolds
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome <- scanFaIndex("../../refgenome/Gmax_508_v4.0_mit_chlp.fasta")
chrs <- paste0("Gm", sprintf("%0.2d", 1:20))
seqlevels(refgenome, pruning.mod = "coarse") <- chrs
gmax4_3Mb_bins <- tileGenome(seqinfo(refgenome), tilewidth = 3000000, cut.last.tile.in.chrom = TRUE)

# Saving to an RData object
# OUTPUT: figures/figure_3_circos/gmax4_3Mb_bins.RData
save(gmax4_3Mb_bins, file = "gmax4_3Mb_bins.RData")

