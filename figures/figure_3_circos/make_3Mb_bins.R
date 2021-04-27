#!/usr/bin/Rscript

# This script simply creates a GRanges object that represents evenly-spaced
# 3-Mb bins along the reference chromosomes of the soybean v4 genome.

# This will make it easier to bin data and generate input files for plotting
# with Circos

# Loading the required libraries
library(Rsamtools)
library(GenomicRanges)

# Reading the fasta index of the reference genome and pruning the unanchored scaffolds
refgenome <- scanFaIndex("/home/malem420/refgenome/Gmax_v4/Gmax_508_v4.0_mit_chlp.fasta")
chrs <- paste0("Gm", sprintf("%0.2d", 1:20))
seqlevels(refgenome, pruning.mod = "coarse") <- chrs
gmax4_3Mb_bins <- tileGenome(seqinfo(refgenome), tilewidth = 3000000, cut.last.tile.in.chrom = TRUE)

# Saving to an RData object
save(gmax4_3Mb_bins, file = "gmax4_3Mb_bins.RData")

