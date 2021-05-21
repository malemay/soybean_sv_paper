#!/usr/bin/Rscript

# Loading the required packages
library(rtracklayer)
library(GenomicRanges)
library(Rsamtools)

# Reading the repeats file so we can output a bed file that includes only non-repeat regions
# The purpose of this will be to use that bed file to analyze SV recall and precision
# in high-confidence (non-repeat) regions
# DEPENDENCY : refgenome/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3
repeat_file <- "../Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3"
repeats <- import(repeat_file)

# Defining a vector of standard chromosome names
chrs <- paste0("Gm", sprintf('%0.2d', 1:20))

# We use this vector to drop the ranges on unassembled scaffolds and redefine the levels
seqlevels(repeats, pruning.mode = "coarse") <- chrs

# Now let us reduce this range, ignoring the strand information
repeats <- reduce(repeats, ignore.strand = TRUE)

# Now we want to get all the regions that are NOT repeats
# For this, we need # the information from the reference genome
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome <- scanFaIndex("../Gmax_508_v4.0_mit_chlp.fasta")

# Let us keep only the standard chromosomes
seqlevels(refgenome, pruning.mode = "coarse") <- chrs

# We can use the setdiff function to get the ranges that are in the reference
# genome but not in the repeats, effectively getting the non-repeated regions
non_repeats <- setdiff(refgenome, repeats)

# We can now save this to a bed file
# OUTPUT : refgenome/repeat_regions/non_repeated_regions.bed
export(non_repeats, "non_repeated_regions.bed")

