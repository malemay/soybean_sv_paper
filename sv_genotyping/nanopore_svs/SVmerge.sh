#!/bin/bash

# Using SVmerge for merging the SVs from all realigned vcf files

# I have used rather stringent parameters so as not to lose information
# from the set of variants; alleles must be sufficiently precise to
# enable genotyping from Illumina data, so we won't merge variants
# that are too different

# WARNING : SVmerge may require some programs to be added to the PATH
# Loading modules that may be required
# module load samtools/1.8
# module load bedtools/2.26.0
# module load mummer/3.23
# The edlib-aligner is in the path (/home/malem420/.local/bin/edlib-aligner)

# Creating a variable for the SVmerge executable
SVmerge=$1

# Running the command
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY : normalized SV files from Oxford Nanopore
$SVmerge -ref ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	-fof svmerge_files.txt \
	-prefix svmerged_preliminary \
	-maxdist 15 \
	-reldist 0.2 \
	-relsizediff 0.1 \
	-relshift 0.1 \
	-seqspecific

