#!/bin/bash

# Using SVmerge for merging the SVs used as input to Paragraph

# WARNING : This command may require some programs to be available in $PATH
# Loading modules that may be required
# module load samtools/1.8
# module load bedtools/2.26.0
# module load mummer/3.23
# The edlib-aligner is in the path (/home/malem420/.local/bin/edlib-aligner)

# Getting the path to the SVmerge executable from the command line
SVmerge=$1

# Running the command
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY : sv_genotyping/illumina_svs/svmerge_files.txt
$SVmerge -ref ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	-fof svmerge_files.txt \
	-prefix svmerged \
	-maxdist 15 \
	-reldist 0.2 \
	-relsizediff 0.1 \
	-relshift 0.1 

