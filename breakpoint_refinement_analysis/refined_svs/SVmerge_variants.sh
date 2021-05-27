#!/bin/bash

# I have used rather stringent parameters so as not to lose information
# from the set of variants; alleles must be sufficiently precise to
# enable genotyping from Illumina data, so we won't merge variants
# that are too different

# WARNING : some programs may need to be available in $PATH for SVmerge to work properly
# Loading modules that may be required
# module load samtools/1.8
# module load bedtools/2.26.0
# module load mummer/3.23
# The edlib-aligner is in the path (/home/malem420/.local/bin/edlib-aligner)

# Getting the path to the SVmerge executable from the command line
SVmerge=$1

# Running the command
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY : breakpoint_refinement_analysis/refined_svs/files.txt
# DEPENDENCY : normalized refined SVs called by Sniffles
$SVmerge -ref ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	-fof files.txt \
	-prefix svmerged_preliminary \
	-maxdist 15 \
	-reldist 0.1 \
	-relsizediff 0.1 \
	-relshift 0.1 \
	-seqspecific

