#!/bin/bash

# First we will prepare the input Oxford Nanopore and Illumina SV vcfs
# so that they are given proper IDs and INFO values (SVTYPE + ACO) to
# be used as input to SVmerge

# Getting the path to the executables from the command line
bcftools=$1
SVmerge=$2

# WARNING : Some of the programs needed by SVmerge may need to be added to $PATH
# Loading modules that may be required
# module load samtools/1.8
# module load bedtools/2.26.0
# module load mummer/3.23
# The edlib-aligner is in the path (/home/malem420/.local/bin/edlib-aligner)

# First processing the Illumina input file
# DEPENDENCY : sv_genotyping/illumina_svs/svmerged.clustered.vcf
illumina_input=../illumina_svs/svmerged.clustered.vcf

# We keep only the ACO, SVTYPE and ClusterIDs INFO ; CLusterIDs is renamed to initial IDs because SVmerge will overwrite it
# We also change the ID to illumina_{row number} to have unique IDs for merging with the Oxford Nanopore dataset
# It will also make it easier to identify which variants come from the Illumina dataset, and which come from the Nanopore dataset
$bcftools annotate -x "^INFO/ACO,INFO/SVTYPE,INFO/ClusterIDs" -Ov $illumina_input | sed 's/ClusterIDs/InitialIDs/g' | \
	awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {$3 = "illumina_" NR; print}' > illumina_svs.vcf

# Then processing the Oxford Nanopore variants
# DEPENDENCY : sv_genotyping/nanopore_svs/svmerged_clustered_sorted.vcf
nanopore_input=../nanopore_svs/svmerged_clustered_sorted.vcf

# For these, we first begin by removing all annotations except SVTYPE
# Then over a single pass with awk we move the original ID to the INFO field
# Then we add the annotation header line with bcftools annotate
# DEPENDENCY : sv_genotyping/combined_svs/header_lines.txt
$bcftools annotate -x "^INFO/SVTYPE" -Ov $nanopore_input | \
	awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {$0=$0 ";ACO=sniffles;InitialIDs=" $3; $3="nanopore_" NR; print}' | \
	$bcftools annotate --header-lines header_lines.txt > nanopore_svs.vcf

# Using SVmerge for merging the Illumina and Oxford Nanopore SVs used as input to Paragraph
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
$SVmerge -ref ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	-fof svmerge_files.txt \
	-prefix svmerged_preliminary \
	-maxdist 15 \
	-reldist 0.2 \
	-relsizediff 0.1 \
	-relshift 0.1 

