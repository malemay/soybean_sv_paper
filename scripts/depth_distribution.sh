#!/bin/bash

# Getting the samtools executable from the command line
samtools=$3

# Giving more meaningful names to positional arguments
input_file=$1
output_file=$2

# Creating a string of chromsome names to extract the chrosomes Gm01 through Gm20 from the bam file
chrom="Gm01 Gm02 Gm03 Gm04 Gm05 Gm06 Gm07 Gm08 Gm09 Gm10 Gm11 Gm12 Gm13 Gm14 Gm15 Gm16 Gm17 Gm18 Gm19 Gm20"

# Computing the sequencing depths at all position and processing them through awk to get the distribution
$samtools view -bh $input_file $chrom | $samtools depth -a - | \
	awk '{depth_counts[$3] += 1} END {for (i in depth_counts) {print i, depth_counts[i]} }' > $output_file

