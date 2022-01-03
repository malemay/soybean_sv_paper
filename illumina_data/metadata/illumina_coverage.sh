#!/bin/bash

# Getting the samtools executable from the command line
samtools=$1

# Looping over all samples
# DEPENDENCY : utilities/all_lines.txt
for sample in $(cat ../../utilities/all_lines.txt)
do
	# DEPENDENCY : illumina_data/ILLUMINA_ALIGNMENT
	bam_file=../aligned_reads/${sample}/${sample}_all.sort.bam
	$samtools coverage $bam_file > coverage/${sample}_coverage.txt
done

