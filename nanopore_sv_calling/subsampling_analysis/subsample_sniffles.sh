#!/bin/bash

# Getting the path to the executables from the command line
samtools=$1
sniffles=$2

# Looping over the 80 rows in the subsample_parameters.txt file
for i in $(seq 1 80)
do
	params=$(head -n $i subsample_parameters.txt | tail -n 1)
	sample=$(echo $params | cut -d " " -f1)
	frac=$(echo $params | cut -d " " -f2)
	seed=$(echo $params | cut -d " " -f3)

	# Creating the output directory (if not already created) and assigning the input bam file to a variable
	mkdir -p $sample
	# DEPENDENCY : nanopore_data/NANOPORE_ALIGNMENT (Oxford Nanopore .bam files)
	bam_file=../../nanopore_data/${sample}_porechopped_aligned.sort.bam

	# Subsampling and indexing the file
	$samtools view -bh --subsample ${frac} --subsample-seed $seed $bam_file > ${sample}/${sample}_frac${frac}_seed${seed}.bam
	$samtools index ${sample}/${sample}_frac${frac}_seed${seed}.bam

	# Getting the bam input file name for this run
	sniffles_input=${sample}/${sample}_frac${frac}_seed${seed}.bam
	sniffles_output=$(echo "$sniffles_input" | sed "s:\.bam$:.vcf:")

	$sniffles -t 1 --min_support 3 --min_seq_size 1000 --min_homo_af 0.7 -m $sniffles_input -v $sniffles_output
done

