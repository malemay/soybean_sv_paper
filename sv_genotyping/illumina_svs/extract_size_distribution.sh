#!/bin/bash

# Takes the vcf file merged with SVmerge and outputs a table with the SV types
# and sizes discovered by each tool. Will be used for plotting in figure S7
# of the project manuscript

# Creating a variable for the bcftools executable
bcftools=$1

# First writing a (tab-separated) header with the following columns
# program svtype size
printf "program\tsvtype\tsize\n" > size_distribution.tsv

# Then parsing only the input we need with bcftools and processing with awk
# DEPENDENCY : sv_genotyping/illumina_svs/svmerged.clustered.vcf
$bcftools query --format "%REF\t%ALT\t%INFO/SVTYPE\t%INFO/ClusterIDs\n" svmerged.clustered.vcf |
	awk '
	# This first block determines the size of the sv
	{
		if ($3 == "INV") 
			size = length($1);
		else 
			size = length($2)- length($1);
	}

	# Then we print a line for each matching tool
	/asmvar/ {printf "asmvar\t" $3 "\t" size "\n"}
	/smoove/ {printf "smoove\t" $3 "\t" size "\n"}
	/manta/ {printf "manta\t" $3 "\t" size "\n"}
	/svaba/ {printf "svaba\t" $3 "\t" size "\n"}
	' >> size_distribution.tsv

# OUTPUT sv_genotyping/illumina_svs/size_distribution.tsv

