#!/bin/bash

# Getting the path to the executables from the command line
Rscript=$1
plink=$2
phylip=$3
vcftools=$4
bcftools=$5

# We use the SV file as the input file
# DEPENDENCY : sv_genotyping/combined_svs/combined_paragraph_filtered.vcf
input=../../../../sv_genotyping/combined_svs/combined_paragraph_filtered.vcf

# Performing the analysis for 100 iterations
for iter in $(seq 1 100)
do
	# Creating a directory for this iteration and moving to it
	mkdir -p bootstrap/sv_bootstrap_iter${iter}
	cd bootstrap/sv_bootstrap_iter${iter}

	# Resampling with replacement
	$bcftools view -h $input > sv_bootstrap_iter${iter}.vcf
	nrecords=$(grep -cv "^#" $input)
	$bcftools view -H $input | shuf -rn $nrecords >> sv_bootstrap_iter${iter}.vcf

	# Using vcftools and plink to convert the vcf to a format suitable for input to fastStructure
	$vcftools --vcf sv_bootstrap_iter${iter}.vcf --out sv_bootstrap_iter${iter} --plink
	$plink --noweb --ped sv_bootstrap_iter${iter}.ped --map sv_bootstrap_iter${iter}.map --make-bed --out sv_bootstrap_iter${iter}
	rm sv_bootstrap_iter${iter}.vcf

	# Computing the distance matrix with plink
	$plink --bfile sv_bootstrap_iter${iter} --distance square 1-ibs --out sv_bootstrap_iter${iter}_distance

	# Converting to the format required by PHYLIP using an R script
	# DEPENDENCY : structure_analysis/trees/convert_to_phylip.R
	$Rscript -e 'source("../../convert_to_phylip.R"); to_phylip("sv_bootstrap_iter'${iter}'_distance.mdist")'

	# Prepraring the input file for PHYLIP
	echo "sv_bootstrap_iter${iter}_distance.mdist.phyl" > phylip_params.txt
	echo "Y" >> phylip_params.txt

	# Running the neighbor program of PHYLIP
	$phylip < phylip_params.txt > phylip.log

	# Going back up two directories
	cd ../..
done

# Updating the timestamp of the SV_BOOTSTRAP file so the Makefile works properly
touch SV_BOOTSTRAP

