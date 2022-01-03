#!/bin/bash

# Getting the paths to the executables from the command line
Rscript=$1
plink=$2
phylip=$3

# Using plink to generate distance matrices from SNPs
# DEPENDENCY : structure_analysis/snp_pca/SNP_PCA (because of the files in PLINK format created by the PCA script)
$plink --bfile ../snp_pca/platypus_filtered_snps --distance square 1-ibs --out snp_distance

# Converting the distance matrix to a format that PHYLIP can use
$Rscript -e "source('convert_to_phylip.R') ; to_phylip('snp_distance.mdist')"

# Generating the input commands for PHYLIP and running the neighbour-joining tree analysis
echo "snp_distance.mdist.phyl" > snp_phylip_params.txt
echo "Y" >> snp_phylip_params.txt

$phylip < snp_phylip_params.txt > snp_phylip.log

# Changing the name of the output tree
mv outtree snp_outtree
rm outfile

