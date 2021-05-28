#!/bin/bash

# WARNING : python 2.7 may need to be added to $PATH for fastStructure to work
# Loading the required modules
# module load python/2.7

# Getting the path to the executables from the command line
faststructure=$1
vcftools=$2
plink=$3

# Using vcftools and plink to convert the vcf to a format suitable for input to fastStructure
# DEPENDENCY : structure_analysis/platypus_filtered_snps.vcf
$vcftools --vcf platypus_filtered_snps.vcf --out platypus_filtered_snps --plink
$plink --noweb --ped platypus_filtered_snps.ped --map platypus_filtered_snps.map --make-bed --out platypus_filtered_snps

# Running fastStructure with K=5 as in Torkamaneh et al. 2018
$faststructure -K 5 --input=platypus_filtered_snps --output=structure --full 

