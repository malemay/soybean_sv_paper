#!/bin/bash

# Getting the path to the executables from the command line
vcftools=$1
plink=$2

# Using vcftools and plink to convert the vcf to a format suitable for input to fastStructure
# DEPENDENCY : structure_analysis/platypus_filtered_snps.vcf
$vcftools --vcf ../platypus_filtered_snps.vcf --out platypus_filtered_snps --plink
$plink --noweb --ped platypus_filtered_snps.ped --map platypus_filtered_snps.map --make-bed --out platypus_filtered_snps

# Using plink to generate pca
plink --bfile platypus_filtered_snps --pca --out snp_pca

