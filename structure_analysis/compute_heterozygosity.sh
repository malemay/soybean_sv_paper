#!/bin/bash

# Getting the path to vcftools from the command line
vcftools=$1

# DEPENDENCY : structure_analysis/platypus_snps.vcf
$vcftools --vcf platypus_filtered_snps.vcf --het --stdout > snp_heterozygosity_rates.txt

