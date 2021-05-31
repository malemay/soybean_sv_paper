#!/bin/bash

# The filtered dataset is further filtered by vcftools to keep only variants with
# a maximum of 60% missing data

# PLINK is then used to format the data for the PCA and to actually compute the PCA

# Getting the path to the executables from the command line
vcftools=$1
plink=$2

# Filtering for missing data
# DEPENDENCY : sv_genotyping/combined_svs/combined_paragraph_filtered.vcf
vcftools --vcf ../../sv_genotyping/combined_svs/combined_paragraph_filtered.vcf --max-missing 0.4 --recode --stdout > combined_filtered_maxmiss.vcf

# Changing all alleles to REF = A and ALT = T because the format modification won't deal properly with SVs (a few seconds)
# DEPENDENCY : scripts/recode_alleles.awk
../../scripts/recode_alleles.awk combined_filtered_maxmiss.vcf > sv_pca_input.vcf

# Using vcftools and plink to convert the vcf to a format suitable for input plink --pca (a few seconds)
vcftools --vcf sv_pca_input.vcf --out sv_pca_input --plink
plink --noweb --ped sv_pca_input.ped --map sv_pca_input.map --make-bed --out sv_pca_input

# Using plink to generate pca
plink --bfile sv_pca_input --pca --out sv_pca

