#!/bin/bash

# File initially created on November 4, 2020

# The purpose of this code is to filter the vcf file obtained from the
# analysis on genome version 4 .

# The following filters will be applied:
# - Only sites on chromosomes Gm01 ... Gm20 are kept
# - Only biallelic SNPs are kept
# - Only records with FILTER == PASS or .
# - Maximum missing data of 40%
# - Minimum minor allele frequency of 0.05
# - Both alleles have more homozygous calls than number of heterozygous calls
# - Maximum heterozygosity rate of 0.1

# Getting the path to the executables from the command line
bcftools=$1
vcftools=$2
rscript=$3

# Creating a variable for the reference genome
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Extracting only biallelic SNPs from this file and removing sites with FILTER != PASS or .
# DEPENDENCY : structure_analysis/platypus_snps.vcf
$bcftools view --apply-filters .,PASS --regions Gm01,Gm02,Gm03,Gm04,Gm05,Gm06,Gm07,Gm08,Gm09,Gm10,Gm11,Gm12,Gm13,Gm14,Gm15,Gm16,Gm17,Gm18,Gm19,Gm20 \
	platypus_snps.vcf -Ou | \
	$bcftools norm -f $refgenome --multiallelics -any -Ou - | \
       	$bcftools view -v snps -Ou - | \
	$bcftools norm --multiallelics +snps -Ou - | \
       	$bcftools view -m2 -M2 -Ov - > platypus_biallelic.vcf

# Getting only the SNPs with MAF > 0.05 and no more than 40% missing data
$vcftools --vcf platypus_biallelic.vcf --maf 0.05 --max-missing 0.6 --recode --stdout > platypus_maf_maxmissing.vcf

# Filtering for heterozygosity
# First we count the number of observations for each genotype
# DEPENDENCY : scripts/het_counts.awk
../scripts/het_counts.awk platypus_maf_maxmissing.vcf > hetrates.txt

# Then we use this file to compute a table of rates to use for filtering out
# DEPENDENCY : scripts/het_stats.awk
gawk '{print $5}' hetrates.txt | ../scripts/het_stats.awk 102 > wgs_hetstats.txt

# Calling a R script to print a list of sites passing the filter; will be used by vcftools
# DEPENDENCY : scripts/filter_hetsites.R
$rscript ../scripts/filter_hetsites.R > kept_sites.txt

# Calling vcftools to filter the vcf
$vcftools --vcf platypus_maf_maxmissing.vcf --positions kept_sites.txt --recode --stdout > platypus_filtered_snps.vcf

