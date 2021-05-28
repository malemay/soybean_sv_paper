#!/bin/bash

# This code uses bcftools to:
# 1- Filter the population-scale SV vcf by setting to missing the genotype calls with DP < 2
# 2- Annotate the vcf file with various measures, including the number of alternate alleles
#    in homozygous ALT allele genotype calls (INFO/AC_Hom) which we will filter on
# 3- Filter out the variants for which AC_Hom < 4

# Duplications and inversions are also filtered out using grep

# Getting the path to executables/plugins from the command line
bcftools=$1
export BCFTOOLS_PLUGINS=$2

# Applying the filter and only then computing the tags to add
# DEPENDENCY : sv_genotyping/combined_svs/combined_paragraph_merged.vcf
$bcftools filter --exclude "FORMAT/DP < 2" --set-GTs . -Ou combined_paragraph_merged.vcf | \
	$bcftools plugin fill-tags -Ou - | \
	$bcftools view --exclude "INFO/AC_Hom < 4" -Ov - | \
	grep -v "SVTYPE=DUP" | grep -v "SVTYPE=INV" > combined_paragraph_filtered.vcf

