#!/bin/bash

# This code uses bcftools to:
# 1- Filter the population-scale SV vcf by setting to missing the genotype calls with DP < 2
# 2- Annotate the vcf file with various measures, including the number of alternate alleles
#    in homozygous ALT allele genotype calls (INFO/AC_Hom) which we will filter on as a quality
#    measure for evaluation with svevalOl

# Getting the path to the bcftools executable from the command line
bcftools=$1

# Exporting the variable BCFTOOLS_PLUGINS so that bcftools can find the fill-tags plugin
export BCFTOOLS_PLUGINS=$2

# Path of the input file
# DEPENDENCY : sv_genotyping/illumina_svs/illumina_paragraph_merged.vcf
# Applying the filter and only then computing the tags to add
$bcftools filter --exclude "FORMAT/DP < 2" --set-GTs . -Ou illumina_paragraph_merged.vcf | \
	$bcftools plugin fill-tags -Ov - > paragraph_svs_minDP2_annotated.vcf

