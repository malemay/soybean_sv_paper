#!/bin/bash

# Getting the path to executables from the command line
bcftools=$1
tabix=$2

# Sorting the vcf files
# DEPENDENCY : utilities/all_lines.txt
# DEPENDENCY : VCF files of Illumina SVs genotyped by Paragraph
cat ../../utilities/all_lines.txt | parallel -j8 "$bcftools sort -m 1000M -Oz {}_results/genotypes.vcf.gz > {}_results/sorted_genotypes.vcf.gz"

# Indexing the vcf files
ls CAD*/sorted_genotypes.vcf.gz | parallel -j8 "$tabix {}"

# Merging the genotype calls of the 102 samples genotyped with Paragraph using bcftools (about 1h30)
$bcftools merge -Ov $(ls CAD1???_results/sorted_genotypes.vcf.gz) > paragraph_merged.vcf

