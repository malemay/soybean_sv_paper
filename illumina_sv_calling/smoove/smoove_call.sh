#!/bin/bash

# WARNING : This code might require some dependencies to be added to $PATH

# Making some required modules available
# module load samtools/1.10
# module load bcftools/1.8
# module load htslib/1.10.2

# Getting the paths for the smoove executable from the command line
smoove=$1
tabix=$2

#Creating a variable for the reference genome location
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# DEPENDENCY : the Illumina reads aligned with BWA
# Calling the genotypes in parallel with smoove call
cat ../../utilities/all_lines.txt | parallel -j 8 "$smoove call --outdir results-smoove/ --excludechroms mitochondrion,chloroplast --name {} --fasta $refgenome -p 1 --genotype ../../illumina_data/aligned_reads/{}/{}_all.sort.bam " 

# Merging the genotypes using smoove merge
$smoove merge --name merged -f $refgenome --outdir ./results-smoove results-smoove/*.genotyped.vcf.gz

# Calling the genotypes with smoove call
cat ../../utilities/all_lines.txt | parallel -j 24 "$smoove genotype -d -x -p 1 --name {}_joint --outdir results-genotyped/ --fasta $refgenome --vcf results-smoove/merged.sites.vcf.gz ../../illumina_data/aligned_reads/{}/{}_all.sort.bam"

# Pasting the results of all samples in a single vcf file
$smoove paste --name all_samples results-genotyped/*.vcf.gz

# Indexing the result with tabix
$tabix all_samples.smoove.square.vcf.gz

