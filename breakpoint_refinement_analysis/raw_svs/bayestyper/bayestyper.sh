#!/bin/bash

# Getting the path to the executables from the command line
bt=$1
btools=$2
bcftools=$3

# Creating a variable for the reference genome
# DEPENDENCY : refgenome/Gmax_508_v4.0.fa
refgenome=../../../refgenome/Gmax_508_v4.0.fa

# Running BayesTyperTools combine on the input dataset to make it suitable for genotyping with bayestyper
# DEPENDENCY : breakpoint_refinement_analysis/raw_svs/svmerged.clustered.vcf
$bt combine -o combined_variants -v all_variants:../svmerged.clustered.vcf 

# DEPENDENCY : breakpoint_refinement_analysis/raw_svs/bayestyper/samples.tsv
# DEPENDENCY : refgenome/bt_decoy_sequences.fasta
# DEPENDENCY : Bayestyper k-mers
# Generating the clusters prior to genotyping
$bt cluster -v combined_variants.vcf -s samples.tsv -g $refgenome -d ../../../refgenome/bt_decoy_sequences.fasta \
	--min-number-of-unit-variants 100000 -p 8 -o bt_cluster

# Genotyping the SVs
$bt genotype -v bt_cluster_unit_1/variant_clusters.bin -c bt_cluster_cluster_data -s samples.tsv -g $refgenome \
	-d ../../../refgenome/bt_decoy_sequences.fasta -o bt_cluster_unit_1/bt_genotypes -p 10 --min-genotype-posterior 0.95

# Re-filtering the BayesTyper file such that calls require a minimum GPP of 0
$btools filter -v bt_cluster_unit_1/bt_genotypes.vcf -o bayestyper_gpp0 --min-genotype-posterior 0 \
	--kmer-coverage-file bt_cluster_unit_1/bt_genotypes_genomic_parameters.txt

# Splitting multiallelic sites because sveval does not input them all if they are joined
$bcftools norm --multiallelics -any -Ov bayestyper_gpp0.vcf > bayestyper_split.vcf

