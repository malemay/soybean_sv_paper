#!/bin/bash

# Getting the path to the BayesTyper executable from the command line
bt=$1

# Creating a variable for the reference genome
# DEPENDENCY : refgenome/Gmax_508_v4.0.fa
refgenome=../../../refgenome/Gmax_508_v4.0.fa

# Running BayesTyperTools combine on the input dataset to make it suitable for genotyping with bayestyper
# DEPENDENCY : breakpoint_refinement_analysis/refined_svs/svmerged_clustered_sorted.vcf
$bt combine -o combined_variants -v all_variants:../svmerged_clustered_sorted.vcf

# DEPENDENCY : breakpoint_refinement_analysis/refined_svs/bayestyper/samples.tsv
# DEPENDENCY : refgenome/bt_decoy_sequences.fasta
# DEPENDENCY : Bayestyper k-mers
# Generating the clusters prior to genotyping
$bt cluster -v combined_variants.vcf -s samples.tsv -g $refgenome -d ../../../refgenome/bt_decoy_sequences.fasta \
	--min-number-of-unit-variants 100000 -p 8 -o bt_cluster

# Genotyping the SVs
$bt genotype -v bt_cluster_unit_1/variant_clusters.bin -c bt_cluster_cluster_data -s samples.tsv -g $refgenome \
	-d ../../../refgenome/bt_decoy_sequences.fasta -o bt_cluster_unit_1/bt_genotypes -p 10 --min-genotype-posterior 0.95

