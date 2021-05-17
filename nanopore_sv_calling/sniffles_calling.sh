#!/bin/bash

# Creating a variable for the sniffles executable
sniffles=$1
rscript=$2

# DEPENDENCY : nanopore_data/AC2001_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/ALTA_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/MAPLE_ISLE_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/MAPLE_PRESTO_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_09_35C_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_CARMAN_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_DRAYTON_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_EMBRO_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_LAKEVIEW_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_MADOC_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_OXFORD_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_PETREL_2_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_PRUDENCE_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_STRATFORD_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OT09-03_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/QS5091.50J_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/ROLAND_porechopped_aligned.sort.bam

# Using parallel to process all 17 lines with Sniffles
tail -n+2 ../utilities/line_ids.txt | cut -f1 | parallel -j10 "$sniffles -t 1 --min_support 3 --min_seq_size 1000 --min_homo_af 0.7 -m ../nanopore_data/{}_porechopped_aligned.sort.bam -v {}_sniffles_minreads3_hom70.vcf"

# Filtering the results using an R script
# DEPENDENCY : nanopore_sv_calling/filter_sniffles_vcfs.R
$rscript filter_sniffles_vcfs.R

