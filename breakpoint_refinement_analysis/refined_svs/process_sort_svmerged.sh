#!/bin/bash

# Getting the paths to the executables from the command line
rscript=$1
bcftools=$2

# Merging the variants by systematically favouring those that have been refined
# DEPENDENCY : breakpoint_refinement_analysis/refined_svs/select_svs.R
$rscript select_svs.R

# Sort vcf entries
$bcftools sort -m 500.0M -Ov svmerged.clustered.vcf > svmerged_clustered_sorted.vcf

