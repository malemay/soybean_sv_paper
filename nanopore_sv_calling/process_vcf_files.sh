#!/bin/bash

# DEPENDENCY : parallel
# DEPENDENCY : bcftools
# DEPENDENCY : R
bcftools=$1
r_exec=$2

# Creating a variable for the reference genome
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Adding the metainfo to all the samples
# DEPENDENCY : nanopore_sv_calling/add_metainfo_all.R
$r_exec add_metainfo_all.R

# Fixing the the reference sequences based on the SV positions
# DEPENDENCY : nanopore_sv_calling/fix_vcfs.R
$r_exec fix_vcfs.R

# Sort vcf entries
# DEPENDENCY : utilities/line_ids.txt
tail -n+2 ../utilities/line_ids.txt | cut -f 1 | parallel -j1 "$bcftools sort -m 500.0M -Ov {}_fixed.vcf > {}_sorted.vcf"

# Normalize the variants and remove duplicates
tail -n+2 ../utilities/line_ids.txt | cut -f 1 | parallel -j1 "$bcftools norm --rm-dup exact -f $refgenome {}_sorted.vcf > {}_normalized.vcf"

# Add IDs such that SVmerge does not complain
tail -n+2 ../utilities/line_ids.txt | cut -f 1 | parallel -j1 "$bcftools annotate --set-id '{}\_%CHROM\_%POS\_%INFO/SVLEN' -Ov {}_normalized.vcf > {}_normalized_ids.vcf"

