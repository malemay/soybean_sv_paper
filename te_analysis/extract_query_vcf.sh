#!/bin/bash

# Variants smaller than 100 nucleotides are filtered out prior to TE analysis
# Only variants in this file will be tested for matches to the transposable
# element database

# Using awk to only extract variants larger than 100 nucleotides
# DEPENDENCY : sv_genotyping/combined_svs/combined_paragraph_filtered.vcf
awk '
BEGIN {OFS = "\t"}

/^#/ {print $0}

!/^#/ {
	if((length($4) - length($5)) >= 100) {
		print $0
	} else if ((length($5) - length($4)) >= 100) {
		print $0
        }
}
' ../sv_genotyping/combined_svs/combined_paragraph_filtered.vcf > query_all.vcf

