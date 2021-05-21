#!/bin/bash

# This code uses bcftools to:
# 1- Filter the population-scale SV vcf by removing the variants with a AC_Hom < 4
# 2- Add an INFO tag with the number of independent calling programs that have called
#    that variant. We use awk to count the number of such programs for each line based
#    on the "ClusterIDs" INFO tag

# Creating a variable for the bcftools executable
bcftools=$1

# Applying the filter and only then computing the tags to add
# DEPENDENCY : sv_genotyping/illumina_svs/sveval_benchmarks/paragraph_svs_minDP2_annotated.vcf
# DEPENDENCY : sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_header_line.txt
$bcftools view --exclude "INFO/AC_Hom < 4" -Ou paragraph_svs_minDP2_annotated.vcf | \
       	$bcftools annotate --header-lines ncallers_header_line.txt -Ov - |
	awk '
	BEGIN {OFS = "\t"}
	/^#/ {print $0}
	{x = 0}
	!/^#/ && $8 ~ "asmvar" {x+=1}
	!/^#/ && $8 ~ "manta"  {x+=1}
	!/^#/ && $8 ~ "smoove" {x+=1}
	!/^#/ && $8 ~ "svaba"  {x+=1}
	!/^#/ {$8 = $8 ";NCALLERS=" x; print $0}
' > paragraph_svs_minAC4_ncallers.vcf

