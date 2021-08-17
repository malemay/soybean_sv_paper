#!/bin/bash

# Extracting metadata from the vcf files of each of the SV calling programs
#  in order to make table 1, which will summarize the number of variants of
#  each SV type and size class for each program and for the merged data set
#  This table is valid for the programs calling SVs from Illumina data.

# Looping over the files that contain all the SVs from every tool
# DEPENDENCY : illumina_sv_calling/asmvar/asmvar_filtering/asmvar_svs.vcf
# DEPENDENCY : illumina_sv_calling/manta/manta_svs.vcf
# DEPENDENCY : illumina_sv_calling/smoove/smoove_svs.vcf
# DEPENDENCY : illumina_sv_calling/svaba/svaba_svs.vcf
# DEPENDENCY : sv_genotyping/illumina_svs/svmerged.clustered.vcf
asmvar=../illumina_sv_calling/asmvar/asmvar_filtering/asmvar_svs.vcf
manta=../illumina_sv_calling/manta/manta_svs.vcf
smoove=../illumina_sv_calling/smoove/smoove_svs.vcf
svaba=../illumina_sv_calling/svaba/svaba_svs.vcf
all=../sv_genotyping/illumina_svs/svmerged.clustered.vcf

# Printing the header of the file
printf "program smalldel mediumdel largedel hugedel smallins mediumins largeins hugeins smalldup mediumdup largedup hugedup smallinv mediuminv largeinv hugeinv\n" > table_1_data.txt

# DEPENDENCY : genotyped_lines.txt
for i in $asmvar $manta $smoove $svaba $all
do
	# DEPENDENCY : scripts/count_svtypes_svsizes.awk
	printf "$i " >> table_1_data.txt
	../scripts/count_svtypes_svsizes.awk $i >> table_1_data.txt
done

# Then reprocessing the file to change the values of the first column
gawk 'BEGIN {OFS = ","} 
	/program/ {$1 = $1}
	/asmvar/ {$1 = "asmvar"} 
	/manta/ {$1 = "manta"} 
	/smoove/ {$1 = "smoove"} 
	/svaba/ {$1 = "svaba"} 
	/svmerged/ {$1 = "merged\\tnote{f}"} 
	{print}' table_1_data.txt > table_1.csv

